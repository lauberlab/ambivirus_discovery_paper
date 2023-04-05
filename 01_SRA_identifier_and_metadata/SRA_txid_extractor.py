#######################################################################
# Fx:                                                                 #
# Last update: 18 Nov 2022                                             #
# Dependencies:                                                       #
# 1. Biopython - pip install biopython                                #
# 2. pysradb web - pip install git+https://github.com/saketkc/pysradb #
# 3. ete3 - pip install ete3                                          #
#######################################################################

import os
import argparse
import logging
import subprocess
from Bio import SeqIO
import pandas as pd
from pysradb.sraweb import SRAweb
from ete3 import NCBITaxa
import queue
import threading

SRR_list_df = None
df_sraDB = None

#----------#
# Argument #
#----------#
def get_args():
    
   parser = argparse.ArgumentParser(description='Please paste script description here', formatter_class=argparse.HelpFormatter)
   parser.add_argument('-i', '--input', help='input file')
   parser.add_argument('-o', '--output', default= '.', help='output directory')
   parser.add_argument('-t', '--threads_num', default=8, help='number of threads')
   
   return parser.parse_args()

def txt2list(input):
    global SRR_list_df
    
    # read the txt file 
    SRR_list_df = pd.read_csv(input, names=['SRA_id'])

def split_dataframe(df, chunk_size = 10000): 
    chunks = list()
    num_chunks = len(df) // chunk_size + 1
    for i in range(num_chunks):
        chunks.append(df[i*chunk_size:(i+1)*chunk_size])
    return chunks

def search_metadata(ids: queue.Queue):
    global df_sraDB
    while True:
        try:
            SRR_list = ids.get_nowait()
            sra_db = SRAweb()
            alldf_sraDB = []
            a = SRR_list.SRA_id.tolist()
            sra_metadata_matched = sra_db.sra_metadata(a, detailed=True)

            for i, row in sra_metadata_matched.iterrows():
                if row['run_accession'] not in a:
                    sra_metadata_matched.drop(i, inplace=True)
            alldf_sraDB.append(sra_metadata_matched)

            for i in alldf_sraDB:
                df_sraDB = pd.concat([df_sraDB, i])
        except queue.Empty:
            break
        else:
            pass
    return True

class myThread(threading.Thread):
    def __init__(self, threadID, queues):
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.queues = queues
    def run(self):
        search_metadata(self.queues)

def fetch_txid(num_of_threads, output):
    logging.info('Fetching the Taxon ID based on the wanted SRA ID')
    global df_sraDB
    
    SRR_list_split = split_dataframe(SRR_list_df, chunk_size=50)
    srr_queues = queue.Queue()
    
    for i in SRR_list_split:
        srr_queues.put(i)
    
    df_sraDB = pd.DataFrame(data = None, columns= SRR_list_df.columns)
    num_of_threads = int(num_of_threads)
    threads = []

    for i in range(num_of_threads):
        threads_asd = myThread(i, srr_queues)
        threads_asd.start()
        threads.append(threads_asd)
    
    for t in threads:
        t.join()

    df_sraDB = df_sraDB[['run_accession', 'organism_taxid ']]
    df_sraDB = df_sraDB.rename(columns={'run_accession': 'SRA_id', 'organism_taxid ': 'taxid'})
    df_sraDB = df_sraDB.dropna()
    df_sraDB = df_sraDB.reset_index(drop=True)
    outfile_df_sraDB = output + "/df_sraDB.csv"
    df_sraDB.to_csv(outfile_df_sraDB, index = False)

logging.info('Importing NCBI Taxon folder')
ncbi = NCBITaxa()

def get_ranks(taxid, desired_ranks):
    lineage = ncbi.get_lineage(taxid)   
    names = ncbi.get_taxid_translator(lineage)
    lineage2ranks = ncbi.get_rank(names)
    ranks2lineage = dict((rank,taxid) for (taxid, rank) in lineage2ranks.items())
    return{'{}_id'.format(rank): ranks2lineage.get(rank, '<not present>') for rank in desired_ranks}

def get_full_lineages(output):
    logging.info('Fetching full lineages of the respective taxon ID')
    results = list()
    desired_ranks = ['superkingdom', 'clade', 'kingdom', 'subkingdom', 'superphylum', 'phylum', 'subphylum', 'superclass', 'class', 
    'subclass', 'superorder', 'order', 'suborder', 'superfamily', 'family', 'subfamily', 'genus', 'subgenus', 'species']
    wanted_taxid = df_sraDB.taxid.tolist()
    
    for taxid in wanted_taxid:
        results.append(list())
        results[-1].append(str(taxid))
        ranks = get_ranks(taxid, desired_ranks)
        for key, rank in ranks.items():
            if rank != '<not present>': 
                results[-1].append(list(ncbi.get_taxid_translator([rank]).values())[0])
            else: 
                results[-1].append(rank)

    i = 0 
    outfile_t2l = output + "/txid2lineage.tsv"
    txid2lineage = open(outfile_t2l,'a')
    for result in results:
        txid2lineage.write('\t'.join(result) + '\n')
    
        i += 1

    txid2lineage.close()

def process_lineage_table(output):
    logging.info('Group all information into a dataframe')

    tax_data_1 = pd.merge(SRR_list_df, df_sraDB, on='SRA_id', how='left')
    tax_data_1['taxid'] = tax_data_1['taxid'].astype('Int64')
    
    taxonomy_file = output + "/txid2lineage.tsv"
    tax_data = pd.read_csv(taxonomy_file, sep='\t', names=['taxid', 'superkingdom', 'clade', 'kingdom', 'subkingdom', 'superphylum', 'phylum', 'subphylum', 'superclass', 'class', 'subclass', 'superorder', 'order', 'suborder', 'superfamily', 'family', 'subfamily', 'genus', 'subgenus', 'species'])
    tax_data.drop_duplicates(keep='first', inplace=True)
    tax_data['taxid'] = tax_data['taxid'].astype('Int64')
    full_tax_data = pd.merge(tax_data_1, tax_data, on='taxid', how='left')

    outfile_full_tax_data = output + "/lineage_table.csv"
    full_tax_data.to_csv(outfile_full_tax_data, index = False)

if __name__ == '__main__':
    args = get_args()
    
    os.mkdir(args.output)

    logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)
    txt2list(args.input)
    fetch_txid(int(args.threads_num), args.output)
    get_full_lineages(args.output)
    process_lineage_table(args.output)
