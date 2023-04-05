# packages
library( tidyverse )
library( cowplot )
library( phytools )
library( ggtree )

# input
treefile  <- "Ambi_RdRp_phyml_SH_tree.nwk"
taxfile   <- "Fungi_SRA_lineage_table.tsv"
minPAT    <- 0.0
ColAtRank <- "order"

# set seed
set.seed(4711)

# condense tree
tree0 <- read.tree( treefile )
system( paste("Rscript tree_condense_by_distance.R", treefile, "tree_condensed.nwk", minPAT), intern=F )

# load tree
tree <- read.tree( "tree_condensed.nwk" )
write( tree$tip.label, file="tree_condensed.ids" )
tree0$tip.label[ ! tree0$tip.label %in% tree$tip.label ]
print( paste( length(tree0$tip.label)-length(tree$tip.label), "tips removed" )  )

# midpoint-pseudoroot tree
tree <- midpoint.root(tree)

# read taxonomy information
tax <- read.delim( taxfile, sep="\t", header=F )
colnames(tax) <- c("accession","taxon","rank")

# color tree branches by host
hosts <- tree %>% as_tibble() %>% mutate( scientific_name = "metatranscriptome" )
for ( id in tax %>% filter(rank==ColAtRank) %>% pull(accession) ){
	nam   <- tax %>% filter( accession==id, rank==ColAtRank ) %>% pull(taxon)
	nam   <- if_else( nam=="", "unclassified", nam )
	hosts <- hosts %>% mutate( scientific_name = if_else( grepl(id,label), nam, scientific_name ) )
}
head(hosts)
table(hosts$scientific_name)

grouping <- list()
for ( hh in unique(hosts$scientific_name) ){
	grouping[[ hh ]] <- hosts %>% filter(scientific_name==hh) %>% pull(label)
}

tree <- groupOTU( tree, grouping )

# define coloring
scinames <- hosts %>% pull(scientific_name) %>% sort() %>% unique()
col_groups <- rainbow( length(scinames) )
col_groups <- sample( col_groups, length(col_groups) )
names(col_groups) <- scinames
col_groups["metatranscriptome"] <- "gray80"

# which nodes to highlight
nodes_lab_num   <- as.numeric( tree$node.label )
nodes_supported <- which( nodes_lab_num >  0.9 )
nodes_supported <- nodes_supported + length(tree$tip.label)

# produce tree plot
p1 <- ggtree(tree, aes(color=group), layout="circular" ) 
p1 <- p1 	+ guides( color=guide_legend(title=ColAtRank) )  
p1 <- p1	+ scale_color_manual( values=col_groups )
p1 <- p1	+ geom_point2( aes(subset=(node %in% nodes_supported)), color=rgb(0,0,0,0.5), size=0.5 )
p1 <- p1	+ geom_treescale( x=1, y=690, offset=5 )

p2 <- ggtree(tree, aes(color=group), layout="circular", branch.length='none' ) 
p2 <- p2 	+ guides( color=guide_legend(title=ColAtRank) )  
p2 <- p2	+ scale_color_manual( values=col_groups )
p2 <- p2	+ geom_treescale()

# save tree plot
ggsave( p1, file=paste0("tree_condensed_circular_",ColAtRank,"_phylo.pdf"),  limitsize=F, width=7, height=6 )
ggsave( p2, file=paste0("tree_condensed_circular_",ColAtRank,"_dendro.pdf"), limitsize=F, width=7, height=6 )

