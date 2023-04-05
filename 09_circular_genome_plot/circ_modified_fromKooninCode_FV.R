# install Biostrings
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("Biostrings")

library(circlize)
library(tidyverse)
library(GenomicRanges)
library(ggplotify)
library(patchwork)
library(Biostrings)
library(ape)
library(readr)

my_data <- read_csv("circ_genome_df_filtered_16contigs.csv")

plotter <- function(merged, id, strains_name = NULL) {
  ambi <- merged %>% filter(seq_id == id)
  ambi$seq <- gsub("\n", "", ambi$seq)
  ambi$structure_plus <- gsub("\n", "", ambi$structure_plus)
  if(nchar(ambi$structure_plus) != nchar(ambi$seq)){
    print('Error in differet sequence length')
  }
  seq <- RNAStringSet(ambi$seq)
  names(seq) <- "vd"
  #orfs <- findORFsFasta(DNAStringSet(seq), is.circular = TRUE, minimumLength = 300)
  #seqlengths(orfs) <- length(seq[[1]])
  #orfs <- data.frame(id = seqnames(orfs), start = start(orfs), end = end(orfs), strand = strand(orfs))
  orfs <- data.frame('id' = c('vd','vd'), 'start' = c(ambi$orf_plus_from, ambi$orf_minus_to), 
                     'end' = c(ambi$orf_plus_to, ambi$orf_minus_from), 'strand' = c('+','-'))
  split_structure <- strsplit(ambi$structure_plus[1], "")[[1]]
  pairs <- data.frame(to = double(), from = double())
  
  open <- c()
  for (i in 1:length(split_structure)) {
    if (split_structure[i] == "(") {
      open <- c(open, i)
    } else if (split_structure[i] == ")") {
      pairs <- add_row(pairs, from = tail(open, 1), to = i)
      open <- open[1:length(open) - 1]
    }
  }
  circos.initializeCircularGenome("vd", genome_size = length(seq[[1]]))
  circos.track(
    ylim = c(0, 1),
    bg.border = "white",
    panel.fun = function(region, value, ...) {
      if (nrow(orfs) > 0) {
        for (i in 1:nrow(orfs)) {
          circos.arrow(
            x1 = orfs[i, "start"],
            x2 = orfs[i, "end"],
            arrow.position = if (orfs[i, "strand"] == "-") "start" else "end",
            col = if (orfs[i, "strand"] == "-") "#3333FF" else "#33FF33",
          )
        }
      }
      # check if rz_plus is NA
      if (!is.na(ambi$rz_plus[1])){
        circos.rect(
          xleft = ambi$rz_plus_from[1],
          xright = ambi$rz_plus_to[1],
          ybottom = 0,
          ytop = 1,
          col = "#FFCCCC"
        )
      }
      # check if rz_minus is NA
      if (!is.na(ambi$rz_minus[1])){
        circos.rect(
          xleft = ambi$rz_minus_to[1],
          xright = ambi$rz_minus_from[1],
          ybottom = 0,
          ytop = 1,
          col = "#FFCCCC"
        )
      }
      # rdrp
      circos.rect(
        xleft = ambi$rdrp_from,
        xright = ambi$rdrp_to,
        ybottom = 0,
        ytop = 1,
        col = "#66666699"
      )
      if(!is.null(strains_name)){
        title(strains_name)
      }
      for (i in seq(1, nrow(pairs), 1)) {
        circos.link("vd", pairs[i, 1], "vd", pairs[i, 2], lwd = 0.1, h.ratio = 0.9)
      }
    }
  )
}


id_indices_to_run <- c(1:16)
#id_indices_to_run <- 1:30	# adjust as needed
outpath <- "/cloud/project/4k_5k_filtered/"		# adjust as needed

for ( i in id_indices_to_run ){
  
  id_to_visualize <- my_data$seq_id[i]
  
  tiff(  paste0(outpath,id_to_visualize,".tiff"), units="in", width=5, height=5, res=300)
  plotter(my_data, id_to_visualize, id_to_visualize) 
  dev.off()
  
}
