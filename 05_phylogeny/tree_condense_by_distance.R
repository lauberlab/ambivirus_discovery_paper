# packages
library(ape)


# input
args = commandArgs( trailingOnly=T )
if (length(args) != 3 ){
	cat("\nusage: tree_condense_by_distance.R <tree_in.nwk> <tree_out.nwk> <min_distance>\n\n")
	q()
}
IFILE = as.character( args[1] )
OFILE = as.character( args[2] )
MIND  = as.numeric(   args[3] )


# read tree
tree <- read.tree( file=IFILE )
#tree


# look at distribution of distances, which may help to define the distance threshold used for removal
dm <- cophenetic( tree ); dm <- dm[upper.tri(dm)]
pdf( "tree_condensed.pdf" )
hist( dm,          breaks=50, main="all"   )
hist( dm[dm<1],    breaks=50, main="<1"    )
hist( dm[dm<0.1],  breaks=50, main="<0.1"  )
hist( dm[dm<0.01], breaks=50, main="<0.01" )
dev.off()


# iteratively remove one sequence for closest pair unless all pairs have minimum or higher required distance
tree2   <- tree
removed <- character()
while ( T ){
	# calculate PATs
	pat <- cophenetic( tree2 )
	pat[ lower.tri(pat) ] <- NA
	diag(pat) <- NA

	# get row, col and value of smallest distance
	patmin <- which(pat == min(pat,na.rm=T), arr.ind=TRUE)[1,]
	ri     <- as.integer( patmin[1] )
	ci     <- as.integer( patmin[2] )
	di     <- pat[ri,ci]
	#cat( paste( ri, ci, di,"\n", sep="," ) )

	# stop if all PATs large enough
	if( is.na(di) || di >= MIND ){
		break
	}

	# remove sequence
	#cat( paste("remove sequence", rownames(pat)[ri]),"\n" )
	seqrm   <- rownames(pat)[ri]
	seqj    <- colnames(pat)[ci]
	if( grepl("_DRR",seqj) | grepl("_ERR",seqj) | grepl("_SRR",seqj) ){
		seqrm <- seqj
	}
	tree2   <- drop.tip( tree2, seqrm )
	removed <- c( removed, seqrm )
}
#tree2


# write resulting tree to file
write.tree( tree2, file=OFILE )


# report removed sequence names
removed

