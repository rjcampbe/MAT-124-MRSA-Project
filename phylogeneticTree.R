library(seqinr)
library(adegenet)
library(ape)
library(ggtree)
library(DECIPHER)
library(viridis)
library(ggplot2)

# setwd("I:/My Drive/Spring 2022/MAT 124/midterm2/MAT-124-MRSA-Project")

# load sequences
seqs <- readDNAStringSet("pnas.1702472114.sd04.txt", format = "fasta")

#look at some of the sequences
seqs

#view (aligned) sequences in a browser
BrowseSeqs(seqs, highlight = 0)

# perform alignment
# (ours are already aligned? but not of class alignment,
# wich is required as input for dist.alignment() function)
# i guess if they are already aligned nothing will change?
seqs_aligned <- AlignSeqs(seqs)

# compare to previous sequences in browser
# (looks the same)
BrowseSeqs(seqs_aligned, highlight = 0)

# write the alignment to a new FASTA file
# (maybe could have started here?)
writeXStringSet(seqs_aligned, file = "mrsa_aligned.fasta")

# read in the aligned data
mrsa <- read.alignment("mrsa_aligned.fasta", format = "fasta")

# create a distance matrix for the alignment
D <- dist.alignment(mrsa, matrix = "similarity")
dist_df <- as.data.frame(as.matrix(D))

# darker shades of gray mean a larger distance 
# you can also make cool color plots 
# but they're much more complicated because they use the image() function
table.paint(dist_df, cleg=0, clabel.row=.5, clabel.col=.5)+
  scale_color_viridis()

# this uses the neighbor joining method (not MLE)
# i.e. Saitou and Nei (1987)
tre_nj <- ape::nj(D)

# all trees created using {ape} package will be of class phylo
class(tre_nj) 

# This function reorganizes the internal structure of the tree 
# to get the ladderized effect when plotted.
tre_nj <- ladderize(tre_nj) 

# ~~~~~~~~~~~~~~~~~~~~~~ Base R plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(tre_nj)

# method = average is used for UPGMA, 
# members can be equal to NULL or a vector with a length of size D
h_cluster <- hclust(D, method = "average", members = NULL) 
plot(h_cluster, cex = 0.6)
# =======================================================================

# ~~~~~~~~~~~~~~~~~~~~~~ Using ggtree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# circular layout (hard to read)
ggtree(tre_nj, layout = "circular")+
  geom_tiplab(size = 3)

# circular unrooted layout (easier to read)
ggtree(tre_nj,branch.length = "none", layout = "circular")+
  geom_tiplab(size = 3)

# Cladogram: rectangular layout
ggtree(tre_nj, branch.length = "none")+
  geom_tiplab(size = 3)

# top-down rectangular layout (Dendrogram)
ggtree(tre_nj)+ 
  layout_dendrogram()
















