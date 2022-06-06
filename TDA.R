#                       Header material (setup)
# =====================================================================
# Package Installation:
# Use the commands:
# install.packages("BiocManager")
# BiocManager::install("Biostrings")
# BiocManager::install("RSQLite")
# BiocManager::install("DECIPHER")
# install.packages("TDAstats")
library(BiocManager)
library(Biostrings)
library(RSQLite)


# Packages used
library(seqinr)
library(adegenet)
library(ape)
library(ggtree)
library(DECIPHER)
library(RSQLite)
library(viridis)
library(ggplot2)
library(BiocManager)
library(Biostrings)
library(RSQLite)
library(TDAstats)

# TDA package
# install.packages("TDAstats")

#                       Beginning of Code (setup)
# =====================================================================
# load sequences
seqs <- readDNAStringSet("pnas.1702472114.sd04.txt", format = "fasta")

#look at some of the sequences
seqs

#view (aligned) sequences in a browser (uncomment to browse)
# BrowseSeqs(seqs, highlight = 0)

# perform alignment (if necessary)
seqs_aligned <- AlignSeqs(seqs)

# compare to previous sequences in browser (uncomment to browse)
# BrowseSeqs(seqs_aligned, highlight = 0)

# write the alignment to a new FASTA file
writeXStringSet(seqs_aligned, file = "mrsa_aligned.fasta")

#                       DISTANCE MATRIX + TDA
# =====================================================================
# read in the aligned data
mrsa <- read.alignment("mrsa_aligned.fasta", format = "fasta")
mrsa

# create a distance matrix for the alignment
D <- dist.alignment(mrsa, matrix = "similarity")
dist_df <- as.matrix(D)

# The dist_df variable now contains the distance matrix in R. 
# To find the homologies, we use the following commands: 
p <- calculate_homology(dist_df, dim =1, threshold = -1, format = "distmat", standardize = FALSE, return_df = FALSE)

# Plot the TDA graphs
plot_barcode(p)
plot_persist(p)
