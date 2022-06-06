#                       Header material (setup)
# =====================================================================
# Package Installation:
# Use the commands:
# install.packages("BiocManager")
library(BiocManager)
# BiocManager::install("Biostrings")
library(Biostrings)
# BiocManager::install("RSQLite")
library(RSQLite)
# BiocManager::install("DECIPHER")
library(DECIPHER)


# Packages used
library(seqinr)
library(adegenet)
library(ape)
library(ggtree) # Warning message received: 'ggtree' was built under R version 4.0.5. Manually added in packages tab.
library(DECIPHER)
library(RSQLite)
library(viridis)
library(ggplot2)

# TDA package
# install.packages("TDAstats")
library(TDAstats)

# To Check Working Directory:
getwd()

# To Set your Working Directory:
# Riley's:
# setwd("I:/My Drive/Spring 2022/MAT 124/midterm2/MAT-124-MRSA-Project")
# 
# Ryan's:
# setwd("/Users/ryancampbell/Documents/GitHub/MAT-124-MRSA-Project")
# 
# Aditya's:
# setwd("C:/Users/adity/Documents/GitHub/MAT-124-MRSA-Project")
# 
# A Relative path for a cloned repo:
# setwd("~/MAT-124-MRSA-Project")

#                       Beginning of Code
# =====================================================================
# load sequences
seqs <- readDNAStringSet("pnas.1702472114.sd04.txt", format = "fasta")

#look at some of the sequences
seqs

#view (aligned) sequences in a browser
# BrowseSeqs(seqs, highlight = 0)

# perform alignment
# (ours are already aligned? but not of class alignment,
# wich is required as input for dist.alignment() function)
# i guess if they are already aligned nothing will change?
seqs_aligned <- AlignSeqs(seqs)

# compare to previous sequences in browser
# (looks the same)
# BrowseSeqs(seqs_aligned, highlight = 0)

# write the alignment to a new FASTA file
# (maybe could have started here?)
writeXStringSet(seqs_aligned, file = "mrsa_aligned.fasta")

#                       DISTANCE MATRIX + TDA
# =====================================================================
# read in the aligned data
mrsa <- read.alignment("mrsa_aligned.fasta", format = "fasta")
mrsa

# create a distance matrix for the alignment
D <- dist.alignment(mrsa, matrix = "similarity")
dist_df <- as.matrix(D)

# copy distance matrix into a matrix variable num
# num <- matrix(unlist(dist_df), ncol = 224, nrow = 224) 
# I don't think this is necessary... Seems to be 1-1 with dist_df

# The num variable now contains the distance matrix in R. 
# To find the homologies, we use the following commands: 
# install.packages("TDAstats")
library("TDAstats")
p <- calculate_homology(dist_df, dim =1, threshold = -1, format = "distmat", standardize = FALSE, return_df = FALSE)
plot_barcode(p)
plot_persist(p)
