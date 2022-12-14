#!/usr/bin/env Rscript
library(GENIE3)
library(doParallel)
library(doRNG)
library(data.table)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("-e", "--expr")
parser$add_argument("-m", "--meta")
parser$add_argument("--known_tfs")
parser$add_argument("-t", "--threads", type="integer")
parser$add_argument("-o", "--output")


args <- parser$parse_args()

input_file_path <- args$expr
meta_file_path <- args$meta
tf_file_path <- args$known_tfs
n_cores <- args$threads
out_file_path <- args$output

print("Session info:")
print(sessionInfo())

print("Arguments:")
print(args)

meta <- fread(meta_file_path)
exprMatr <- as.matrix(fread(input_file_path),rownames=1)
exprMatr <- exprMatr[,meta$condition==0]

print("Matrix format:")
print(dim(exprMatr))


set.seed(123) # For reproducibility of results

reg <-  read.delim(tf_file_path, header = FALSE)
print("Number of known TFs:")
print(length(reg$V1))

print("Number of known TFs in expr data:")
print(length(intersect(rownames(exprMatr), reg$V1)))

weightMat <- GENIE3(exprMatr, nCores=n_cores, verbose = TRUE, regulators = intersect(rownames(exprMatr), reg$V1))
linkList <- getLinkList(weightMat)
print(head(linkList))

write.csv(linkList, out_file_path, row.names = FALSE, quote=F)
