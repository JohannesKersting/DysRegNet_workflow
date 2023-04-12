suppressPackageStartupMessages(library(lionessR))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(arrow))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()

parser$add_argument("-e", "--expr",
                    help="Path to a .csv file containing the counts")
parser$add_argument("-m", "--meta",
                    help="Path to a file containing patient meta data")
parser$add_argument("-g", "--grn",
                    help="Path to a reference network")
parser$add_argument("-o", "--output",
                    help="Path to the output file")

args <- parser$parse_args()

print("Session info:")
print(sessionInfo())

print("args:")
print(args)


print("read data and correct format")
expression_data1 <- as.data.frame(fread(args$expr))
expression_data <- expression_data1[,-1]
rownames(expression_data) <- expression_data1[,1]
expression_data1 <- NULL
meta_data <- as.data.frame(fread(args$meta))
reference_network <- as.data.frame(fread(args$grn))

print("make sure the genes in the reference file are in the expression data")
names(reference_network)[1] <- "SYMBOL_TF" 
names(reference_network)[2] <- "SYMBOL_TG" 
reference_network <- reference_network %>% filter (SYMBOL_TF %in% rownames(expression_data) & SYMBOL_TG %in% rownames(expression_data))

print("only use those genes that are in the reference network")
expression_data_filtered <- expression_data %>% filter(rownames(expression_data) %in% reference_network$SYMBOL_TF | rownames(expression_data) %in% reference_network$SYMBOL_TG)
expression_data_filtered <- expression_data_filtered %>% select(meta_data$sample)

print("run lioness")
se <- SummarizedExperiment(assays = list(counts = as.matrix(expression_data_filtered)), 
                           colData = meta_data, rowData = expression_data_filtered)
lioness_result <- lioness(se, netFun)
print(lioness_result)
result <- assay(lioness_result)
lioness_result <- NULL


print("filter for reference network")
mergedNetwork <- reference_network %>% unite("edges",SYMBOL_TF:SYMBOL_TG, sep = "_", remove = TRUE)
filtered_result <- data.matrix(as.data.frame(result) %>% filter(rownames(result) %in% mergedNetwork$edges))
print(str(filtered_result))
result <- NULL

print("calculate Z scores and set unsignificant edges 0")

result_case <- filtered_result[ , meta_data$condition == 1]
result_control <- filtered_result[ , meta_data$condition == 0]

res <- data.frame(matrix( nrow = ncol(result_case)), row.names = colnames(as.data.frame(result_case)))
res[, rownames(result_case)] <- t((result_case-rowMeans(result_control))/ rowSds(result_control)) 
res <- as.matrix(res[,-1])
res[2*pnorm(q= abs(res), lower.tail=FALSE) * ncol(result_case) > 0.01] <- 0

#correct format of colnames
colnames(res) <- paste0('("', str_replace(colnames(res), "_", '", "'), '")')

print(str(res))


res<- as.data.table(res, keep.rownames = "patient id")
arrow::write_feather(res, args$output)
