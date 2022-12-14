suppressMessages(library(DESeq2))
suppressMessages(library(data.table))

# input file is the first cmd argument
args = commandArgs(trailingOnly=TRUE)
input_file_path <- args[1]
output_file_path <- args[2]

print("Session info:")
print(sessionInfo())

print("Arguments:")
print(args)
# read input
counts <- as.matrix(fread(input_file_path),rownames=1)

# create DESeq2 object
col_data <- data.frame(sample=c(rep("sample",times = ncol(counts))))
dds <- DESeq2::DESeqDataSetFromMatrix(counts,col_data, ~ 1)

# normalize and turn into data.table
dds <- DESeq2::estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts <- as.data.table(normalized_counts, keep.rownames="sample")

# write output
fwrite(normalized_counts, output_file_path)