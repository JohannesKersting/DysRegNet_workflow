suppressMessages(library(DESeq2))
suppressMessages(library(data.table))

# input file is the first cmd argument
args = commandArgs(trailingOnly=TRUE)
input_file_path <- args[1]
meta_file_path <- args[2]
output_file_path <- args[3]

print("Session info:")
print(sessionInfo())

print("Arguments:")
print(args)
# read input
counts <- as.matrix(fread(input_file_path),rownames=1)
meta <- fread(meta_file_path)

# create DESeq2 object
col_data <- data.frame(sample=c(rep("sample",times = ncol(counts))))
dds <- DESeq2::DESeqDataSetFromMatrix(counts,col_data, ~ 1)

# normalize
dds <- DESeq2::estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)

# log2(x+1) transform
normalized_counts <- log2(normalized_counts+1)

# standardize based on control samples
means <- apply(normalized_counts[,meta$condition==0],1,mean)
sds <- apply(normalized_counts[,meta$condition==0],1,sd)
normalized_counts <- apply(normalized_counts,2,function(col){ (col-means)/sds })

# drop rows with 0 sd (produced nans and infs)
if(any(sds==0)){
  print(paste("Dropping ", sum(sds==0), "genes with sd 0 in the control samples..."))
  normalized_counts <- normalized_counts[-which(sds==0),]
}

# turn into data.table
normalized_counts <- as.data.table(normalized_counts, keep.rownames="sample")

# write output
fwrite(normalized_counts, output_file_path)