library(DESeq2)

# Get arguments
args <- commandArgs(trailingOnly = TRUE)

raw_counts_file <- args[1]
metadata_file <- args[2]
condition_column <- args[3]
output_file <- args[4]

# Read in raw count data
raw_counts <- read.csv(raw_counts_file, row.names = 1, check.names = FALSE)

# Filter out genes with low counts across all samples
keep_genes <- rowSums(raw_counts) >= 100
raw_counts <- raw_counts[keep_genes, ]

# Read metadata
metadata <- read.csv(metadata_file, row.names = 1)

# Ensure the metadata matches the order of samples in the counts data
metadata <- metadata[colnames(raw_counts), ]

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(
  countData = raw_counts,
  colData = metadata,
  design = formula(paste("~", condition_column))
)

# Differential expression analysis
dds <- DESeq(dds)

# Getting results
res <- results(dds)

# Save results to CSV
write.csv(as.data.frame(res), file = output_file)
