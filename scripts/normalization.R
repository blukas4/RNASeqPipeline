library(DESeq2)

# Get arguments
args <- commandArgs(trailingOnly = TRUE)

raw_counts_file <- args[1]
output_file <- args[2]

# Read in raw count data
raw_counts <- read.csv(raw_counts_file, row.names = 1, check.names = FALSE)

# Filter out genes with low counts across all samples
keep_genes <- rowSums(raw_counts) >= 100
raw_counts <- raw_counts[keep_genes, ]

# Create dummy metadata
sample_names <- colnames(raw_counts)
dummy_metadata <- data.frame(row.names = sample_names)

# Create DESeq2 object with a null model
dds <- DESeqDataSetFromMatrix(
  countData = raw_counts,
  colData = dummy_metadata,
  design = ~NULL
)

# Perform normalization
dds <- DESeq(dds)

# Extract normalized counts
normalized_counts <- counts(dds, normalized = TRUE)

# Write normalized counts to CSV
write.csv(normalized_counts, file = output_file)
