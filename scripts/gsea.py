import pandas as pd
import gseapy as gp
import sys

# Read command line arguments
input_file = sys.argv[1]
metadata_file = sys.argv[2]
condition_column = sys.argv[3]
gene_sets = sys.argv[4]
output_dir = sys.argv[5]

# Load the normalized counts data
data = pd.read_csv(input_file, index_col=0)

# Load metadata
metadata = pd.read_csv(metadata_file, index_col=0)

# Ensure that the metadata matches the order of samples in the counts data
metadata = metadata.reindex(data.columns)

# Extract conditions (class vector)
conditions = metadata[condition_column].values

# Run GSEA
gp.gsea(
    data=data,
    gene_sets=gene_sets,
    cls=conditions,
    outdir=output_dir,
)
