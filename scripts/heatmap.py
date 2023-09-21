import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from utils.color_conversion import get_hex
import sys

# Number of top variable genes to retain
N = 50

# Read command line arguments
counts_file = sys.argv[1]
metadata_file = sys.argv[2]
condition_column = sys.argv[3]
modules_file = sys.argv[4]
output_file = sys.argv[5]

# Load the data
data = pd.read_csv(counts_file, index_col=0)
metadata = pd.read_csv(metadata_file, index_col=0)
modules = pd.read_csv(modules_file, index_col=0)

# Ensure the metadata matches the order of samples in the counts data
metadata = metadata.reindex(data.columns)

# Create column colors based on sample conditions
conditions = metadata[condition_column]
unique_conditions = sorted(set(conditions))
palette = sns.color_palette("tab10", len(unique_conditions))
condition_colors = conditions.map(dict(zip(unique_conditions, palette)))

# Calculate variance for each gene
gene_variability = data.var(axis=1)

# Get the top N most variable genes' indices
top_genes = gene_variability.nlargest(N).index

# Subset the data
data = data.loc[top_genes]

# Create row colors based on gene module membership
modules = modules.reindex(data.index)
module_colors1 = modules["Dynamic_Color"].map(get_hex)
module_colors2 = modules["Merged_Color"].map(get_hex)

# Concatenate the colors into a DataFrame
module_colors = pd.DataFrame({"Dynamic": module_colors1, "Merged": module_colors2})

# Create the clustermap
sns.clustermap(
    data,
    method="average",
    metric="euclidean",
    z_score=0,  # Standardize gene expression for better visualization
    col_colors=condition_colors,
    row_colors=module_colors,
    figsize=(10, 10),
    cmap="vlag",
    yticklabels=1,
    xticklabels=1,
)

# Save the clustermap to a file
plt.savefig(output_file, dpi=300, bbox_inches="tight")
