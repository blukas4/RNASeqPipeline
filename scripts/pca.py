import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import sys

# Read command line arguments
input_file = sys.argv[1]
metadata_file = sys.argv[2]
condition_column = sys.argv[3]
output_file = sys.argv[4]

# Load the normalized counts data
data = pd.read_csv(input_file, index_col=0)

# Load metadata
metadata = pd.read_csv(metadata_file, index_col=0)

# Ensure that the metadata matches the order of samples in the counts data
metadata = metadata.reindex(data.columns)

# Extract conditions
conditions = metadata[condition_column].values

# Perform PCA
pca = PCA(
    n_components=2
)  # We only need the first two principal components for our plot
principal_components = pca.fit_transform(
    data.T
)  # Transpose since PCA operates on samples as rows
pc_df = pd.DataFrame(data=principal_components, columns=["PC1", "PC2"])

# Plotting
plt.figure(figsize=(5, 5), dpi=300)

# Create a scatter plot for each unique condition (sorted alphabetically)
for condition in sorted(set(conditions)):
    subset = pc_df[conditions == condition]
    plt.scatter(subset["PC1"], subset["PC2"], s=50, label=condition)

plt.title("PCA of RNA-seq Data")
plt.xlabel(f"PC1: {pca.explained_variance_ratio_[0]*100:.2f}%")
plt.ylabel(f"PC2: {pca.explained_variance_ratio_[1]*100:.2f}%")
plt.legend()
plt.tight_layout()

# Save plot
plt.savefig(output_file)
