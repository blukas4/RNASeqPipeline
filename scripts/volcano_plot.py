import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

# Read command line arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

# Load differential expression results
results = pd.read_csv(input_file, index_col=0)

# Calculate -log10(adjusted p-value)
results["-log10(adjusted pvalue)"] = -np.log10(results["padj"])

# Defining the significance criteria
alpha = 0.05
log2fc_threshold = 1

# Categorize genes based on significance and fold change
# Not Statistically Significant (NSS)
nss = results["padj"] >= alpha
# Statistically Significant, but Fold Change not exceeding threshold (SS-FC)
ss_fc = (results["padj"] < alpha) & (abs(results["log2FoldChange"]) <= log2fc_threshold)
# Statistically Significant and Fold Change exceeds threshold (SS+FC)
ss_plus_fc = (results["padj"] < alpha) & (
    abs(results["log2FoldChange"]) > log2fc_threshold
)

# Plotting
plt.figure(figsize=(5, 5), dpi=300)

# Plotting each category
plt.scatter(
    results["log2FoldChange"][nss],
    results["-log10(adjusted pvalue)"][nss],
    color="gray",
    s=5,
    label="NSS",
)
plt.scatter(
    results["log2FoldChange"][ss_fc],
    results["-log10(adjusted pvalue)"][ss_fc],
    color="blue",
    s=5,
    label="SS-FC",
)
plt.scatter(
    results["log2FoldChange"][ss_plus_fc],
    results["-log10(adjusted pvalue)"][ss_plus_fc],
    color="red",
    s=5,
    label="SS+FC",
)

# Labelling axes and title
plt.title("Volcano plot")
plt.xlabel("log2(Fold Change)")
plt.ylabel("-log10(Adjusted p-value)")

# Adding vertical lines to denote the log2 fold change thresholds
plt.axvline(x=log2fc_threshold, linestyle="--", color="black")
plt.axvline(x=-log2fc_threshold, linestyle="--", color="black")

# Adding horizontal line to denote the alpha threshold
plt.axhline(y=-np.log10(alpha), linestyle="--", color="black")

plt.legend()
plt.tight_layout()

# Save the plot
plt.savefig(output_file)
