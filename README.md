# RNASeqPipeline

A comprehensive RNA-seq data analysis pipeline that normalizes raw count matrices, performs PCA, differential expression analysis, visualizes results with volcano plots, identifies biological pathways with GSEA, and clusters genes using WGCNA.

## Features

- **Normalization**: Convert raw counts matrix to normalized counts.
- **PCA**: Use Principal Component Analysis to check if samples can be differentiated based on conditions.
- **Differential Expression Analysis**: Identify genes that are differentially expressed between conditions.
- **Volcano Plot**: Visualize the differentially expressed genes.
- **GSEA**: Identify which biological pathways are activated or repressed.
- **WGCNA**: Cluster genes based on co-expression into modules.
- **Heatmap**: Hierarchical clustering of gene expression with module and condition annotations.

## Quick Start

1. Clone the repository:
   ```bash
   git clone https://github.com/buriedsand/RNASeqPipeline.git
   ```

2. Navigate to the directory:
   ```bash
   cd RNASeqPipeline
   ```

3. Install the required packages using the provided environment file (requires conda):
   ```bash
   conda env create -f envs/environment.yaml
   ```

4. Activate the environment:
   ```bash
   conda activate rnaseq-env
   ```

5. Configure your analysis by editing the `config.yaml` file.

6. Run the pipeline:
   ```bash
   snakemake --use-conda --configfile path/to/yourconfig.yaml
   ```

## Configuration

Edit the `config.yaml` file to set up your analysis.

## Contribute

Feel free to fork, open pull requests, or submit issues. We appreciate any feedback or contributions!