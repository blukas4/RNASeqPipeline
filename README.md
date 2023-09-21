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

## Requirements

**Snakemake**: This project is a Snakemake pipeline, and you need to have Snakemake installed to run the workflows. To install Snakemake, you can use `conda`:
  ```bash
  conda install snakemake
  ```

**R**: Some steps in this pipeline use R scripts. Ensure you have [R installed](https://cran.r-project.org/mirrors.html). 

  R packages required:
  - `WGCNA`
  - `DESeq2`

  To install these R packages, open an R session and run:

  ```R
  install.packages(c("BiocManager"))
  BiocManager::install(c("WGCNA", "DESeq2"))
  ```

This project has several dependencies which need to be installed for proper functionality. Here's a list of the required libraries:

- `sklearn`
- `gseapy`
- `matplotlib`
- `numpy`
- `pandas`
- `seaborn`

To install these requirements, you can use `pip`:

```bash
pip install -r requirements.txt
```


## Quick Start

1. Clone the repository:
   ```bash
   git clone https://github.com/blukas4/RNASeqPipeline.git
   ```

2. Navigate to the directory:
   ```bash
   cd RNASeqPipeline
   ```

3. Configure your analysis by editing the `config.yaml` file.

4. Run the pipeline:
   ```bash
   snakemake --cores all
   ```

## Configuration

Edit the `config.yaml` file to set up your analysis.

## Contribute

Feel free to fork, open pull requests, or submit issues. We appreciate any feedback or contributions!
