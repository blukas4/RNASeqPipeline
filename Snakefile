# Load configuration
configfile: "config.yaml"

rule all:
    input:
        f"{config['output_dir']}/normalized_matrix.csv",
        f"{config['output_dir']}/pca_plot.png",
        f"{config['output_dir']}/differential_expression_results.csv",
        f"{config['output_dir']}/volcano_plot.png",
        f"{config['output_dir']}/gsea_results/.done",
        f"{config['output_dir']}/wgcna_modules.csv",
        f"{config['output_dir']}/heatmap.png"

rule normalize:
    input:
        counts = config["raw_counts_file"]
    output:
        f"{config['output_dir']}/normalized_matrix.csv"
    shell:
        "Rscript scripts/normalization.R {input.counts} {output}"

rule pca:
    input:
        counts = f"{config['output_dir']}/normalized_matrix.csv",
        metadata = config["metadata_file"]
    output:
        f"{config['output_dir']}/pca_plot.png"
    params:
        condition_column = config["condition_column"]
    shell:
        "python scripts/pca.py {input.counts} {input.metadata} {params.condition_column} {output}"

rule differential_expression:
    input:
        counts = config["raw_counts_file"],
        metadata = config["metadata_file"]
    output:
        f"{config['output_dir']}/differential_expression_results.csv"
    params:
        condition_column = config["condition_column"]
    shell:
        "Rscript scripts/differential_expression.R {input.counts} {input.metadata} {params.condition_column} {output}"

rule volcano_plot:
    input:
        differential_expression_results = f"{config['output_dir']}/differential_expression_results.csv"
    output:
        f"{config['output_dir']}/volcano_plot.png"
    shell:
        "python scripts/volcano_plot.py {input.differential_expression_results} {output}"

rule gsea:
    input:
        counts = f"{config['output_dir']}/normalized_matrix.csv",
        metadata = config["metadata_file"]
    output:
        sentinel = f"{config['output_dir']}/gsea_results/.done"
    params:
        condition_column = config["condition_column"],
        gene_sets = ",".join(config["gene_sets"]),
        output_dir = f"{config['output_dir']}"
    shell:
        """
        python scripts/gsea.py "{input.counts}" "{input.metadata}" "{params.condition_column}" "{params.gene_sets}" "{params.output_dir}/gsea_results"
        touch {output.sentinel}
        """

rule wgcna:
    input:
        counts = f"{config['output_dir']}/normalized_matrix.csv"
    output:
        f"{config['output_dir']}/wgcna_modules.csv"
    shell:
        "Rscript scripts/wgcna.R {input.counts} {output}"

rule heatmap:
    input:
        counts = f"{config['output_dir']}/normalized_matrix.csv",
        metadata = config["metadata_file"],
        modules = f"{config['output_dir']}/wgcna_modules.csv"
    params:
        condition_column = config["condition_column"]
    output:
        f"{config['output_dir']}/heatmap.png"
    shell:
        "python scripts/heatmap.py {input.counts} {input.metadata} {params.condition_column} {input.modules} {output}"
