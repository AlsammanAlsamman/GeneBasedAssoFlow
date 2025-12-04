"""
Snakemake rule: add_gene_names
Add gene names/symbols to MAGMA gene-based analysis results
Maps Entrez Gene IDs to gene symbols using NCBI gene info file
"""

import sys
sys.path.append("utils")
from bioconfigme import get_results_dir, get_analysis_value

# Initialize configs
from bioconfigme import init_configs
init_configs("configs/software.yml", "configs/analysis.yml", "configs/config.yml")

# Get results directory
RESULTS_DIR = get_results_dir()

rule add_gene_names:
    input:
        genes_out = RESULTS_DIR + "/{target}/magma/gene_analysis/{target}.genes.out",
        done_marker = RESULTS_DIR + "/{target}/magma/gene_analysis/{target}_gene_analysis.done"
    output:
        genes_with_names = RESULTS_DIR + "/{target}/magma/gene_analysis/{target}_genes_with_names.txt",
        gene_mapping = RESULTS_DIR + "/{target}/magma/gene_analysis/gene_id_to_name.txt",
        done = RESULTS_DIR + "/{target}/magma/gene_analysis/{target}_add_names.done"
    params:
        gene_info = get_analysis_value('reference.ncbi_gene_info'),
        script = "scripts/add_gene_names.sh"
    log:
        RESULTS_DIR + "/{target}/magma/log/add_gene_names_{target}.log"
    resources:
        mem_mb = 32000,
        time = "00:30:00"
    threads: 2
    shell:
        """
        # Create output directories
        mkdir -p $(dirname {output.genes_with_names})
        mkdir -p $(dirname {log})
        
        # Run gene name addition script
        bash {params.script} \
            --genes-out {input.genes_out} \
            --gene-info {params.gene_info} \
            --output {output.genes_with_names} \
            --mapping-output {output.gene_mapping} \
            > {log} 2>&1
        
        # Create done marker
        touch {output.done}
        echo "Added gene names to results for {wildcards.target}" >> {output.done}
        """
