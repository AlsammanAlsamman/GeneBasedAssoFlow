"""
Snakemake rule: create_loci_table
Create loci-based Excel table mapping significant genes to genomic loci
Maps genes with P < threshold to loci regions with distance buffer
"""

import sys
sys.path.append("utils")
from bioconfigme import get_results_dir, get_software_module, get_target_analysis, get_dataset, get_analysis_value

# Initialize configs
from bioconfigme import init_configs
init_configs("configs/software.yml", "configs/analysis.yml", "configs/config.yml")

# Get results directory
RESULTS_DIR = get_results_dir()

rule create_loci_table:
    input:
        genes_with_names = RESULTS_DIR + "/{target}/magma/gene_analysis/{target}_genes_with_names.txt",
        done_marker = RESULTS_DIR + "/{target}/magma/gene_analysis/{target}_add_names.done"
    output:
        loci_table = RESULTS_DIR + "/{target}/magma/loci_analysis/{target}_loci_genes.xlsx",
        done = RESULTS_DIR + "/{target}/magma/loci_analysis/{target}_loci_table.done"
    params:
        gene_loc = get_analysis_value('reference.gene_loc'),
        loci_file = lambda wildcards: get_dataset(get_target_analysis(wildcards.target)['dataset'])['loci_file'],
        p_threshold = lambda wildcards: get_target_analysis(wildcards.target)['loci_gene_p_thereshold'],
        distance_kb = lambda wildcards: get_target_analysis(wildcards.target)['loci_distance_thereshold_kb'],
        r_module = get_software_module('r'),
        script = "scripts/create_loci_table.R"
    log:
        RESULTS_DIR + "/{target}/magma/log/create_loci_table_{target}.log"
    resources:
        mem_mb = 32000,
        time = "00:30:00"
    threads: 2
    shell:
        """
        # Create output directories
        mkdir -p $(dirname {output.loci_table})
        mkdir -p $(dirname {log})
        
        # Load R module
        module load {params.r_module}
        
        # Run loci table creation script
        Rscript {params.script} \
            --genes-file {input.genes_with_names} \
            --gene-loc {params.gene_loc} \
            --loci-file {params.loci_file} \
            --p-threshold {params.p_threshold} \
            --distance-kb {params.distance_kb} \
            --output {output.loci_table} \
            > {log} 2>&1
        
        # Create done marker
        touch {output.done}
        echo "Created loci-based gene table for {wildcards.target}" >> {output.done}
        """
