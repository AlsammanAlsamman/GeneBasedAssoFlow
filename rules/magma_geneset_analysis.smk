"""
Snakemake rule: magma_geneset_analysis
Perform gene-set enrichment analysis using MAGMA
Tests for enrichment of associations in predefined gene-sets
"""

import sys
sys.path.append("utils")
from bioconfigme import get_results_dir, get_software_module, get_software_path, get_software_command, get_target_analysis, get_database

# Initialize configs
from bioconfigme import init_configs
init_configs("configs/software.yml", "configs/analysis.yml", "configs/config.yml")

# Get results directory
RESULTS_DIR = get_results_dir()

rule magma_geneset_analysis:
    input:
        genes_raw = RESULTS_DIR + "/{target}/magma/gene_analysis/{target}.genes.raw",
        done_marker = RESULTS_DIR + "/{target}/magma/gene_analysis/{target}_gene_analysis.done"
    output:
        gsa_out = RESULTS_DIR + "/{target}/magma/geneset_analysis/{target}_{geneset}.gsa.out",
        done = RESULTS_DIR + "/{target}/magma/geneset_analysis/{target}_{geneset}_geneset_analysis.done"
    params:
        geneset_file = lambda wildcards: get_database(wildcards.geneset),
        magma_module = get_software_module('magma'),
        magma_path = get_software_path('magma'),
        magma_command = get_software_command('magma'),
        output_prefix = lambda wildcards: f"{RESULTS_DIR}/{wildcards.target}/magma/geneset_analysis/{wildcards.target}_{wildcards.geneset}",
        script = "scripts/magma_geneset_analysis.sh"
    log:
        RESULTS_DIR + "/{target}/magma/log/geneset_analysis_{target}_{geneset}.log"
    resources:
        mem_mb = 32000,
        time = "00:30:00"
    threads: 2
    shell:
        """
        # Create output directories
        mkdir -p $(dirname {output.gsa_out})
        mkdir -p $(dirname {log})
        
        # Run geneset analysis script
        bash {params.script} \
            --gene-results {input.genes_raw} \
            --set-annot {params.geneset_file} \
            --output-prefix {params.output_prefix} \
            --magma-module "{params.magma_module}" \
            --magma-path "{params.magma_path}" \
            --magma-command "{params.magma_command}" \
            > {log} 2>&1
        
        # Create done marker
        touch {output.done}
        echo "Gene-set analysis completed for {wildcards.target} - {wildcards.geneset}" >> {output.done}
        """
