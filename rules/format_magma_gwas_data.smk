"""
Snakemake rule: magma_format_data
Format GWAS summary statistics for MAGMA analysis
Creates SNP location and p-value files with optional filtering
"""

import sys
sys.path.append("utils")
from bioconfigme import get_results_dir, get_dataset, get_software_module, get_target_analysis
import json

# Initialize configs
from bioconfigme import init_configs
init_configs("configs/software.yml", "configs/analysis.yml", "configs/config.yml")

# Get all target analyses
from bioconfigme import get_analysis_value
target_analyses = get_analysis_value('target_analysis')
TARGET_NAMES = list(target_analyses.keys()) if target_analyses else []

# Get results directory
RESULTS_DIR = get_results_dir()

rule magma_format_data:
    input:
        gwas_file = lambda wildcards: get_dataset(get_target_analysis(wildcards.target)['dataset'])['file']
    output:
        snploc = RESULTS_DIR + "/{target}/magma/formatted/{target}_snploc.tsv",
        pval = RESULTS_DIR + "/{target}/magma/formatted/{target}_pval.tsv",
        filtered = RESULTS_DIR + "/{target}/magma/formatted/{target}_filtered.tsv",
        done = RESULTS_DIR + "/{target}/magma/formatted/{target}_format.done"
    params:
        dataset_config = lambda wildcards: json.dumps(get_dataset(get_target_analysis(wildcards.target)['dataset'])),
        r_module = get_software_module('r'),
        script = "scripts/format_gwas_for_magma.R"
    log:
        RESULTS_DIR + "/{target}/magma/log/format_data_{target}.log"
    resources:
        mem_mb = 32000,
        time = "00:30:00"
    threads: 2
    shell:
        """
        # Create output directories
        mkdir -p $(dirname {output.snploc})
        mkdir -p $(dirname {log})
        
        # Load R module
        module load {params.r_module}
        
        # Run R script
        Rscript {params.script} \
            --input {input.gwas_file} \
            --output-snploc {output.snploc} \
            --output-pval {output.pval} \
            --output-filtered {output.filtered} \
            --dataset-config '{params.dataset_config}' \
            > {log} 2>&1
        
        # Create done marker
        touch {output.done}
        echo "Formatted GWAS data for {wildcards.target}" >> {output.done}
        """
