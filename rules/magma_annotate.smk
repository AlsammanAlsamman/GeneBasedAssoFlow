"""
Snakemake rule: magma_annotate
Annotate SNPs to genes using MAGMA
Maps SNPs to genes based on genomic location with configurable window
"""

import sys
sys.path.append("utils")
from bioconfigme import get_results_dir, get_software_module, get_software_path, get_software_command, get_target_analysis, get_analysis_value

# Initialize configs
from bioconfigme import init_configs
init_configs("configs/software.yml", "configs/analysis.yml", "configs/config.yml")

# Get results directory
RESULTS_DIR = get_results_dir()

rule magma_annotate:
    input:
        snploc = RESULTS_DIR + "/{target}/magma/formatted/{target}_snploc.tsv",
        done_marker = RESULTS_DIR + "/{target}/magma/formatted/{target}_format.done"
    output:
        annot = RESULTS_DIR + "/{target}/magma/annotation/{target}.genes.annot",
        done = RESULTS_DIR + "/{target}/magma/annotation/{target}_annotate.done"
    params:
        gene_loc = lambda wildcards: get_analysis_value('reference.gene_loc'),
        window = lambda wildcards: get_target_analysis(wildcards.target)['annotation_window'],
        magma_module = get_software_module('magma'),
        magma_path = get_software_path('magma'),
        magma_command = get_software_command('magma'),
        output_prefix = lambda wildcards: f"{RESULTS_DIR}/{wildcards.target}/magma/annotation/{wildcards.target}",
        script = "scripts/magma_annotate.sh"
    log:
        RESULTS_DIR + "/{target}/magma/log/annotate_{target}.log"
    resources:
        mem_mb = 32000,
        time = "00:30:00"
    threads: 2
    shell:
        """
        # Create output directories
        mkdir -p $(dirname {output.annot})
        mkdir -p $(dirname {log})
        
        # Run annotation script
        bash {params.script} \
            --snploc {input.snploc} \
            --gene-loc {params.gene_loc} \
            --window {params.window} \
            --output-prefix {params.output_prefix} \
            --magma-module "{params.magma_module}" \
            --magma-path "{params.magma_path}" \
            --magma-command "{params.magma_command}" \
            > {log} 2>&1
        
        # Create done marker
        touch {output.done}
        echo "Annotated SNPs to genes for {wildcards.target}" >> {output.done}
        """
