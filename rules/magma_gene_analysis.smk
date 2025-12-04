"""
Snakemake rule: magma_gene_analysis
Perform gene-based association analysis using MAGMA
Aggregates SNP p-values to gene-level associations using LD reference panel
"""

import sys
sys.path.append("utils")
from bioconfigme import get_results_dir, get_software_module, get_software_path, get_software_command, get_target_analysis, get_dataset, get_ref_panel

# Initialize configs
from bioconfigme import init_configs
init_configs("configs/software.yml", "configs/analysis.yml", "configs/config.yml")

# Get results directory
RESULTS_DIR = get_results_dir()

rule magma_gene_analysis:
    input:
        pval = RESULTS_DIR + "/{target}/magma/formatted/{target}_pval.tsv",
        annot = RESULTS_DIR + "/{target}/magma/annotation/{target}.genes.annot",
        done_marker = RESULTS_DIR + "/{target}/magma/annotation/{target}_annotate.done"
    output:
        genes_out = RESULTS_DIR + "/{target}/magma/gene_analysis/{target}.genes.out",
        done = RESULTS_DIR + "/{target}/magma/gene_analysis/{target}_gene_analysis.done"
    params:
        ref_panel = lambda wildcards: get_ref_panel(get_target_analysis(wildcards.target)['ref_panel']),
        sample_size = lambda wildcards: get_dataset(get_target_analysis(wildcards.target)['dataset'])['samples'],
        gene_model = lambda wildcards: get_target_analysis(wildcards.target)['gene_model'],
        magma_module = get_software_module('magma'),
        magma_path = get_software_path('magma'),
        magma_command = get_software_command('magma'),
        output_prefix = lambda wildcards: f"{RESULTS_DIR}/{wildcards.target}/magma/gene_analysis/{wildcards.target}",
        script = "scripts/magma_gene_analysis.sh"
    log:
        RESULTS_DIR + "/{target}/magma/log/gene_analysis_{target}.log"
    resources:
        mem_mb = 32000,
        time = "00:30:00"
    threads: 2
    shell:
        """
        # Create output directories
        mkdir -p $(dirname {output.genes_out})
        mkdir -p $(dirname {log})
        
        # Run gene analysis script
        bash {params.script} \
            --bfile {params.ref_panel} \
            --pval {input.pval} \
            --gene-annot {input.annot} \
            --sample-size {params.sample_size} \
            --gene-model "{params.gene_model}" \
            --output-prefix {params.output_prefix} \
            --magma-module "{params.magma_module}" \
            --magma-path "{params.magma_path}" \
            --magma-command "{params.magma_command}" \
            > {log} 2>&1
        
        # Create done marker
        touch {output.done}
        echo "Gene-based analysis completed for {wildcards.target}" >> {output.done}
        """
