"""
Snakemake rule: create_geneset_report
Create Excel report from MAGMA gene-set enrichment analysis results
Generates two sheets: (1) Geneset details with individual genes, (2) Geneset-by-Loci matrix
"""

import sys
sys.path.append("utils")
from bioconfigme import get_results_dir, get_software_module, get_target_analysis, get_dataset, get_analysis_value

# Initialize configs
from bioconfigme import init_configs
init_configs("configs/software.yml", "configs/analysis.yml", "configs/config.yml")

# Get results directory
RESULTS_DIR = get_results_dir()

rule create_geneset_report:
    input:
        genes_with_names = RESULTS_DIR + "/{target}/magma/gene_analysis/{target}_genes_with_names.txt",
        geneset_done = RESULTS_DIR + "/{target}/magma/geneset_analysis/{target}_ALL_geneset_analysis.done",
        add_names_done = RESULTS_DIR + "/{target}/magma/gene_analysis/{target}_add_names.done"
    output:
        report = RESULTS_DIR + "/{target}/magma/geneset_analysis/{target}_geneset_report.xlsx",
        done = RESULTS_DIR + "/{target}/magma/geneset_analysis/{target}_geneset_report.done"
    params:
        results_base = lambda wildcards: f"{RESULTS_DIR}/{wildcards.target}/magma/geneset_analysis",
        target = "{target}",
        gene_loc = get_analysis_value('reference.gene_loc'),
        loci_file = lambda wildcards: get_dataset(get_target_analysis(wildcards.target)['dataset'])['loci_file'],
        genesets = lambda wildcards: ",".join(get_target_analysis(wildcards.target).get('genesets', [])),
        geneset_p_threshold = lambda wildcards: get_target_analysis(wildcards.target).get('geneset_p_threshold', 0.05),
        geneset_gene_p_condition = lambda wildcards: get_target_analysis(wildcards.target).get('geneset_gene_p_condition', 1.0),
        loci_gene_p_threshold = lambda wildcards: get_target_analysis(wildcards.target)['loci_gene_p_thereshold'],
        loci_distance_kb = lambda wildcards: get_target_analysis(wildcards.target)['loci_distance_thereshold_kb'],
        r_module = get_software_module('r'),
        script = "scripts/create_geneset_report.R"
    log:
        RESULTS_DIR + "/{target}/magma/log/create_geneset_report_{target}.log"
    resources:
        mem_mb = 32000,
        time = "00:30:00"
    threads: 2
    shell:
        """
        # Create output directories
        mkdir -p $(dirname {output.report})
        mkdir -p $(dirname {log})
        
        # Load R module
        module load {params.r_module}
        
        # Run geneset report creation script
        Rscript {params.script} \
            --genes-file {input.genes_with_names} \
            --gene-loc {params.gene_loc} \
            --loci-file {params.loci_file} \
            --results-base {params.results_base} \
            --target {params.target} \
            --genesets {params.genesets} \
            --geneset-p-threshold {params.geneset_p_threshold} \
            --geneset-gene-p-condition {params.geneset_gene_p_condition} \
            --loci-gene-p-threshold {params.loci_gene_p_threshold} \
            --loci-distance-kb {params.loci_distance_kb} \
            --output {output.report} \
            > {log} 2>&1
        
        # Create done marker
        touch {output.done}
        
        echo "Gene-set report created: {output.report}"
        """
