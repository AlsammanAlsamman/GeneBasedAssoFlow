"""
Snakemake rule for comparing SNP-based GWAS with gene-based MAGMA results
Generates Manhattan plots and Miami plots (mirrored comparison)
Creates both genome-wide and per-chromosome visualizations
"""

import sys
sys.path.append("utils")
from bioconfigme import (
    init_configs,
    get_results_dir,
    get_analysis_value,
    get_target_analysis,
    get_dataset,
    get_software_module
)

# Initialize configurations
init_configs("configs/software.yml", "configs/analysis.yml", "configs/config.yml")

RESULTS_DIR = get_results_dir()

rule plot_snp_gene_comparison:
    input:
        gwas = lambda wildcards: get_dataset(get_target_analysis(wildcards.target)['dataset'])['file'],
        genes = RESULTS_DIR + "/{target}/magma/gene_analysis/{target}.genes.out",
        gene_loc = get_analysis_value('reference.gene_loc'),
        loci = lambda wildcards: get_dataset(get_target_analysis(wildcards.target)['dataset'])['loci_file'],
        report_done = RESULTS_DIR + "/{target}/magma/geneset_analysis/{target}_geneset_report.done"
    output:
        # Genome-wide plots (manhattantwin adds _pub.png suffix)
        manhattan_gw = RESULTS_DIR + "/{target}/plots/genome_wide/{target}_snp_pub.png",
        miami_gw = RESULTS_DIR + "/{target}/plots/genome_wide/{target}_miami.png",
        # Filtered data files
        snp_filtered = RESULTS_DIR + "/{target}/plots/filtered_data/{target}_snps_filtered.txt",
        gene_filtered = RESULTS_DIR + "/{target}/plots/filtered_data/{target}_genes_filtered.txt",
        # Done marker
        done = RESULTS_DIR + "/{target}/plots/{target}_plots.done"
    params:
        analysis = "{target}",
        snp_p = lambda wildcards: get_target_analysis(wildcards.target).get('snp_plot_p_threshold', 5e-2),
        gene_p = lambda wildcards: get_target_analysis(wildcards.target).get('gene_plot_p_threshold', 5e-8),
        cluster_dist = lambda wildcards: get_target_analysis(wildcards.target).get('snp_cluster_distance', 250000),
        n_cases = lambda wildcards: get_dataset(get_target_analysis(wildcards.target)['dataset']).get('cases', 'NA'),
        n_controls = lambda wildcards: get_dataset(get_target_analysis(wildcards.target)['dataset']).get('controls', 'NA'),
        plot_per_chr = lambda wildcards: get_target_analysis(wildcards.target).get('plot_per_chromosome', True),
        plot_only_sig = lambda wildcards: get_target_analysis(wildcards.target).get('plot_only_sig_chromosomes', True),
        genome_wide_dir = RESULTS_DIR + "/{target}/plots/genome_wide",
        per_chr_dir = RESULTS_DIR + "/{target}/plots/per_chromosome",
        filtered_dir = RESULTS_DIR + "/{target}/plots/filtered_data",
        r_module = get_software_module('r'),
        script = "scripts/plot_snp_gene_comparison.R"
    log:
        RESULTS_DIR + "/{target}/magma/log/plot_snp_gene_comparison_{target}.log"
    resources:
        mem_mb = 16000,
        time = "01:00:00"
    threads: 2
    shell:
        """
        # Create output directories
        mkdir -p {params.genome_wide_dir}
        mkdir -p {params.per_chr_dir}
        mkdir -p {params.filtered_dir}
        mkdir -p $(dirname {log})
        
        # Load R module
        module load {params.r_module}
        
        # Run plotting script
        Rscript {params.script} \
            --gwas {input.gwas} \
            --genes {input.genes} \
            --gene_loc {input.gene_loc} \
            --loci {input.loci} \
            --analysis {params.analysis} \
            --snp_p_threshold {params.snp_p} \
            --gene_p_threshold {params.gene_p} \
            --cluster_distance {params.cluster_dist} \
            --n_cases {params.n_cases} \
            --n_controls {params.n_controls} \
            --plot_per_chromosome {params.plot_per_chr} \
            --plot_only_sig_chromosomes {params.plot_only_sig} \
            --genome_wide_dir {params.genome_wide_dir} \
            --per_chromosome_dir {params.per_chr_dir} \
            --filtered_data_dir {params.filtered_dir} \
            --output_done {output.done} \
            > {log} 2>&1
        
        # Create done marker
        touch {output.done}
        
        echo "SNP-Gene comparison plots created successfully"
        """
