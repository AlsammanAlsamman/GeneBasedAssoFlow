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
        gene_cond_p = lambda wildcards: get_target_analysis(wildcards.target).get('geneset_gene_p_condition', 1.0),
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
    run:
        import os
        import subprocess
        
        # Create output directories
        os.makedirs(os.path.dirname(output.gsa_out), exist_ok=True)
        os.makedirs(os.path.dirname(log[0]), exist_ok=True)
        
        # Check if geneset file exists
        if not os.path.exists(params.geneset_file):
            # Create empty output files
            with open(output.gsa_out, 'w') as f:
                f.write(f"# Gene-set file not found: {params.geneset_file}\n")
                f.write("# Skipping analysis for this gene-set\n")
            
            with open(output.done, 'w') as f:
                f.write(f"Gene-set analysis skipped for {wildcards.target} - {wildcards.geneset}\n")
                f.write(f"Reason: Gene-set file not found at {params.geneset_file}\n")
            
            with open(log[0], 'w') as f:
                f.write(f"WARNING: Gene-set file not found: {params.geneset_file}\n")
                f.write("Creating empty output files and skipping analysis.\n")
        else:
            # Run geneset analysis script
            cmd = [
                "bash", params.script,
                "--gene-results", input.genes_raw,
                "--set-annot", params.geneset_file,
                "--output-prefix", params.output_prefix,
                "--gene-cond-p", str(params.gene_cond_p),
                "--magma-module", params.magma_module,
                "--magma-path", params.magma_path,
                "--magma-command", params.magma_command
            ]
            
            with open(log[0], 'w') as log_file:
                result = subprocess.run(cmd, stdout=log_file, stderr=subprocess.STDOUT)
            
            # Check if analysis failed due to gene ID incompatibility
            if result.returncode != 0:
                # Read log to check for gene ID incompatibility
                with open(log[0], 'r') as f:
                    log_content = f.read()
                
                if "no valid gene sets" in log_content or "contains no genes defined in genotype data" in log_content:
                    # Gene-set uses incompatible gene IDs (symbols instead of Entrez IDs)
                    with open(output.gsa_out, 'w') as f:
                        f.write(f"# Gene-set incompatible with MAGMA: {params.geneset_file}\n")
                        f.write("# This gene-set uses gene symbols, but MAGMA requires Entrez Gene IDs\n")
                        f.write("# Analysis skipped - no valid gene sets found\n")
                    
                    with open(output.done, 'w') as f:
                        f.write(f"Gene-set analysis skipped for {wildcards.target} - {wildcards.geneset}\n")
                        f.write(f"Reason: Gene-set uses gene symbols instead of Entrez IDs\n")
                        f.write(f"MAGMA requires gene-sets with Entrez Gene IDs (.entrez.gmt format)\n")
                    
                    # Log is already written, just append summary
                    with open(log[0], 'a') as f:
                        f.write("\n" + "="*50 + "\n")
                        f.write("SUMMARY: Gene-set incompatible with MAGMA\n")
                        f.write("Gene-set uses gene symbols, but MAGMA requires Entrez Gene IDs\n")
                        f.write("Creating empty output files and marking as skipped.\n")
                else:
                    # Real error - raise exception
                    raise Exception(f"Gene-set analysis failed for {wildcards.geneset}")
            else:
                # Create done marker for successful analysis
                with open(output.done, 'w') as f:
                    f.write(f"Gene-set analysis completed for {wildcards.target} - {wildcards.geneset}\n")

ruleorder: all_geneset_analyses > magma_geneset_analysis

rule all_geneset_analyses:
    input:
        lambda wildcards: expand(
            RESULTS_DIR + "/{target}/magma/geneset_analysis/{target}_{geneset}_geneset_analysis.done",
            target=wildcards.target,
            geneset=get_target_analysis(wildcards.target).get('genesets', [])
        )
    output:
        done = RESULTS_DIR + "/{target}/magma/geneset_analysis/{target}_ALL_geneset_analysis.done"
    shell:
        """
        touch {output.done}
        echo "All gene-set analyses completed for {wildcards.target}" >> {output.done}
        """
