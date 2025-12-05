#!/usr/bin/env python3
"""
Format gene-set GMT files: convert gene symbols to Entrez IDs
Creates backups of original files with .original suffix
"""

import sys
sys.path.append("utils")
from bioconfigme import get_config, get_results_dir

# Load configuration
config_data = get_config()

# Get all GMT database paths
databases = config_data.get("databases", {})
gmt_files = list(databases.values())

# Get gene reference file
gene_ref = config_data["reference"]["gene_loc"]

# Output
results_dir = get_results_dir()

rule format_genesets:
    """
    Convert gene-set GMT files from symbols to Entrez IDs
    Backs up originals and creates conversion report
    """
    input:
        gmt_files = gmt_files,
        gene_ref = gene_ref
    output:
        done = f"{results_dir}/format_genesets.done",
        report = f"{results_dir}/log/format_genesets_report.txt"
    params:
        gmt_list = ",".join(gmt_files)
    log:
        f"{results_dir}/log/format_genesets.log"
    resources:
        mem_mb = 8000,
        time = "00:20:00"
    threads: 1
    shell:
        """
        python scripts/format_genesets.py \
            --gmt-files {params.gmt_list} \
            --gene-ref {input.gene_ref} \
            --report {output.report} \
            --done {output.done} \
            2>&1 | tee {log}
        """
