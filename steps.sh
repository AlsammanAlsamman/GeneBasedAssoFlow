#!/bin/bash
# Pipeline execution steps

# Step 1: Format GWAS data for MAGMA
# Structure: results/{target_analysis}/magma/formatted/
./submit.sh --snakefile rules/format_magma_gwas_data.smk /s/nath-lab/alsamman/___Analysis___/Hispanic_2025/MAGMA_PipiLIne/results/example_analysis/magma/formatted/example_analysis_format.done

# Step 2: Annotate SNPs to genes
./submit.sh --snakefile rules/magma_annotate.smk /s/nath-lab/alsamman/___Analysis___/Hispanic_2025/MAGMA_PipiLIne/results/example_analysis/magma/annotation/example_analysis_annotate.done

# Step 3: Gene-based association analysis
./submit.sh --snakefile rules/magma_gene_analysis.smk /s/nath-lab/alsamman/___Analysis___/Hispanic_2025/MAGMA_PipiLIne/results/example_analysis/magma/gene_analysis/example_analysis_gene_analysis.done

# Step 4: Gene-set enrichment analysis (hallmark)
./submit.sh --snakefile rules/magma_geneset_analysis.smk /s/nath-lab/alsamman/___Analysis___/Hispanic_2025/MAGMA_PipiLIne/results/example_analysis/magma/geneset_analysis/example_analysis_msigdb_hallmark_geneset_analysis.done

# Step 4b: Gene-set enrichment analysis (GO biological process)
./submit.sh --snakefile rules/magma_geneset_analysis.smk /s/nath-lab/alsamman/___Analysis___/Hispanic_2025/MAGMA_PipiLIne/results/example_analysis/magma/geneset_analysis/example_analysis_msigdb_go_bp_geneset_analysis.done

# Step 5: Add gene names to results
./submit.sh --snakefile rules/add_gene_names.smk /s/nath-lab/alsamman/___Analysis___/Hispanic_2025/MAGMA_PipiLIne/results/example_analysis/magma/gene_analysis/example_analysis_add_names.done

# Step 6: Create loci-based gene table
./submit.sh --snakefile rules/create_loci_table.smk /s/nath-lab/alsamman/___Analysis___/Hispanic_2025/MAGMA_PipiLIne/results/example_analysis/magma/loci_analysis/example_analysis_loci_table.done
