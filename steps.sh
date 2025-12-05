#!/bin/bash
# Pipeline execution steps

# Step 0: Format gene-set GMT files (convert symbols to Entrez IDs)
./submit.sh --snakefile rules/format_genesets.smk /s/nath-lab/alsamman/___Analysis___/Hispanic_2025/MAGMA_PipiLIne/results/format_genesets.done

# Step 1: Format GWAS data for MAGMA
# Structure: results/{target_analysis}/magma/formatted/
./submit.sh --snakefile rules/format_magma_gwas_data.smk /s/nath-lab/alsamman/___Analysis___/Hispanic_2025/MAGMA_PipiLIne/results/example_analysis/magma/formatted/example_analysis_format.done

# Step 2: Annotate SNPs to genes
./submit.sh --snakefile rules/magma_annotate.smk /s/nath-lab/alsamman/___Analysis___/Hispanic_2025/MAGMA_PipiLIne/results/example_analysis/magma/annotation/example_analysis_annotate.done

# Step 3: Gene-based association analysis
./submit.sh --snakefile rules/magma_gene_analysis.smk /s/nath-lab/alsamman/___Analysis___/Hispanic_2025/MAGMA_PipiLIne/results/example_analysis/magma/gene_analysis/example_analysis_gene_analysis.done

# Step 4: Gene-set enrichment analysis (all databases)
# This will run all configured gene-set databases: hallmark, GO BP, BioCarta, KEGG, Reactome, ImmuneSigDB, Cell Type
./submit.sh --snakefile rules/magma_geneset_analysis.smk /s/nath-lab/alsamman/___Analysis___/Hispanic_2025/MAGMA_PipiLIne/results/example_analysis/magma/geneset_analysis/example_analysis_ALL_geneset_analysis.done

# Step 5: Add gene names to results
./submit.sh --snakefile rules/add_gene_names.smk /s/nath-lab/alsamman/___Analysis___/Hispanic_2025/MAGMA_PipiLIne/results/example_analysis/magma/gene_analysis/example_analysis_add_names.done

# Step 6: Create loci-based gene table
./submit.sh --snakefile rules/create_loci_table.smk /s/nath-lab/alsamman/___Analysis___/Hispanic_2025/MAGMA_PipiLIne/results/example_analysis/magma/loci_analysis/example_analysis_loci_table.done

# Step 7: Create gene-set enrichment report
./submit.sh --snakefile rules/create_geneset_report.smk /s/nath-lab/alsamman/___Analysis___/Hispanic_2025/MAGMA_PipiLIne/results/example_analysis/magma/geneset_analysis/example_analysis_geneset_report.done
