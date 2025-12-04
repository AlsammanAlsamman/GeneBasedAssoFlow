#!/bin/bash
# Pipeline execution steps

# Step 1: Format GWAS data for MAGMA
./submit.sh --snakefile rules/format_gwas_data.smk /s/nath-lab/alsamman/___Analysis___/Hispanic_2025/MAGMA_PipiLIne/results/formatted/Hisp_0.1_OneMissing_format.done
