#!/bin/bash
# MAGMA Gene-Set Analysis Script
# Tests for enrichment of associations in predefined gene-sets

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --gene-results)
            GENE_RESULTS="$2"
            shift 2
            ;;
        --set-annot)
            SET_ANNOT="$2"
            shift 2
            ;;
        --output-prefix)
            OUTPUT_PREFIX="$2"
            shift 2
            ;;
        --gene-cond-p)
            GENE_COND_P="$2"
            shift 2
            ;;
        --magma-module)
            MAGMA_MODULE="$2"
            shift 2
            ;;
        --magma-path)
            MAGMA_PATH="$2"
            shift 2
            ;;
        --magma-command)
            MAGMA_COMMAND="$2"
            shift 2
            ;;
        *)
            echo "Unknown argument: $1"
            exit 1
            ;;
    esac
done

# Set default for GENE_COND_P if not provided
if [[ -z "$GENE_COND_P" ]]; then
    GENE_COND_P="1.0"
fi

# Validate required arguments
if [[ -z "$GENE_RESULTS" || -z "$SET_ANNOT" || -z "$OUTPUT_PREFIX" ]]; then
    echo "Error: Missing required arguments"
    echo "Usage: $0 --gene-results FILE --set-annot FILE --output-prefix PREFIX --magma-module MODULE --magma-path PATH --magma-command CMD"
    exit 1
fi

# Determine MAGMA executable
if [[ -n "$MAGMA_PATH" && "$MAGMA_PATH" != "None" ]]; then
    # Use direct path if provided
    MAGMA_EXEC="$MAGMA_PATH"
    echo "Using MAGMA from path: $MAGMA_EXEC"
else
    # Load module and use command name
    if [[ -n "$MAGMA_MODULE" && "$MAGMA_MODULE" != "None" ]]; then
        echo "Loading MAGMA module: $MAGMA_MODULE"
        module load "$MAGMA_MODULE"
    fi
    MAGMA_EXEC="$MAGMA_COMMAND"
    echo "Using MAGMA command: $MAGMA_EXEC"
fi

# Print configuration
echo "========================================="
echo "MAGMA Gene-Set Enrichment Analysis"
echo "========================================="
echo "Gene results file: $GENE_RESULTS"
echo "Gene-set annotation: $SET_ANNOT"
echo "Output prefix: $OUTPUT_PREFIX"
echo "Gene P-value condition: $GENE_COND_P"
echo "========================================="

# Check if we need to filter genes by P-value
FILTERED_GENES="$GENE_RESULTS"
if (( $(echo "$GENE_COND_P < 1.0" | bc -l) )); then
    echo "Filtering genes: only genes with P < $GENE_COND_P will be included"
    
    # Create filtered gene results file
    FILTERED_GENES="${OUTPUT_PREFIX}_filtered.genes.raw"
    
    # Read header and data, filter by P-value
    # .genes.raw format: GENE CHR START STOP NSNPS NPARAM N ZSTAT P
    head -1 "$GENE_RESULTS" > "$FILTERED_GENES"
    awk -v threshold="$GENE_COND_P" 'NR>1 && $9 < threshold' "$GENE_RESULTS" >> "$FILTERED_GENES"
    
    NGENES_BEFORE=$(tail -n +2 "$GENE_RESULTS" | wc -l)
    NGENES_AFTER=$(tail -n +2 "$FILTERED_GENES" | wc -l)
    echo "Genes before filtering: $NGENES_BEFORE"
    echo "Genes after filtering (P < $GENE_COND_P): $NGENES_AFTER"
    
    # Check if any genes pass the threshold
    if [ "$NGENES_AFTER" -eq 0 ]; then
        echo "WARNING: No genes pass the P-value threshold!"
        echo "Gene-set analysis will likely find no significant results."
    fi
else
    echo "Including all genes in gene-set analysis (no P-value conditioning)"
fi

# Run MAGMA gene-set analysis
echo "Running MAGMA gene-set analysis..."
echo "Command: $MAGMA_EXEC --gene-results $FILTERED_GENES --set-annot $SET_ANNOT --out $OUTPUT_PREFIX"

$MAGMA_EXEC \
    --gene-results "$FILTERED_GENES" \
    --set-annot "$SET_ANNOT" \
    --out "$OUTPUT_PREFIX"

# Check if analysis was successful
MAGMA_EXIT_CODE=$?

# Clean up filtered file if created
if [ "$FILTERED_GENES" != "$GENE_RESULTS" ] && [ -f "$FILTERED_GENES" ]; then
    rm -f "$FILTERED_GENES"
fi

if [[ $MAGMA_EXIT_CODE -eq 0 ]]; then
    echo "========================================="
    echo "Gene-set analysis completed successfully"
    echo "Output: ${OUTPUT_PREFIX}.gsa.out"
    echo "========================================="
else
    echo "Error: MAGMA gene-set analysis failed with exit code $MAGMA_EXIT_CODE"
    exit 1
fi
