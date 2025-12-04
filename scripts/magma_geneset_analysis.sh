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
echo "========================================="

# Run MAGMA gene-set analysis
echo "Running MAGMA gene-set analysis..."
$MAGMA_EXEC \
    --gene-results "$GENE_RESULTS" \
    --set-annot "$SET_ANNOT" \
    --out "$OUTPUT_PREFIX"

# Check if analysis was successful
if [[ $? -eq 0 ]]; then
    echo "========================================="
    echo "Gene-set analysis completed successfully"
    echo "Output: ${OUTPUT_PREFIX}.gsa.out"
    echo "========================================="
else
    echo "Error: MAGMA gene-set analysis failed"
    exit 1
fi
