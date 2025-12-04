#!/bin/bash
# MAGMA Gene-Based Analysis Script
# Aggregates SNP p-values to gene-level associations using LD reference panel

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --bfile)
            BFILE="$2"
            shift 2
            ;;
        --pval)
            PVAL="$2"
            shift 2
            ;;
        --gene-annot)
            GENE_ANNOT="$2"
            shift 2
            ;;
        --sample-size)
            SAMPLE_SIZE="$2"
            shift 2
            ;;
        --gene-model)
            GENE_MODEL="$2"
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
if [[ -z "$BFILE" || -z "$PVAL" || -z "$GENE_ANNOT" || -z "$SAMPLE_SIZE" || -z "$GENE_MODEL" || -z "$OUTPUT_PREFIX" ]]; then
    echo "Error: Missing required arguments"
    echo "Usage: $0 --bfile FILE --pval FILE --gene-annot FILE --sample-size INT --gene-model MODEL --output-prefix PREFIX --magma-module MODULE --magma-path PATH --magma-command CMD"
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
echo "MAGMA Gene-Based Analysis"
echo "========================================="
echo "Reference panel (bfile): $BFILE"
echo "P-value file: $PVAL"
echo "Gene annotation: $GENE_ANNOT"
echo "Sample size: $SAMPLE_SIZE"
echo "Gene model: $GENE_MODEL"
echo "Output prefix: $OUTPUT_PREFIX"
echo "========================================="

# Run MAGMA gene-based analysis
echo "Running MAGMA gene-based analysis..."
$MAGMA_EXEC \
    --bfile "$BFILE" \
    --pval "$PVAL" N=$SAMPLE_SIZE \
    --gene-annot "$GENE_ANNOT" \
    --gene-model "$GENE_MODEL" \
    --out "$OUTPUT_PREFIX"

# Check if analysis was successful
if [[ $? -eq 0 ]]; then
    echo "========================================="
    echo "Gene-based analysis completed successfully"
    echo "Output: ${OUTPUT_PREFIX}.genes.out"
    echo "========================================="
else
    echo "Error: MAGMA gene-based analysis failed"
    exit 1
fi
