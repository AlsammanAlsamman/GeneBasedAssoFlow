#!/bin/bash
# MAGMA SNP Annotation Script
# Annotates SNPs to genes based on genomic location

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --snploc)
            SNPLOC="$2"
            shift 2
            ;;
        --gene-loc)
            GENE_LOC="$2"
            shift 2
            ;;
        --window)
            WINDOW="$2"
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
if [[ -z "$SNPLOC" || -z "$GENE_LOC" || -z "$WINDOW" || -z "$OUTPUT_PREFIX" ]]; then
    echo "Error: Missing required arguments"
    echo "Usage: $0 --snploc FILE --gene-loc FILE --window INT --output-prefix PREFIX --magma-module MODULE --magma-path PATH --magma-command CMD"
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
echo "MAGMA SNP Annotation"
echo "========================================="
echo "SNP location file: $SNPLOC"
echo "Gene location file: $GENE_LOC"
echo "Annotation window: ${WINDOW}kb"
echo "Output prefix: $OUTPUT_PREFIX"
echo "========================================="

# Run MAGMA annotation
echo "Running MAGMA annotation..."
$MAGMA_EXEC \
    --annotate window=$WINDOW \
    --snp-loc "$SNPLOC" \
    --gene-loc "$GENE_LOC" \
    --out "$OUTPUT_PREFIX"

# Check if annotation was successful
if [[ $? -eq 0 ]]; then
    echo "========================================="
    echo "SNP annotation completed successfully"
    echo "Output: ${OUTPUT_PREFIX}.genes.annot"
    echo "========================================="
else
    echo "Error: MAGMA annotation failed"
    exit 1
fi
