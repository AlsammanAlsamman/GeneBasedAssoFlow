#!/bin/bash
# Add Gene Names Script
# Maps Entrez Gene IDs to gene symbols using NCBI gene info file

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --genes-out)
            GENES_OUT="$2"
            shift 2
            ;;
        --gene-info)
            GENE_INFO="$2"
            shift 2
            ;;
        --output)
            OUTPUT="$2"
            shift 2
            ;;
        --mapping-output)
            MAPPING_OUTPUT="$2"
            shift 2
            ;;
        *)
            echo "Unknown argument: $1"
            exit 1
            ;;
    esac
done

# Validate required arguments
if [[ -z "$GENES_OUT" || -z "$GENE_INFO" || -z "$OUTPUT" || -z "$MAPPING_OUTPUT" ]]; then
    echo "Error: Missing required arguments"
    echo "Usage: $0 --genes-out FILE --gene-info FILE --output FILE --mapping-output FILE"
    exit 1
fi

# Check if input files exist
if [[ ! -f "$GENES_OUT" ]]; then
    echo "Error: Gene results file not found: $GENES_OUT"
    exit 1
fi

if [[ ! -f "$GENE_INFO" ]]; then
    echo "Error: Gene info file not found: $GENE_INFO"
    echo "Please ensure reference.ncbi_gene_info is correctly configured in analysis.yml"
    exit 1
fi

# Print configuration
echo "========================================="
echo "Adding Gene Names to MAGMA Results"
echo "========================================="
echo "Gene results file: $GENES_OUT"
echo "Gene info file: $GENE_INFO"
echo "Output file: $OUTPUT"
echo "Mapping file: $MAPPING_OUTPUT"
echo "========================================="

# Extract gene ID and symbol (columns 1 and 6) from gene location file
echo "Creating gene ID to name mapping..."
awk '{print $1"\t"$6}' "$GENE_INFO" > "$MAPPING_OUTPUT"

if [[ $? -ne 0 ]]; then
    echo "Error: Failed to create gene mapping file"
    exit 1
fi

# Count mappings
MAPPING_COUNT=$(wc -l < "$MAPPING_OUTPUT" | tr -d ' ')
echo "Created mapping for $MAPPING_COUNT genes"

# Add gene names to MAGMA results
echo "Adding gene names to results..."
awk -v mapping_file="$MAPPING_OUTPUT" '
BEGIN {
    OFS="\t"
    # Read gene ID to name mapping
    while ((getline < mapping_file) > 0) {
        id_to_name[$1] = $2
    }
    close(mapping_file)
}
NR==1 {
    # Print header with GENE_NAME added as second column
    print "GENE", "GENE_NAME", "CHR", "START", "STOP", "NSNPS", "NPARAM", "N", "ZSTAT", "P"
    next
}
{
    gene_id = $1
    gene_name = (gene_id in id_to_name) ? id_to_name[gene_id] : "NA"
    # Print all fields tab-separated
    print gene_id, gene_name, $2, $3, $4, $5, $6, $7, $8, $9
}' "$GENES_OUT" > "$OUTPUT"

if [[ $? -eq 0 ]]; then
    # Count results
    RESULT_COUNT=$(tail -n +2 "$OUTPUT" | wc -l | tr -d ' ')
    NAMED_COUNT=$(tail -n +2 "$OUTPUT" | awk '$2 != "NA"' | wc -l | tr -d ' ')
    NA_COUNT=$(tail -n +2 "$OUTPUT" | awk '$2 == "NA"' | wc -l | tr -d ' ')
    
    echo "========================================="
    echo "Gene names added successfully"
    echo "Total genes: $RESULT_COUNT"
    echo "Genes with names: $NAMED_COUNT"
    echo "Genes without names (NA): $NA_COUNT"
    echo "Output: $OUTPUT"
    echo "Mapping saved to: $MAPPING_OUTPUT"
    echo "========================================="
else
    echo "Error: Failed to add gene names to results"
    exit 1
fi
