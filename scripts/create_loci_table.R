#!/usr/bin/env Rscript
# Create Loci-Based Gene Table
# Maps significant genes to genomic loci with distance buffer

# Load required libraries
suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx2)
  library(optparse)
})

# Parse command line arguments
option_list <- list(
  make_option("--genes-file", type="character", help="Path to genes_with_names.txt file"),
  make_option("--gene-loc", type="character", help="Path to gene location file"),
  make_option("--loci-file", type="character", help="Path to loci info file"),
  make_option("--p-threshold", type="numeric", help="P-value threshold for gene inclusion"),
  make_option("--distance-kb", type="numeric", help="Distance threshold in kb for loci matching"),
  make_option("--output", type="character", help="Output Excel file path")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Validate arguments
if (is.null(opt$`genes-file`) || is.null(opt$`gene-loc`) || is.null(opt$`loci-file`) || 
    is.null(opt$`p-threshold`) || is.null(opt$`distance-kb`) || is.null(opt$output)) {
  stop("All arguments are required!")
}

cat("========================================\n")
cat("Creating Loci-Based Gene Table\n")
cat("========================================\n")
cat("Genes file:", opt$`genes-file`, "\n")
cat("Gene location file:", opt$`gene-loc`, "\n")
cat("Loci file:", opt$`loci-file`, "\n")
cat("P-value threshold:", opt$`p-threshold`, "\n")
cat("Distance threshold:", opt$`distance-kb`, "kb\n")
cat("Output file:", opt$output, "\n")
cat("========================================\n\n")

# Read gene results with names
cat("Reading gene results...\n")
genes <- fread(opt$`genes-file`, header=TRUE, sep="\t")
cat("Total genes:", nrow(genes), "\n")

# Print raw column names for debugging
cat("Raw columns:", paste(names(genes), collapse=" | "), "\n")

# Clean column names (remove extra whitespace and standardize)
old_names <- names(genes)
new_names <- trimws(old_names)
new_names <- gsub("\\s+", "_", new_names)  # Replace multiple spaces with underscore
setnames(genes, old_names, new_names)

cat("Cleaned columns:", paste(names(genes), collapse=", "), "\n")

# Find P-value column (could be just "P" or with extra text)
p_col_idx <- which(names(genes) == "P")
if (length(p_col_idx) == 0) {
  # Try to find any column ending with P
  p_col_idx <- grep("P$", names(genes))
  if (length(p_col_idx) == 0) {
    stop("Cannot find P-value column!")
  }
}

p_col_name <- names(genes)[p_col_idx[1]]
cat("Using P-value column:", p_col_name, "\n")

# Filter by p-value threshold
genes_sig <- genes[get(p_col_name) <= opt$`p-threshold`]
cat("Significant genes (P <=", opt$`p-threshold`, "):", nrow(genes_sig), "\n")

if (nrow(genes_sig) == 0) {
  cat("WARNING: No significant genes found. Creating empty output.\n")
  wb <- wb_workbook()
  wb$add_worksheet("Loci_Genes")
  wb$add_data(sheet=1, x=data.frame(Message="No significant genes found at specified threshold"))
  wb$save(opt$output)
  cat("Empty Excel file created.\n")
  quit(save="no", status=0)
}

# Read gene location file
cat("\nReading gene location file...\n")
gene_loc <- fread(opt$`gene-loc`, header=FALSE, 
                  col.names=c("GENE", "CHR", "START", "END", "STRAND", "GENE_NAME"))
cat("Gene locations loaded:", nrow(gene_loc), "\n")

# Merge gene results with locations
cat("\nMerging gene results with locations...\n")
# Note: genes_sig already has CHR, START, STOP from MAGMA output
# We need gene locations from gene_loc for precise positioning
# Use suffix to distinguish if there are duplicates
genes_with_pos <- merge(genes_sig, 
                        gene_loc[, .(GENE, CHR_loc=CHR, START_loc=START, END_loc=END)], 
                        by="GENE", all.x=TRUE)

# Use gene_loc coordinates (more accurate)
genes_with_pos[, ':='(CHR_final = CHR_loc, START_final = START_loc, END_final = END_loc)]
genes_with_pos[, ':='(CHR_loc = NULL, START_loc = NULL, END_loc = NULL)]

cat("Genes with positions:", nrow(genes_with_pos[!is.na(CHR_final)]), "\n")

# Read loci file
cat("\nReading loci file...\n")
loci <- fread(opt$`loci-file`, header=TRUE)
cat("Loci loaded:", nrow(loci), "\n")

# Ensure column names (handle variations)
if ("locus" %in% names(loci)) setnames(loci, "locus", "Locus", skip_absent=TRUE)
if (!"Locus" %in% names(loci) && "Locus" %in% names(loci)) {
  # Already correct
} else if ("locus" %in% tolower(names(loci))) {
  setnames(loci, grep("locus", names(loci), ignore.case=TRUE)[1], "Locus")
}

# Calculate distance buffer in bp
distance_bp <- opt$`distance-kb` * 1000

# Match genes to loci with distance buffer
cat("\nMatching genes to loci (with ±", opt$`distance-kb`, "kb buffer)...\n")

matched_list <- list()

for (i in 1:nrow(loci)) {
  locus_name <- loci$Locus[i]
  locus_chr <- loci$chr[i]
  locus_start <- loci$start[i] - distance_bp
  locus_end <- loci$end[i] + distance_bp
  
  # Find genes overlapping this locus region (with buffer)
  matched_genes <- genes_with_pos[
    CHR_final == locus_chr & 
    !is.na(START_final) & !is.na(END_final) &
    ((START_final >= locus_start & START_final <= locus_end) |
     (END_final >= locus_start & END_final <= locus_end) |
     (START_final <= locus_start & END_final >= locus_end))
  ]
  
  if (nrow(matched_genes) > 0) {
    # Add all locus information from loci file
    matched_genes[, Locus := locus_name]
    matched_genes[, Locus_Chr := locus_chr]
    matched_genes[, Locus_Start := loci$start[i]]
    matched_genes[, Locus_End := loci$end[i]]
    
    # Add any additional columns from loci file
    loci_cols <- setdiff(names(loci), c("Locus", "chr", "start", "end"))
    for (col in loci_cols) {
      matched_genes[, paste0("Locus_", col) := loci[[col]][i]]
    }
    
    matched_list[[i]] <- matched_genes
  }
}

# Combine all matches
if (length(matched_list) > 0) {
  result <- rbindlist(matched_list, fill=TRUE)
  
  # Sort by locus and P-value
  result <- result[order(Locus_Chr, Locus_Start, get(p_col_name))]
  
  # Reorder columns: Locus info first (including all Locus_ columns), then gene info
  locus_info_cols <- c("Locus", "Locus_Chr", "Locus_Start", "Locus_End")
  additional_locus_cols <- grep("^Locus_", names(result), value=TRUE)
  additional_locus_cols <- setdiff(additional_locus_cols, locus_info_cols)
  
  gene_info_cols <- c("GENE", "GENE_NAME", "CHR_final", "START_final", "END_final", p_col_name)
  other_cols <- setdiff(names(result), c(locus_info_cols, additional_locus_cols, gene_info_cols))
  
  setcolorder(result, c(locus_info_cols, additional_locus_cols, gene_info_cols, other_cols))
  
  cat("Total gene-locus matches:", nrow(result), "\n")
  
  # Create Excel workbook
  cat("\nCreating Excel file...\n")
  wb <- wb_workbook()
  wb$add_worksheet("Loci_Genes")
  
  # Add data
  wb$add_data(sheet=1, x=result, start_col=1, start_row=1)
  
  # Style header
  wb$add_fill(sheet=1, dims="A1:Z1", color=wb_color(hex="4472C4"))
  wb$add_font(sheet=1, dims="A1:Z1", color=wb_color(hex="FFFFFF"), bold=TRUE)
  
  # Merge cells for loci information
  cat("Applying merged cells for loci regions...\n")
  unique_loci <- unique(result$Locus)
  
  # Find how many Locus columns we have
  locus_col_count <- length(c(locus_info_cols, additional_locus_cols))
  
  current_row <- 2  # Start from row 2 (after header)
  for (locus in unique_loci) {
    locus_rows <- which(result$Locus == locus)
    n_genes <- length(locus_rows)
    
    if (n_genes > 1) {
      end_row <- current_row + n_genes - 1
      # Merge all locus info columns
      for (col_idx in 1:locus_col_count) {
        wb$merge_cells(sheet=1, rows=current_row:end_row, cols=col_idx)
      }
    }
    
    current_row <- current_row + n_genes
  }
  
  # Auto-size columns
  wb$set_col_widths(sheet=1, cols=1:ncol(result), widths="auto")
  
  # Save workbook
  wb$save(opt$output)
  
  cat("========================================\n")
  cat("Loci-based gene table created successfully!\n")
  cat("Output:", opt$output, "\n")
  cat("Unique loci:", length(unique_loci), "\n")
  cat("Total genes:", nrow(result), "\n")
  cat("========================================\n")
  
} else {
  cat("WARNING: No genes matched to any loci.\n")
  wb <- wb_workbook()
  wb$add_worksheet("Loci_Genes")
  wb$add_data(sheet=1, x=data.frame(Message="No genes matched to loci regions"))
  wb$save(opt$output)
  cat("Empty Excel file created.\n")
}
