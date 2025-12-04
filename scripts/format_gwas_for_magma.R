#!/usr/bin/env Rscript
# Format GWAS summary statistics for MAGMA analysis
# Creates SNP location file and p-value file with optional filtering

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  arg_list <- list()
  i <- 1
  while (i <= length(args)) {
    if (startsWith(args[i], "--")) {
      key <- gsub("^--", "", args[i])
      key <- gsub("-", "_", key)
      if (i + 1 <= length(args) && !startsWith(args[i + 1], "--")) {
        arg_list[[key]] <- args[i + 1]
        i <- i + 2
      } else {
        arg_list[[key]] <- TRUE
        i <- i + 1
      }
    } else {
      i <- i + 1
    }
  }
  return(arg_list)
}

args_parsed <- parse_args(args)

# Required arguments
input_file <- args_parsed$input
output_snploc <- args_parsed$output_snploc
output_pval <- args_parsed$output_pval
output_filtered <- args_parsed$output_filtered
dataset_config_json <- args_parsed$dataset_config

# Check required arguments
if (is.null(input_file) || is.null(output_snploc) || is.null(output_pval) || is.null(output_filtered) || is.null(dataset_config_json)) {
  cat("Usage: Rscript format_gwas_for_magma.R --input FILE --output-snploc FILE --output-pval FILE --output-filtered FILE --dataset-config JSON\n")
  stop("Missing required arguments")
}

# Parse dataset config
dataset_config <- fromJSON(dataset_config_json)

# Extract configuration
sample_size <- dataset_config$samples
columns <- dataset_config$columns
p_threshold <- if (!is.null(dataset_config$p_threshold)) dataset_config$p_threshold else 1.0
meta_missing_threshold <- if (!is.null(dataset_config$meta_missing_threshold)) dataset_config$meta_missing_threshold else Inf

cat("=== GWAS Data Formatting for MAGMA ===\n")
cat("Input file:", input_file, "\n")
cat("Sample size:", sample_size, "\n")
cat("P-value threshold:", p_threshold, "\n")
cat("Meta missing threshold:", meta_missing_threshold, "\n\n")

# Read GWAS data
cat("Reading GWAS summary statistics...\n")
gwas_data <- fread(input_file, 
                   sep = "\t",
                   header = TRUE,
                   stringsAsFactors = FALSE,
                   data.table = TRUE)

cat("Total SNPs read:", nrow(gwas_data), "\n")
cat("Columns found:", paste(colnames(gwas_data), collapse = ", "), "\n\n")

# Check required columns
required_cols <- c(columns$snp, columns$chr, columns$pos, columns$pval)
cat("Looking for required columns:", paste(required_cols, collapse = ", "), "\n")
missing_cols <- setdiff(required_cols, colnames(gwas_data))
if (length(missing_cols) > 0) {
  cat("ERROR: Missing required columns:", paste(missing_cols, collapse = ", "), "\n")
  cat("Available columns:", paste(colnames(gwas_data), collapse = ", "), "\n")
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Rename columns for easier processing
setnames(gwas_data, 
         old = c(columns$snp, columns$chr, columns$pos, columns$pval),
         new = c("SNP", "CHR", "BP", "P"))

# Add Direction column if exists
has_direction <- FALSE
if (!is.null(columns$direction) && columns$direction %in% colnames(gwas_data)) {
  setnames(gwas_data, old = columns$direction, new = "DIRECTION")
  has_direction <- TRUE
  cat("Direction column found:", columns$direction, "\n")
} else {
  cat("No direction column specified\n")
}

# Keep only needed columns - create a copy to avoid reference issues
keep_cols <- c("SNP", "CHR", "BP", "P")
if (has_direction) keep_cols <- c(keep_cols, "DIRECTION")
gwas_clean <- gwas_data[, ..keep_cols]

# DIAGNOSTIC: Check what columns we have
cat("\nDIAGNOSTIC - Columns in gwas_clean after selection:\n")
cat("  Columns:", paste(colnames(gwas_clean), collapse=", "), "\n")

# Convert numeric columns
gwas_clean[, CHR := as.integer(CHR)]
gwas_clean[, BP := as.integer(BP)]
gwas_clean[, P := as.numeric(P)]

# Remove rows with missing essential values
initial_count <- nrow(gwas_clean)
gwas_clean <- gwas_clean[!is.na(SNP) & !is.na(CHR) & !is.na(BP) & !is.na(P)]
cat("SNPs after removing NA values:", nrow(gwas_clean), 
    "(removed:", initial_count - nrow(gwas_clean), ")\n")

# Apply p-value threshold
if (p_threshold < 1.0) {
  before_pval <- nrow(gwas_clean)
  gwas_clean <- gwas_clean[P <= p_threshold]
  cat("SNPs after p-value threshold (P <=", p_threshold, "):", nrow(gwas_clean),
      "(removed:", before_pval - nrow(gwas_clean), ")\n")
}

# Apply Direction filtering if column exists
if (has_direction) {
  before_direction <- nrow(gwas_clean)
  
  # DIAGNOSTIC: Check first few Direction values before conversion
  cat("\nDIAGNOSTIC - Direction values before conversion:\n")
  cat("  Class:", class(gwas_clean$DIRECTION), "\n")
  cat("  First 10 values:", paste(head(gwas_clean$DIRECTION, 10), collapse=", "), "\n")
  cat("  Unique values (first 20):", paste(head(unique(gwas_clean$DIRECTION), 20), collapse=", "), "\n")
  
  # Convert Direction to character and trim whitespace
  gwas_clean[, DIRECTION := as.character(trimws(DIRECTION))]
  
  # Convert all '0' characters in Direction strings to '?' for clarity
  gwas_clean[, DIRECTION := gsub("0", "?", DIRECTION)]
  
  # DIAGNOSTIC: Check after conversion
  cat("\nDIAGNOSTIC - Direction values after conversion (0 -> ?):\n")
  cat("  Class:", class(gwas_clean$DIRECTION), "\n")
  cat("  First 10 values:", paste(head(gwas_clean$DIRECTION, 10), collapse=", "), "\n")
  cat("  Count of strings with '?' (missing):", sum(grepl("\\?", gwas_clean$DIRECTION), na.rm=TRUE), "\n")
  
  # First filter: remove rows where Direction is just "0", "?", empty, or NA
  gwas_clean <- gwas_clean[
    !is.na(DIRECTION) &
    DIRECTION != "" &
    DIRECTION != "0" & 
    DIRECTION != "?"
  ]
  
  cat("\nSNPs after removing Direction='0', '?', empty, or NA:", nrow(gwas_clean), 
      "(removed:", before_direction - nrow(gwas_clean), ")\n")
  
  # Second filter: check signal content for remaining rows
  before_signal_filter <- nrow(gwas_clean)
  
  # Count + and - in direction string
  gwas_clean[, n_plus := nchar(gsub("[^+]", "", DIRECTION))]
  gwas_clean[, n_minus := nchar(gsub("[^-]", "", DIRECTION))]
  gwas_clean[, total_signals := n_plus + n_minus]
  gwas_clean[, direction_length := nchar(DIRECTION)]
  gwas_clean[, n_missing := direction_length - total_signals]
  
  # Filter: exclude if ALL missing (no + or -), OR missing count exceeds threshold
  gwas_clean <- gwas_clean[
    total_signals > 0 & 
    n_missing <= meta_missing_threshold
  ]
  
  # Remove temporary columns
  gwas_clean[, c("n_plus", "n_minus", "total_signals", "direction_length", "n_missing") := NULL]
  
  cat("SNPs after signal filtering (must have +/-, max", meta_missing_threshold, "missing):", 
      nrow(gwas_clean), "(removed:", before_signal_filter - nrow(gwas_clean), ")\n")
}

# Final check
if (nrow(gwas_clean) == 0) {
  stop("No SNPs remaining after filtering!")
}
cat("\n=== Creating Output Files ===\n")

# Create SNP location file (for MAGMA --annotate)
# Format: SNP CHR BP
snploc_data <- gwas_clean[, .(SNP, CHR, BP)]

cat("Writing SNP location file:", output_snploc, "\n")
fwrite(snploc_data, output_snploc, sep = " ", col.names = FALSE, quote = FALSE)
cat("  Lines written:", nrow(snploc_data), "\n")

# Create p-value file (for MAGMA --gene-analysis)
# Format: SNP P
pval_data <- gwas_clean[, .(SNP, P)]

cat("Writing p-value file:", output_pval, "\n")
fwrite(pval_data, output_pval, sep = " ", col.names = FALSE, quote = FALSE)
cat("  Lines written:", nrow(pval_data), "\n")

# Create full filtered table (with all columns, tab-delimited with header)
cat("Writing full filtered table:", output_filtered, "\n")
fwrite(gwas_clean, output_filtered, sep = "\t", col.names = TRUE, quote = FALSE)
cat("  Lines written:", nrow(gwas_clean), "\n")

cat("\n=== Formatting Complete ===\n")
cat("Summary:\n")
cat("  Total SNPs processed:", nrow(gwas_clean), "\n")
cat("  Sample size (N):", sample_size, "\n")
cat("  Output files ready for MAGMA\n")
cat("  Filtered table:", output_filtered, "\n")
