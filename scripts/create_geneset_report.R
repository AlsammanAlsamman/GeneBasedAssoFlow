#!/usr/bin/env Rscript
# Create Gene-Set Enrichment Report
# Parses MAGMA gene-set analysis results and creates Excel report with:
# 1. Geneset_Details sheet: Significant gene-sets with individual gene rows
# 2. Geneset_Loci_Matrix sheet: Gene-sets (rows) by Loci (columns) with gene counts

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
  make_option("--results-base", type="character", help="Base directory for gene-set results"),
  make_option("--target", type="character", help="Target analysis name"),
  make_option("--genesets", type="character", help="Comma-separated list of gene-set database keys"),
  make_option("--geneset-p-threshold", type="numeric", help="P-value threshold for gene-set significance"),
  make_option("--geneset-gene-p-condition", type="numeric", help="P-value threshold for gene filtering (informational)"),
  make_option("--loci-gene-p-threshold", type="numeric", help="P-value threshold for gene-loci matching"),
  make_option("--loci-distance-kb", type="numeric", help="Distance threshold in kb for loci matching"),
  make_option("--output", type="character", help="Output Excel file path")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Validate arguments
required_args <- c("genes-file", "gene-loc", "loci-file", "results-base", "target", 
                   "genesets", "geneset-p-threshold", "geneset-gene-p-condition",
                   "loci-gene-p-threshold", "loci-distance-kb", "output")
missing_args <- required_args[!required_args %in% names(opt) | sapply(opt[required_args], is.null)]
if (length(missing_args) > 0) {
  stop("Missing required arguments: ", paste(missing_args, collapse=", "))
}

cat("========================================\n")
cat("Creating Gene-Set Enrichment Report\n")
cat("========================================\n")
cat("Genes file:", opt$`genes-file`, "\n")
cat("Gene location file:", opt$`gene-loc`, "\n")
cat("Loci file:", opt$`loci-file`, "\n")
cat("Results base:", opt$`results-base`, "\n")
cat("Target:", opt$target, "\n")
cat("Genesets:", opt$genesets, "\n")
cat("Geneset P-value threshold:", opt$`geneset-p-threshold`, "\n")
cat("Loci gene P-value threshold:", opt$`loci-gene-p-threshold`, "\n")
cat("Loci distance threshold:", opt$`loci-distance-kb`, "kb\n")
cat("Output file:", opt$output, "\n")
cat("========================================\n\n")

# Parse geneset list
geneset_list <- strsplit(opt$genesets, ",")[[1]]
if (length(geneset_list) == 0) {
  stop("No gene-sets specified!")
}

# Read gene results with names
cat("Reading gene results...\n")
genes <- fread(opt$`genes-file`, header=TRUE, sep="\t")
cat("Total genes:", nrow(genes), "\n")

# Clean column names
old_names <- names(genes)
new_names <- trimws(old_names)
new_names <- gsub("\\s+", "_", new_names)
setnames(genes, old_names, new_names)

# Find P-value column
p_col_idx <- which(names(genes) == "P")
if (length(p_col_idx) == 0) {
  p_col_idx <- grep("P$", names(genes))
  if (length(p_col_idx) == 0) {
    stop("Cannot find P-value column!")
  }
}
p_col_name <- names(genes)[p_col_idx[1]]
cat("Using P-value column:", p_col_name, "\n")

# Read gene location file
cat("\nReading gene location file...\n")
gene_loc <- fread(opt$`gene-loc`, header=FALSE, 
                  col.names=c("GENE", "CHR", "START", "END", "STRAND", "GENE_NAME"))
cat("Gene locations loaded:", nrow(gene_loc), "\n")

# Read loci file
cat("\nReading loci file...\n")
loci <- fread(opt$`loci-file`, header=TRUE)
cat("Loci loaded:", nrow(loci), "\n")

# Standardize loci column names
if ("locus" %in% names(loci)) setnames(loci, "locus", "Locus", skip_absent=TRUE)
if (!"Locus" %in% names(loci) && any(grepl("locus", names(loci), ignore.case=TRUE))) {
  setnames(loci, grep("locus", names(loci), ignore.case=TRUE)[1], "Locus")
}

# Parse gene-set analysis results
cat("\nParsing gene-set analysis results...\n")
all_genesets <- list()
all_geneset_genes <- list()

database_names <- list(
  msigdb_hallmark_entrez = "MSigDB Hallmark",
  msigdb_hallmark_symbols = "MSigDB Hallmark",
  msigdb_go_bp_entrez = "MSigDB GO Biological Process",
  msigdb_go_bp_symbols = "MSigDB GO Biological Process",
  msigdb_biocarta = "MSigDB BioCarta",
  msigdb_kegg_legacy = "MSigDB KEGG Legacy",
  msigdb_reactome = "MSigDB Reactome",
  msigdb_immunesigdb = "MSigDB ImmuneSigDB",
  msigdb_cell_type = "MSigDB Cell Type"
)

for (db_key in geneset_list) {
  gsa_file <- file.path(opt$`results-base`, paste0(opt$target, "_", db_key, ".gsa.out"))
  
  if (!file.exists(gsa_file)) {
    cat("  WARNING: File not found:", gsa_file, "\n")
    next
  }
  
  cat("  Processing:", db_key, "\n")
  
  # Read entire file
  lines <- readLines(gsa_file)
  
  # Find main gene-set table (before detailed sections)
  header_idx <- grep("^VARIABLE\\s+TYPE\\s+NGENES", lines)
  if (length(header_idx) == 0) {
    cat("    WARNING: No gene-set table found\n")
    next
  }
  
  # Find where detailed sections start (marked by # _SET1_, # _SET2_, etc.)
  detail_start_idx <- grep("^# _SET[0-9]+_", lines)
  if (length(detail_start_idx) > 0) {
    table_end_idx <- detail_start_idx[1] - 1
  } else {
    table_end_idx <- length(lines)
  }
  
  # Extract gene-set summary table
  table_lines <- lines[(header_idx+1):table_end_idx]
  table_lines <- table_lines[!grepl("^#", table_lines) & nzchar(trimws(table_lines))]
  
  if (length(table_lines) == 0) {
    cat("    No gene-sets found\n")
    next
  }
  
  # Parse table
  con <- textConnection(table_lines)
  geneset_table <- tryCatch({
    read.table(con, header=FALSE, sep="", fill=TRUE, stringsAsFactors=FALSE,
               col.names=c("VARIABLE", "TYPE", "NGENES", "BETA", "BETA_STD", "SE", "P", "FULL_NAME"))
  }, error = function(e) {
    NULL
  })
  close(con)
  
  if (is.null(geneset_table) || nrow(geneset_table) == 0) {
    cat("    Failed to parse gene-set table\n")
    next
  }
  
  geneset_table <- as.data.table(geneset_table)
  geneset_table[, Database := ifelse(db_key %in% names(database_names), 
                                      database_names[[db_key]], db_key)]
  geneset_table[, Database_Key := db_key]
  
  # Filter by P-value threshold
  sig_genesets <- geneset_table[P <= opt$`geneset-p-threshold`]
  
  if (nrow(sig_genesets) > 0) {
    cat("    Significant gene-sets:", nrow(sig_genesets), "\n")
    all_genesets[[db_key]] <- sig_genesets
    
    # Parse detailed gene information for significant gene-sets
    for (i in 1:nrow(sig_genesets)) {
      geneset_var <- sig_genesets$VARIABLE[i]
      
      # Find this gene-set's detail section
      set_pattern <- paste0("^# _SET[0-9]+_\\s+VARIABLE = ", gsub("\\.", "\\\\.", geneset_var))
      set_header_idx <- grep(set_pattern, lines)
      
      if (length(set_header_idx) == 0) next
      
      # Find gene table header
      gene_header_pattern <- paste0("^_SET[0-9]+_\\s+GENE\\s+CHR")
      gene_header_idx <- grep(gene_header_pattern, lines)
      gene_header_idx <- gene_header_idx[gene_header_idx > set_header_idx[1]]
      
      if (length(gene_header_idx) == 0) next
      
      # Find end of gene table (next # _SET or end of file)
      next_set_idx <- grep("^# _SET[0-9]+_", lines)
      next_set_idx <- next_set_idx[next_set_idx > gene_header_idx[1]]
      
      if (length(next_set_idx) > 0) {
        gene_end_idx <- next_set_idx[1] - 1
      } else {
        gene_end_idx <- length(lines)
      }
      
      # Extract gene lines
      gene_lines <- lines[(gene_header_idx[1]+1):gene_end_idx]
      gene_lines <- gene_lines[grepl("^_SET[0-9]+_", gene_lines)]
      
      if (length(gene_lines) == 0) next
      
      # Parse gene table
      gene_con <- textConnection(gene_lines)
      gene_table <- tryCatch({
        read.table(gene_con, header=FALSE, sep="", fill=TRUE, stringsAsFactors=FALSE)
      }, error = function(e) {
        NULL
      })
      close(gene_con)
      
      if (!is.null(gene_table) && ncol(gene_table) >= 11) {
        colnames(gene_table) <- c("SET_ID", "GENE", "CHR", "START", "STOP", "NSNPS", 
                                   "NPARAM", "N", "ZSTAT", "P", "ZFITTED_BASE", "ZRESID_BASE")
        gene_table <- as.data.table(gene_table)
        gene_table[, Geneset := geneset_var]
        gene_table[, Database := sig_genesets$Database[i]]
        gene_table[, Geneset_P := sig_genesets$P[i]]
        gene_table[, Geneset_BETA := sig_genesets$BETA[i]]
        gene_table[, Geneset_NGENES := sig_genesets$NGENES[i]]
        
        all_geneset_genes[[paste0(db_key, "_", geneset_var)]] <- gene_table
      }
    }
  }
}

if (length(all_genesets) == 0) {
  cat("\nNo significant gene-sets found at P <=", opt$`geneset-p-threshold`, "\n")
  cat("Creating empty output file.\n")
  
  wb <- wb_workbook()
  wb$add_worksheet("Geneset_Details")
  wb$add_data(sheet=1, x=data.frame(Message="No significant gene-sets found at specified threshold"))
  wb$save(opt$output)
  
  cat("Empty Excel file created.\n")
  quit(save="no", status=0)
}

# Combine all significant gene-sets
cat("\nCombining results...\n")
combined_genesets <- rbindlist(all_genesets, fill=TRUE)
cat("Total significant gene-sets:", nrow(combined_genesets), "\n")

# Combine all genes from significant gene-sets
if (length(all_geneset_genes) > 0) {
  combined_genes <- rbindlist(all_geneset_genes, fill=TRUE)
  cat("Total genes in significant gene-sets:", nrow(combined_genes), "\n")
  
  # Add gene names
  combined_genes <- merge(combined_genes, 
                          gene_loc[, .(GENE, GENE_NAME)], 
                          by="GENE", all.x=TRUE)
  
  # Sort by Geneset and gene P-value
  combined_genes <- combined_genes[order(Geneset_P, Geneset, P)]
  
  # Create Sheet 1: Geneset_Details with merged cells
  cat("\nCreating Geneset_Details sheet...\n")
  
  # Reorder columns for better readability
  geneset_cols <- c("Database", "Geneset", "Geneset_P", "Geneset_BETA", "Geneset_NGENES")
  gene_cols <- c("GENE", "GENE_NAME", "CHR", "START", "STOP", "P", "ZSTAT", "NSNPS")
  other_cols <- setdiff(names(combined_genes), c(geneset_cols, gene_cols))
  
  setcolorder(combined_genes, c(geneset_cols, gene_cols, other_cols))
  
  wb <- wb_workbook()
  wb$add_worksheet("Geneset_Details")
  wb$add_data(sheet=1, x=combined_genes, start_col=1, start_row=1)
  
  # Format P-value columns as scientific notation
  geneset_p_idx <- which(names(combined_genes) == "Geneset_P")
  if (length(geneset_p_idx) > 0) {
    geneset_p_letter <- LETTERS[geneset_p_idx]
    wb$add_numfmt(sheet=1, dims=paste0(geneset_p_letter, "2:", geneset_p_letter, nrow(combined_genes)+1), 
                  numfmt="0.00E+00")
  }
  
  gene_p_idx <- which(names(combined_genes) == "P")
  if (length(gene_p_idx) > 0) {
    gene_p_letter <- LETTERS[gene_p_idx]
    wb$add_numfmt(sheet=1, dims=paste0(gene_p_letter, "2:", gene_p_letter, nrow(combined_genes)+1), 
                  numfmt="0.00E+00")
  }
  
  # Style header
  header_range <- paste0("A1:", LETTERS[ncol(combined_genes)], "1")
  wb$add_fill(sheet=1, dims=header_range, color=wb_color(hex="4472C4"))
  wb$add_font(sheet=1, dims=header_range, color=wb_color(hex="FFFFFF"), bold=TRUE)
  wb$add_cell_style(sheet=1, dims=header_range, horizontal="center", vertical="center", textRotation=90)
  
  # Apply center and middle alignment to all data cells
  data_range <- paste0("A2:", LETTERS[ncol(combined_genes)], nrow(combined_genes)+1)
  wb$add_cell_style(sheet=1, dims=data_range, horizontal="center", vertical="center")
  
  # Merge cells for gene-set information
  cat("Applying merged cells for gene-sets...\n")
  unique_genesets <- unique(combined_genes$Geneset)
  merge_col_count <- length(geneset_cols)
  
  current_row <- 2  # Start from row 2 (after header)
  for (geneset in unique_genesets) {
    geneset_rows <- which(combined_genes$Geneset == geneset)
    n_genes <- length(geneset_rows)
    
    if (n_genes > 1) {
      end_row <- current_row + n_genes - 1
      # Merge geneset info columns
      for (col_idx in 1:merge_col_count) {
        wb$merge_cells(sheet=1, rows=current_row:end_row, cols=col_idx)
      }
    }
    
    current_row <- current_row + n_genes
  }
  
  # Auto-size columns
  wb$set_col_widths(sheet=1, cols=1:ncol(combined_genes), widths="auto")
  
} else {
  cat("\nNo gene-level details found for significant gene-sets\n")
  
  wb <- wb_workbook()
  wb$add_worksheet("Geneset_Details")
  wb$add_data(sheet=1, x=combined_genesets)
  
  # Style header
  header_range <- paste0("A1:", LETTERS[ncol(combined_genesets)], "1")
  wb$add_fill(sheet=1, dims=header_range, color=wb_color(hex="4472C4"))
  wb$add_font(sheet=1, dims=header_range, color=wb_color(hex="FFFFFF"), bold=TRUE)
}

# Create Sheet 2: Geneset_Loci_Matrix
cat("\nCreating Geneset_Loci_Matrix sheet...\n")

# Merge genes with positions
genes_with_pos <- merge(genes, 
                        gene_loc[, .(GENE, CHR_loc=CHR, START_loc=START, END_loc=END)], 
                        by="GENE", all.x=TRUE)
genes_with_pos[, ':='(CHR_final = CHR_loc, START_final = START_loc, END_final = END_loc)]
genes_with_pos[, ':='(CHR_loc = NULL, START_loc = NULL, END_loc = NULL)]

# Filter genes by threshold for loci matching
genes_sig <- genes_with_pos[get(p_col_name) <= opt$`loci-gene-p-threshold`]
cat("Significant genes for loci matching (P <=", opt$`loci-gene-p-threshold`, "):", nrow(genes_sig), "\n")

# Calculate distance buffer
distance_bp <- opt$`loci-distance-kb` * 1000

# Need to read GMT files to get gene-set memberships since MAGMA doesn't output gene details
# Build database_files mapping
database_files <- list()
config_lines <- readLines("configs/analysis.yml")
in_databases <- FALSE
for (line in config_lines) {
  if (grepl("^databases:", line)) {
    in_databases <- TRUE
    next
  }
  if (in_databases && grepl("^\\s+\\w+:", line)) {
    # Parse: "  key: /path/to/file.gmt"
    parts <- strsplit(trimws(line), ":\\s*", perl=TRUE)[[1]]
    if (length(parts) == 2) {
      key <- trimws(parts[1])
      path <- trimws(gsub('"', '', parts[2]))
      database_files[[key]] <- path
    }
  }
  if (in_databases && !grepl("^\\s+", line) && grepl("^\\w+:", line)) {
    in_databases <- FALSE
  }
}

# Function to read GMT file
read_gmt <- function(gmt_file) {
  if (!file.exists(gmt_file)) {
    cat("WARNING: GMT file not found:", gmt_file, "\n")
    return(NULL)
  }
  
  lines <- readLines(gmt_file)
  geneset_list <- list()
  
  for (line in lines) {
    parts <- strsplit(line, "\t")[[1]]
    if (length(parts) >= 3) {
      geneset_name <- parts[1]
      genes <- parts[3:length(parts)]
      genes <- genes[nchar(genes) > 0]
      geneset_list[[geneset_name]] <- genes
    }
  }
  
  return(geneset_list)
}

# Read GMT files for databases we're analyzing
all_gmt_data <- list()
for (db_key in geneset_list) {
  if (db_key %in% names(database_files)) {
    gmt_file <- database_files[[db_key]]
    cat("Reading GMT file for", db_key, ":", gmt_file, "\n")
    gmt_data <- read_gmt(gmt_file)
    if (!is.null(gmt_data)) {
      all_gmt_data[[db_key]] <- gmt_data
    }
  }
}

# Create matrix: rows = gene-sets, columns = loci
matrix_list <- list()

for (geneset in unique(combined_genesets$VARIABLE)) {
  # Find which database this geneset belongs to
  db_key <- combined_genesets[VARIABLE == geneset, Database_Key][1]
  
  # Get genes in this gene-set from GMT file
  geneset_genes <- NULL
  if (db_key %in% names(all_gmt_data)) {
    # Try exact match first
    if (geneset %in% names(all_gmt_data[[db_key]])) {
      geneset_genes <- all_gmt_data[[db_key]][[geneset]]
    } else {
      # Try matching without abbreviations (e.g., "HALLMARK_ADIPOGENESIS" should match)
      for (gmt_name in names(all_gmt_data[[db_key]])) {
        if (grepl(geneset, gmt_name, fixed=TRUE) || grepl(gmt_name, geneset, fixed=TRUE)) {
          geneset_genes <- all_gmt_data[[db_key]][[gmt_name]]
          break
        }
      }
    }
  }
  
  if (is.null(geneset_genes) || length(geneset_genes) == 0) {
    cat("  WARNING: No genes found for", geneset, "\n")
    # Create row with zeros
    row_data <- data.table(
      Database = combined_genesets[VARIABLE == geneset, Database][1],
      Database_Key = combined_genesets[VARIABLE == geneset, Database_Key][1],
      Geneset = geneset,
      Geneset_P = combined_genesets[VARIABLE == geneset, P][1],
      Total_Genes = combined_genesets[VARIABLE == geneset, NGENES][1]
    )
    for (locus in loci$Locus) {
      row_data[, (locus) := 0]
    }
    matrix_list[[geneset]] <- row_data
    next
  }
  
  # Convert to numeric Entrez IDs if needed
  geneset_genes_numeric <- suppressWarnings(as.numeric(geneset_genes))
  geneset_genes_numeric <- geneset_genes_numeric[!is.na(geneset_genes_numeric)]
  
  # Get gene info for these genes that pass threshold
  geneset_gene_info <- genes_sig[GENE %in% geneset_genes_numeric]
  
  if (nrow(geneset_gene_info) == 0) {
    # No genes pass threshold - create row with zeros
    row_data <- data.table(
      Database = combined_genesets[VARIABLE == geneset, Database][1],
      Database_Key = combined_genesets[VARIABLE == geneset, Database_Key][1],
      Geneset = geneset,
      Geneset_P = combined_genesets[VARIABLE == geneset, P][1],
      Total_Genes = combined_genesets[VARIABLE == geneset, NGENES][1]
    )
    
    # Add zero counts for each locus
    for (locus in loci$Locus) {
      row_data[, (locus) := 0]
    }
    
    matrix_list[[geneset]] <- row_data
    next
  }  # Count genes in each locus
  row_data <- data.table(
    Database = combined_genesets[VARIABLE == geneset, Database][1],
    Database_Key = combined_genesets[VARIABLE == geneset, Database_Key][1],
    Geneset = geneset,
    Geneset_P = combined_genesets[VARIABLE == geneset, P][1],
    Total_Genes = combined_genesets[VARIABLE == geneset, NGENES][1]
  )
  
  for (i in 1:nrow(loci)) {
    locus_name <- loci$Locus[i]
    locus_chr <- loci$chr[i]
    locus_start <- loci$start[i] - distance_bp
    locus_end <- loci$end[i] + distance_bp
    
    # Count genes in this locus
    n_genes_in_locus <- nrow(geneset_gene_info[
      CHR_final == locus_chr & 
      !is.na(START_final) & !is.na(END_final) &
      ((START_final >= locus_start & START_final <= locus_end) |
       (END_final >= locus_start & END_final <= locus_end) |
       (START_final <= locus_start & END_final >= locus_end))
    ])
    
    row_data[, (locus_name) := n_genes_in_locus]
  }
  
  matrix_list[[geneset]] <- row_data
}

if (length(matrix_list) > 0) {
  geneset_loci_matrix <- rbindlist(matrix_list, fill=TRUE)
  
  # Add column for genes passing threshold
  cat("\nCalculating genes passing threshold (P <", opt$`geneset-gene-p-condition`, ") per gene-set...\n")
  geneset_loci_matrix[, N_Genes_Pass_Threshold := 0]
  
  for (i in 1:nrow(geneset_loci_matrix)) {
    geneset_name <- geneset_loci_matrix$Geneset[i]
    db_key <- geneset_loci_matrix$Database_Key[i]
    
    # Get genes in this gene-set from GMT data
    geneset_genes <- NULL
    if (db_key %in% names(all_gmt_data)) {
      if (geneset_name %in% names(all_gmt_data[[db_key]])) {
        geneset_genes <- all_gmt_data[[db_key]][[geneset_name]]
      }
    }
    
    if (!is.null(geneset_genes) && length(geneset_genes) > 0) {
      # Convert to numeric Entrez IDs
      geneset_genes_numeric <- suppressWarnings(as.numeric(geneset_genes))
      geneset_genes_numeric <- geneset_genes_numeric[!is.na(geneset_genes_numeric)]
      
      # Count genes passing threshold
      n_pass <- nrow(genes[GENE %in% geneset_genes_numeric & get(p_col_name) < opt$`geneset-gene-p-condition`])
      geneset_loci_matrix[i, N_Genes_Pass_Threshold := n_pass]
    }
  }
  
  # Reorder columns to put N_Genes_Pass_Threshold after Total_Genes
  total_genes_idx <- which(names(geneset_loci_matrix) == "Total_Genes")
  threshold_idx <- which(names(geneset_loci_matrix) == "N_Genes_Pass_Threshold")
  if (length(total_genes_idx) > 0 && length(threshold_idx) > 0) {
    other_cols <- setdiff(names(geneset_loci_matrix), c("N_Genes_Pass_Threshold"))
    insert_pos <- total_genes_idx + 1
    new_order <- c(other_cols[1:total_genes_idx], "N_Genes_Pass_Threshold", other_cols[(total_genes_idx+1):length(other_cols)])
    setcolorder(geneset_loci_matrix, new_order)
  }
  
  # Sort by Geneset_P
  geneset_loci_matrix <- geneset_loci_matrix[order(Geneset_P)]
  
  # Replace NA with 0
  for (col in names(geneset_loci_matrix)) {
    if (is.numeric(geneset_loci_matrix[[col]])) {
      geneset_loci_matrix[is.na(get(col)), (col) := 0]
    }
  }
  
  cat("Gene-set by Loci matrix: ", nrow(geneset_loci_matrix), "gene-sets x", 
      nrow(loci), "loci\n")
  
  # Add to workbook
  wb$add_worksheet("Geneset_Loci_Matrix")
  wb$add_data(sheet="Geneset_Loci_Matrix", x=geneset_loci_matrix, start_col=1, start_row=1)
  
  # Format Geneset_P as scientific notation
  geneset_p_idx <- which(names(geneset_loci_matrix) == "Geneset_P")
  if (length(geneset_p_idx) > 0) {
    geneset_p_letter <- LETTERS[geneset_p_idx]
    wb$add_numfmt(sheet="Geneset_Loci_Matrix", 
                  dims=paste0(geneset_p_letter, "2:", geneset_p_letter, nrow(geneset_loci_matrix)+1), 
                  numfmt="0.00E+00")
  }
  
  # Style header
  header_range <- paste0("A1:", LETTERS[ncol(geneset_loci_matrix)], "1")
  wb$add_fill(sheet="Geneset_Loci_Matrix", dims=header_range, color=wb_color(hex="4472C4"))
  wb$add_font(sheet="Geneset_Loci_Matrix", dims=header_range, color=wb_color(hex="FFFFFF"), bold=TRUE)
  wb$add_cell_style(sheet="Geneset_Loci_Matrix", dims=header_range, 
                    horizontal="center", vertical="center", textRotation=90)
  
  # Apply center and middle alignment to all data cells
  data_range <- paste0("A2:", LETTERS[ncol(geneset_loci_matrix)], nrow(geneset_loci_matrix)+1)
  wb$add_cell_style(sheet="Geneset_Loci_Matrix", dims=data_range, 
                    horizontal="center", vertical="center")
  
  # Auto-size columns
  wb$set_col_widths(sheet="Geneset_Loci_Matrix", cols=1:ncol(geneset_loci_matrix), widths="auto")
} else {
  cat("No gene-sets found to create matrix.\n")
}

# Sheet 3: Analysis Parameters
cat("\nCreating Analysis_Parameters sheet...\n")
params_data <- data.table(
  Parameter = c(
    "Target Analysis",
    "Gene-set Databases",
    "Gene-set P-value Threshold",
    "Gene P-value Threshold (for gene filtering)",
    "Loci Gene P-value Threshold",
    "Loci Distance Threshold (kb)",
    "Number of Significant Gene-sets",
    "Number of Loci",
    "Total Genes Analyzed"
  ),
  Value = c(
    opt$target,
    paste(geneset_list, collapse=", "),
    format(opt$`geneset-p-threshold`, scientific=TRUE),
    format(opt$`geneset-gene-p-condition`, scientific=TRUE),
    format(opt$`loci-gene-p-threshold`, scientific=TRUE),
    opt$`loci-distance-kb`,
    nrow(combined_genesets),
    nrow(loci),
    nrow(genes)
  )
)

wb$add_worksheet("Analysis_Parameters", grid_lines=FALSE)
wb$add_data(sheet="Analysis_Parameters", x=params_data, start_col=1, start_row=1)

# Style header
wb$add_fill(sheet="Analysis_Parameters", dims="A1:B1", color=wb_color(hex="4472C4"))
wb$add_font(sheet="Analysis_Parameters", dims="A1:B1", color=wb_color(hex="FFFFFF"), bold=TRUE)
wb$add_cell_style(sheet="Analysis_Parameters", dims="A1:B1", horizontal="center", vertical="center")

# Style parameter column
wb$add_font(sheet="Analysis_Parameters", dims=paste0("A2:A", nrow(params_data)+1), bold=TRUE)

# Auto-size columns
wb$set_col_widths(sheet="Analysis_Parameters", cols=1:2, widths="auto")

# Sheet 4: Significant Genes (genes passing threshold)
cat("\nCreating Significant_Genes sheet...\n")

# Get genes passing threshold - use genes_with_pos which already has location info
genes_pass <- genes_with_pos[get(p_col_name) < opt$`geneset-gene-p-condition`]
cat("Genes passing threshold (P <", opt$`geneset-gene-p-condition`, "):", nrow(genes_pass), "\n")

if (nrow(genes_pass) > 0) {
  # Rename position columns if they exist as CHR_final, START_final, END_final
  if ("CHR_final" %in% names(genes_pass)) {
    setnames(genes_pass, "CHR_final", "CHR", skip_absent=TRUE)
  }
  if ("START_final" %in% names(genes_pass)) {
    setnames(genes_pass, "START_final", "START", skip_absent=TRUE)
  }
  if ("END_final" %in% names(genes_pass)) {
    setnames(genes_pass, "END_final", "END", skip_absent=TRUE)
  }
  
  # Add gene names from gene_loc if not already present
  if (!"GENE_NAME" %in% names(genes_pass)) {
    genes_pass <- merge(genes_pass, 
                        gene_loc[, .(GENE, GENE_NAME)], 
                        by="GENE", all.x=TRUE)
  }
  
  # Fill missing GENE_NAME with GENE ID
  if ("GENE_NAME" %in% names(genes_pass)) {
    genes_pass[is.na(GENE_NAME), GENE_NAME := as.character(GENE)]
  } else {
    genes_pass[, GENE_NAME := as.character(GENE)]
  }
  
  # For each gene, find which gene-sets contain it
  genes_pass[, Genesets := ""]
  genes_pass[, N_Genesets := 0]
  
  for (i in 1:nrow(genes_pass)) {
    gene_id <- genes_pass$GENE[i]
    gene_id_str <- as.character(gene_id)
    
    # Find all gene-sets containing this gene
    matching_genesets <- c()
    
    for (db_key in names(all_gmt_data)) {
      for (geneset_name in names(all_gmt_data[[db_key]])) {
        geneset_genes <- all_gmt_data[[db_key]][[geneset_name]]
        if (gene_id_str %in% geneset_genes) {
          # Check if this gene-set is significant
          if (geneset_name %in% combined_genesets$VARIABLE) {
            matching_genesets <- c(matching_genesets, geneset_name)
          }
        }
      }
    }
    
    if (length(matching_genesets) > 0) {
      genes_pass[i, Genesets := paste(matching_genesets, collapse=", ")]
      genes_pass[i, N_Genesets := length(matching_genesets)]
    }
  }
  
  # Create output table - build dynamically based on available columns
  sig_genes_table <- data.table(
    Gene_ID = genes_pass$GENE,
    Gene_Name = if("GENE_NAME" %in% names(genes_pass)) genes_pass$GENE_NAME else as.character(genes_pass$GENE),
    Chr = if("CHR" %in% names(genes_pass)) genes_pass$CHR else NA,
    Start = if("START" %in% names(genes_pass)) genes_pass$START else NA,
    End = if("END" %in% names(genes_pass)) genes_pass$END else NA,
    P_Value = genes_pass[[p_col_name]],
    N_Significant_Genesets = genes_pass$N_Genesets,
    Significant_Genesets = genes_pass$Genesets
  )
  
  # Sort by P-value
  sig_genes_table <- sig_genes_table[order(P_Value)]
  
  # Add to workbook
  wb$add_worksheet("Significant_Genes")
  wb$add_data(sheet="Significant_Genes", x=sig_genes_table, start_col=1, start_row=1)
  
  # Format P-value column as scientific notation
  p_col_idx <- which(names(sig_genes_table) == "P_Value")
  if (length(p_col_idx) > 0) {
    p_col_letter <- LETTERS[p_col_idx]
    wb$add_numfmt(sheet="Significant_Genes", 
                  dims=paste0(p_col_letter, "2:", p_col_letter, nrow(sig_genes_table)+1), 
                  numfmt="0.00E+00")
  }
  
  # Style header
  header_range <- paste0("A1:", LETTERS[ncol(sig_genes_table)], "1")
  wb$add_fill(sheet="Significant_Genes", dims=header_range, color=wb_color(hex="4472C4"))
  wb$add_font(sheet="Significant_Genes", dims=header_range, color=wb_color(hex="FFFFFF"), bold=TRUE)
  wb$add_cell_style(sheet="Significant_Genes", dims=header_range, horizontal="center", vertical="center")
  
  # Center align all columns except Genesets
  for (col_idx in 1:(ncol(sig_genes_table)-1)) {
    col_letter <- LETTERS[col_idx]
    wb$add_cell_style(sheet="Significant_Genes", 
                      dims=paste0(col_letter, "2:", col_letter, nrow(sig_genes_table)+1),
                      horizontal="center", vertical="center")
  }
  
  # Auto-size columns
  wb$set_col_widths(sheet="Significant_Genes", cols=1:ncol(sig_genes_table), widths="auto")
  
  cat("Significant genes table created:", nrow(sig_genes_table), "genes\n")
} else {
  cat("No genes pass the threshold. Skipping Significant_Genes sheet.\n")
}

# Save workbook
wb$save(opt$output)

cat("\n========================================\n")
cat("Gene-set enrichment report created successfully!\n")
cat("Output:", opt$output, "\n")
cat("Significant gene-sets:", nrow(combined_genesets), "\n")
if (exists("combined_genes")) {
  cat("Total genes in gene-sets:", nrow(combined_genes), "\n")
}
if (exists("geneset_loci_matrix")) {
  cat("Matrix dimensions:", nrow(geneset_loci_matrix), "gene-sets x", nrow(loci), "loci\n")
}
cat("========================================\n")
