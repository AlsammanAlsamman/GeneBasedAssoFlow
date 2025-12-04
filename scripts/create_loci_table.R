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
loci_without_genes <- c()

for (i in 1:nrow(loci)) {
  locus_name <- loci$Locus[i]
  locus_chr <- loci$chr[i]
  locus_start <- loci$start[i] - distance_bp
  locus_end <- loci$end[i] + distance_bp
  
  # Find genes overlapping this locus region (with buffer) at original threshold
  matched_genes <- genes_with_pos[
    CHR_final == locus_chr & 
    !is.na(START_final) & !is.na(END_final) &
    ((START_final >= locus_start & START_final <= locus_end) |
     (END_final >= locus_start & END_final <= locus_end) |
     (START_final <= locus_start & END_final >= locus_end))
  ]
  
  # If no genes found, try relaxing threshold progressively
  relaxed_threshold <- NA_real_
  if (nrow(matched_genes) == 0) {
    cat("  Locus", locus_name, "has no genes at P <=", opt$`p-threshold`, "- trying relaxed thresholds...\n")
    
    # Try increasing threshold by factors of 10 up to 0.05
    max_threshold <- 0.05
    current_threshold <- opt$`p-threshold`
    while (nrow(matched_genes) == 0 && current_threshold < max_threshold) {
      current_threshold <- current_threshold * 10
      
      # Get genes at relaxed threshold
      genes_relaxed <- genes[get(p_col_name) <= current_threshold]
      genes_relaxed_pos <- merge(genes_relaxed, 
                                  gene_loc[, .(GENE, CHR_loc=CHR, START_loc=START, END_loc=END)], 
                                  by="GENE", all.x=TRUE)
      genes_relaxed_pos[, ':='(CHR_final = CHR_loc, START_final = START_loc, END_final = END_loc)]
      genes_relaxed_pos[, ':='(CHR_loc = NULL, START_loc = NULL, END_loc = NULL)]
      
      matched_genes <- genes_relaxed_pos[
        CHR_final == locus_chr & 
        !is.na(START_final) & !is.na(END_final) &
        ((START_final >= locus_start & START_final <= locus_end) |
         (END_final >= locus_start & END_final <= locus_end) |
         (START_final <= locus_start & END_final >= locus_end))
      ]
      
      if (nrow(matched_genes) > 0) {
        cat("    Found", nrow(matched_genes), "gene(s) at relaxed P <=", current_threshold, "\n")
        relaxed_threshold <- current_threshold
        matched_genes[, Relaxed_Threshold := relaxed_threshold]
      }
    }
  } else {
    # Genes found at original threshold - set Relaxed_Threshold to NA
    matched_genes[, Relaxed_Threshold := NA_real_]
  }
  
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
  } else {
    # No genes found even with relaxed threshold - add locus info only
    cat("  Locus", locus_name, "has no genes even at relaxed threshold - adding locus info only\n")
    loci_without_genes <- c(loci_without_genes, locus_name)
    
    # Create empty row with locus info
    empty_row <- data.table(
      Locus = locus_name,
      Locus_Chr = locus_chr,
      Locus_Start = loci$start[i],
      Locus_End = loci$end[i],
      GENE = NA_character_,
      GENE_NAME = "No significant genes",
      CHR_final = NA_integer_,
      START_final = NA_integer_,
      END_final = NA_integer_
    )
    
    # Add p-value column
    empty_row[, (p_col_name) := NA_real_]
    
    # Add any additional columns from loci file
    loci_cols <- setdiff(names(loci), c("Locus", "chr", "start", "end"))
    for (col in loci_cols) {
      empty_row[, paste0("Locus_", col) := loci[[col]][i]]
    }
    
    matched_list[[i]] <- empty_row
  }
}

# Combine all matches
if (length(matched_list) > 0) {
  result <- rbindlist(matched_list, fill=TRUE)
  
  cat("\nTotal gene-locus matches:", nrow(result), "\n")
  
  if (length(loci_without_genes) > 0) {
    cat("Loci without significant genes:", length(loci_without_genes), "\n")
    cat("  ", paste(loci_without_genes, collapse=", "), "\n")
  }
  
  # Sort by locus and P-value
  result <- result[order(Locus_Chr, Locus_Start, get(p_col_name))]
  
  # Identify representative gene (most significant) for each locus
  cat("Identifying representative genes for each locus...\n")
  result[, Representative_Gene := "", by=Locus]
  result[, Representative_Gene_P := NA_real_, by=Locus]
  
  for (locus in unique(result$Locus)) {
    locus_genes <- result[Locus == locus & !is.na(get(p_col_name))]
    if (nrow(locus_genes) > 0) {
      min_p_idx <- which.min(locus_genes[[p_col_name]])
      rep_gene <- locus_genes$GENE_NAME[min_p_idx]
      rep_gene_p <- locus_genes[[p_col_name]][min_p_idx]
      result[Locus == locus, Representative_Gene := rep_gene]
      result[Locus == locus, Representative_Gene_P := rep_gene_p]
    }
  }
  
  # Reorder columns: Locus info, Representative gene info, then individual gene info
  locus_info_cols <- c("Locus", "Locus_Chr", "Locus_Start", "Locus_End")
  additional_locus_cols <- grep("^Locus_", names(result), value=TRUE)
  additional_locus_cols <- setdiff(additional_locus_cols, locus_info_cols)
  
  representative_cols <- c("Representative_Gene", "Representative_Gene_P")
  
  # Add Relaxed_Threshold column if it doesn't exist and set to original threshold when not relaxed
  if (!"Relaxed_Threshold" %in% names(result)) {
    result[, Relaxed_Threshold := opt$`p-threshold`]
  } else {
    # Replace NA values with original threshold
    result[is.na(Relaxed_Threshold), Relaxed_Threshold := opt$`p-threshold`]
  }
  
  gene_info_cols <- c("GENE", "GENE_NAME", "CHR_final", "START_final", "END_final", p_col_name, "Relaxed_Threshold")
  other_cols <- setdiff(names(result), c(locus_info_cols, additional_locus_cols, representative_cols, gene_info_cols))
  
  setcolorder(result, c(locus_info_cols, additional_locus_cols, representative_cols, gene_info_cols, other_cols))
  
  # Create Excel workbook
  cat("\nCreating Excel file...\n")
  wb <- wb_workbook()
  wb$add_worksheet("Loci_Genes")
  
  # Add data
  wb$add_data(sheet=1, x=result, start_col=1, start_row=1)
  
  # Format P-value columns as scientific notation
  p_col_idx <- which(names(result) == p_col_name)
  if (length(p_col_idx) > 0) {
    p_col_letter <- LETTERS[p_col_idx]
    wb$add_numfmt(sheet=1, dims=paste0(p_col_letter, "2:", p_col_letter, nrow(result)+1), numfmt="0.00E+00")
  }
  
  rep_p_col_idx <- which(names(result) == "Representative_Gene_P")
  if (length(rep_p_col_idx) > 0) {
    rep_p_col_letter <- LETTERS[rep_p_col_idx]
    wb$add_numfmt(sheet=1, dims=paste0(rep_p_col_letter, "2:", rep_p_col_letter, nrow(result)+1), numfmt="0.00E+00")
  }
  
  # Style header
  header_range <- paste0("A1:", LETTERS[ncol(result)], "1")
  wb$add_fill(sheet=1, dims=header_range, color=wb_color(hex="4472C4"))
  wb$add_font(sheet=1, dims=header_range, color=wb_color(hex="FFFFFF"), bold=TRUE)
  wb$add_cell_style(sheet=1, dims=header_range, horizontal="center", vertical="center", textRotation=90)
  
  # Merge cells for loci information and representative gene info
  cat("Applying merged cells for loci regions...\n")
  unique_loci <- unique(result$Locus)
  
  # Find how many columns to merge (locus info + representative gene info)
  locus_col_count <- length(c(locus_info_cols, additional_locus_cols))
  representative_col_count <- length(representative_cols)
  merge_col_count <- locus_col_count + representative_col_count
  
  current_row <- 2  # Start from row 2 (after header)
  for (locus in unique_loci) {
    locus_rows <- which(result$Locus == locus)
    n_genes <- length(locus_rows)
    
    if (n_genes > 1) {
      end_row <- current_row + n_genes - 1
      # Merge all locus info columns and representative gene columns
      for (col_idx in 1:merge_col_count) {
        wb$merge_cells(sheet=1, rows=current_row:end_row, cols=col_idx)
      }
    }
    
    current_row <- current_row + n_genes
  }
  
  # Apply color formatting to Representative_Gene column cells based on P-value thresholds
  cat("Applying color formatting to Representative_Gene column based on P-value thresholds...\n")
  rep_gene_col_idx <- which(names(result) == "Representative_Gene")
  
  if (length(rep_gene_col_idx) > 0) {
    rep_gene_col_letter <- LETTERS[rep_gene_col_idx]
    
    for (i in 1:nrow(result)) {
      row_num <- i + 1  # +1 for header
      p_val <- result$Representative_Gene_P[i]
      
      if (!is.na(p_val)) {
        # Determine color based on P-value
        if (p_val < 5E-8) {
          color_hex <- "FFC0CB"  # Pink
        } else if (p_val < 5E-5) {
          color_hex <- "90EE90"  # Green
        } else if (p_val < 5E-3) {
          color_hex <- "FFFF99"  # Yellow
        } else {
          color_hex <- "FFFFFF"  # White
        }
        
        # Apply fill to Representative_Gene cell only
        cell_range <- paste0(rep_gene_col_letter, row_num)
        wb$add_fill(sheet=1, dims=cell_range, color=wb_color(hex=color_hex))
      }
    }
  }
  
  # Apply center and middle alignment to all data cells
  cat("Applying cell alignment...\n")
  data_range <- paste0("A2:", LETTERS[ncol(result)], nrow(result)+1)
  wb$add_cell_style(sheet=1, dims=data_range, horizontal="center", vertical="center")
  
  # Auto-size columns
  wb$set_col_widths(sheet=1, cols=1:ncol(result), widths="auto")
  
  # Create Summary sheet
  cat("\nCreating Summary sheet...\n")
  wb$add_worksheet("Summary")
  
  # Calculate statistics
  total_loci <- length(unique(result$Locus))
  loci_with_genes <- sum(!is.na(result$Representative_Gene) & result$Representative_Gene != "")
  # Count unique loci with genes (not just rows)
  loci_with_genes <- length(unique(result[!is.na(GENE) & GENE_NAME != "No significant genes", Locus]))
  highest_gene_count <- max(table(result$Locus))
  lowest_gene_count <- min(table(result[GENE_NAME != "No significant genes", Locus]))
  
  summary_data <- data.frame(
    Metric = c("Total Number of Loci", 
               "Number of Loci with Significant Genes",
               "Number of Loci without Significant Genes",
               "Highest Number of Genes per Locus",
               "Lowest Number of Genes per Locus",
               "Total Genes Analyzed",
               "P-value Threshold",
               "Distance Buffer (kb)"),
    Value = c(total_loci,
              loci_with_genes,
              total_loci - loci_with_genes,
              highest_gene_count,
              ifelse(is.finite(lowest_gene_count), lowest_gene_count, 0),
              nrow(result[GENE_NAME != "No significant genes"]),
              opt$`p-threshold`,
              opt$`distance-kb`)
  )
  
  wb$add_data(sheet="Summary", x=summary_data, start_col=1, start_row=1)
  
  # Style summary sheet
  wb$add_fill(sheet="Summary", dims="A1:B1", color=wb_color(hex="4472C4"))
  wb$add_font(sheet="Summary", dims="A1:B1", color=wb_color(hex="FFFFFF"), bold=TRUE)
  wb$add_cell_style(sheet="Summary", dims="A1:B1", horizontal="center", vertical="center")
  wb$add_cell_style(sheet="Summary", dims="A2:B9", horizontal="center", vertical="center")
  wb$set_col_widths(sheet="Summary", cols=1:2, widths="auto")
  
  # Create Loci_Representative sheet
  cat("Creating Loci_Representative sheet...\n")
  wb$add_worksheet("Loci_Representative")
  
  # Extract unique loci with representative gene info - preserve column order from main sheet
  loci_rep_cols <- c(locus_info_cols, additional_locus_cols, representative_cols)
  loci_rep <- unique(result[, loci_rep_cols, with=FALSE])
  
  # Replace empty Representative_Gene with "No significant genes"
  loci_rep[Representative_Gene == "" | is.na(Representative_Gene), Representative_Gene := "No significant genes"]
  
  # Sort by chromosome and start position
  loci_rep <- loci_rep[order(Locus_Chr, Locus_Start)]
  
  wb$add_data(sheet="Loci_Representative", x=loci_rep, start_col=1, start_row=1)
  
  # Format Representative_Gene_P as scientific notation
  rep_p_idx <- which(names(loci_rep) == "Representative_Gene_P")
  if (length(rep_p_idx) > 0) {
    rep_p_letter <- LETTERS[rep_p_idx]
    wb$add_numfmt(sheet="Loci_Representative", 
                  dims=paste0(rep_p_letter, "2:", rep_p_letter, nrow(loci_rep)+1), 
                  numfmt="0.00E+00")
  }
  
  # Style Loci_Representative sheet header
  loci_rep_header_range <- paste0("A1:", LETTERS[ncol(loci_rep)], "1")
  wb$add_fill(sheet="Loci_Representative", dims=loci_rep_header_range, color=wb_color(hex="4472C4"))
  wb$add_font(sheet="Loci_Representative", dims=loci_rep_header_range, color=wb_color(hex="FFFFFF"), bold=TRUE)
  wb$add_cell_style(sheet="Loci_Representative", dims=loci_rep_header_range, 
                    horizontal="center", vertical="center", textRotation=90)
  
  # Style data cells
  loci_rep_data_range <- paste0("A2:", LETTERS[ncol(loci_rep)], nrow(loci_rep)+1)
  wb$add_cell_style(sheet="Loci_Representative", dims=loci_rep_data_range, 
                    horizontal="center", vertical="center")
  
  # Apply color formatting to Representative_Gene column
  rep_gene_idx_loci <- which(names(loci_rep) == "Representative_Gene")
  if (length(rep_gene_idx_loci) > 0) {
    rep_gene_letter_loci <- LETTERS[rep_gene_idx_loci]
    
    for (i in 1:nrow(loci_rep)) {
      row_num <- i + 1
      p_val <- loci_rep$Representative_Gene_P[i]
      
      if (!is.na(p_val)) {
        if (p_val < 5E-8) {
          color_hex <- "FFC0CB"  # Pink
        } else if (p_val < 5E-5) {
          color_hex <- "90EE90"  # Green
        } else if (p_val < 5E-3) {
          color_hex <- "FFFF99"  # Yellow
        } else {
          color_hex <- "FFFFFF"  # White
        }
        
        cell_range <- paste0(rep_gene_letter_loci, row_num)
        wb$add_fill(sheet="Loci_Representative", dims=cell_range, color=wb_color(hex=color_hex))
      }
    }
  }
  
  wb$set_col_widths(sheet="Loci_Representative", cols=1:ncol(loci_rep), widths="auto")
  
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
