#!/usr/bin/env Rscript

# Plot SNP-based GWAS vs Gene-based MAGMA results
# Generates Manhattan and Miami plots for comparison

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
  
  # Set personal library path
  personal_lib <- Sys.getenv("R_LIBS_USER")
  if (personal_lib == "") {
    personal_lib <- path.expand("~/R/library")
  }
  
  # Create personal library if it doesn't exist
  if (!dir.exists(personal_lib)) {
    dir.create(personal_lib, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Add to library paths
  .libPaths(c(personal_lib, .libPaths()))
  
  # Install manhattantwin if not available
  if (!require("manhattantwin", quietly = TRUE)) {
    cat("Installing manhattantwin from GitHub to personal library:", personal_lib, "\n")
    if (!require("devtools", quietly = TRUE)) {
      install.packages("devtools", lib = personal_lib, repos = "https://cloud.r-project.org")
      library(devtools, lib.loc = personal_lib)
    }
    devtools::install_github("AlsammanAlsamman/manhattantwin", lib = personal_lib)
    library(manhattantwin, lib.loc = personal_lib)
  } else {
    library(manhattantwin)
  }
})

# Parse command line arguments
option_list <- list(
  make_option("--gwas", type="character", help="GWAS summary statistics file"),
  make_option("--genes", type="character", help="MAGMA gene results file (.genes.raw)"),
  make_option("--gene_loc", type="character", help="Gene location file (NCBI format)"),
  make_option("--loci", type="character", help="Loci information file"),
  make_option("--analysis", type="character", help="Analysis name"),
  make_option("--snp_p_threshold", type="double", default=5e-2, help="SNP P-value threshold"),
  make_option("--gene_p_threshold", type="double", default=5e-8, help="Gene P-value threshold"),
  make_option("--cluster_distance", type="integer", default=250000, help="SNP clustering distance (bp)"),
  make_option("--n_cases", type="character", help="Number of cases"),
  make_option("--n_controls", type="character", help="Number of controls"),
  make_option("--plot_per_chromosome", type="logical", default=TRUE, help="Generate per-chromosome plots"),
  make_option("--plot_only_sig_chromosomes", type="logical", default=TRUE, help="Only plot chromosomes with data"),
  make_option("--genome_wide_dir", type="character", help="Output directory for genome-wide plots"),
  make_option("--per_chromosome_dir", type="character", help="Output directory for per-chromosome plots"),
  make_option("--filtered_data_dir", type="character", help="Output directory for filtered data"),
  make_option("--output_done", type="character", help="Done marker file")
)

opt <- parse_args(OptionParser(option_list=option_list))

cat("========================================\n")
cat("SNP-Gene Comparison Plot Generation\n")
cat("========================================\n")
cat("Analysis:", opt$analysis, "\n")
cat("SNP P-value threshold:", opt$snp_p_threshold, "\n")
cat("Gene P-value threshold:", opt$gene_p_threshold, "\n")
cat("Clustering distance:", opt$cluster_distance, "bp\n")
cat("N cases:", opt$n_cases, "\n")
cat("N controls:", opt$n_controls, "\n\n")

# Create output directories
dir.create(opt$genome_wide_dir, recursive=TRUE, showWarnings=FALSE)
if (opt$plot_per_chromosome) {
  dir.create(opt$per_chromosome_dir, recursive=TRUE, showWarnings=FALSE)
}
dir.create(opt$filtered_data_dir, recursive=TRUE, showWarnings=FALSE)

# ============================================================================
# 1. Load and process GWAS data
# ============================================================================
cat("Loading GWAS summary statistics...\n")
gwas <- fread(opt$gwas, header=TRUE)
cat("  Original GWAS rows:", nrow(gwas), "\n")

# Standardize column names - convert to lowercase first
setnames(gwas, 
         old = names(gwas), 
         new = tolower(names(gwas)), 
         skip_absent = TRUE)

# Map column names to expected format
col_mapping <- c(
  chr = "chrom",
  bp = "pos", 
  snp = "rsid",
  p = "pvalue"
)

# Check which columns exist and rename them
for (new_name in names(col_mapping)) {
  old_name <- col_mapping[new_name]
  if (old_name %in% names(gwas)) {
    setnames(gwas, old_name, new_name)
  }
}

# Ensure required columns exist after mapping
required_cols <- c("chr", "bp", "snp", "p")
missing_cols <- setdiff(required_cols, names(gwas))
if (length(missing_cols) > 0) {
  stop("GWAS file missing required columns after mapping: ", paste(missing_cols, collapse=", "))
}

# Keep only autosomes (1-22)
gwas <- gwas[chr %in% 1:22]
cat("  After filtering to autosomes:", nrow(gwas), "\n")

# Remove NA P-values
gwas <- gwas[!is.na(p)]
cat("  After removing NA P-values:", nrow(gwas), "\n")

# ============================================================================
# 2. Load loci information and assign locus labels to SNPs
# ============================================================================
cat("\nLoading loci information...\n")
loci <- fread(opt$loci, header=TRUE)
cat("  Loaded", nrow(loci), "loci\n")

# Assign locus to each SNP based on proximity
cat("  Assigning locus labels to SNPs...\n")
gwas[, locus := NA_character_]

for (i in 1:nrow(loci)) {
  locus_chr <- loci$Chr[i]
  locus_start <- loci$Start[i]
  locus_end <- loci$End[i]
  locus_name <- loci$Locus[i]
  
  # Assign locus name to SNPs within this region
  gwas[chr == locus_chr & bp >= locus_start & bp <= locus_end, locus := locus_name]
}

n_snps_with_locus <- gwas[!is.na(locus), .N]
cat("  ", n_snps_with_locus, "SNPs assigned to loci\n")

# Use locus name as the label for SNPs (for grouping in manhattan plot)
gwas[, snp_label := locus]

# ============================================================================
# 4. Load and process gene data
# ============================================================================
cat("\nLoading MAGMA gene results...\n")
# Read .genes.out file with columns: GENE CHR START STOP NSNPS NPARAM N ZSTAT P
genes <- fread(opt$genes, header=TRUE)
cat("  Original gene rows:", nrow(genes), "\n")
cat("  Gene file columns:", paste(names(genes), collapse=", "), "\n")

# Standardize column names to uppercase if needed
setnames(genes, toupper(names(genes)))

# Load gene location file for gene names
cat("Loading gene location file...\n")
gene_loc <- fread(opt$gene_loc, header=FALSE)
setnames(gene_loc, c("GENE", "CHR_LOC", "START_LOC", "END_LOC", "STRAND", "GENE_NAME"))
cat("  Loaded", nrow(gene_loc), "gene locations\n")

# Merge to get gene names (genes.out already has CHR, START, STOP)
genes_with_names <- merge(genes, gene_loc[, .(GENE, GENE_NAME)], by="GENE", all.x=TRUE)
cat("  After merging with gene names:", nrow(genes_with_names), "\n")

# Calculate gene middle position from START and STOP in genes.out
# Note: genes.out uses STOP instead of END
genes_with_names[, gene_middle := as.integer((START + STOP) / 2)]

# Keep only autosomes
genes_with_names <- genes_with_names[CHR %in% 1:22]
cat("  After filtering to autosomes:", nrow(genes_with_names), "\n")

# Remove genes with missing positions or P-values
genes_with_names <- genes_with_names[!is.na(gene_middle) & !is.na(P)]
cat("  After removing NAs:", nrow(genes_with_names), "\n")

# Use GENE_NAME if available, otherwise use GENE
genes_with_names[, gene_label := ifelse(!is.na(GENE_NAME), GENE_NAME, as.character(GENE))]

# Prepare gene data for Manhattan plot format
gene_data <- genes_with_names[, .(
  chr = CHR,
  bp = gene_middle,
  snp = gene_label,
  p = P,
  gene = gene_label,
  label = gene_label
)]

cat("  Final gene data rows:", nrow(gene_data), "\n")

# ============================================================================
# 5. Filter data for plotting
# ============================================================================
cat("\nFiltering data for plotting...\n")

# Filter SNPs by P-value threshold
snp_filtered <- gwas[p <= opt$snp_p_threshold]
cat("  SNPs passing threshold (P <=", opt$snp_p_threshold, "):", nrow(snp_filtered), "\n")

# Filter genes by P-value threshold
gene_filtered <- gene_data[p <= opt$gene_p_threshold]
cat("  Genes passing threshold (P <=", opt$gene_p_threshold, "):", nrow(gene_filtered), "\n")

# Save filtered data with additional columns
snp_out <- snp_filtered[, .(chr, pos=bp, snp, pvalue=p, locus, label=snp_label)]
gene_out <- gene_filtered[, .(chr, pos=bp, gene=snp, pvalue=p, label)]

fwrite(snp_out, paste0(opt$filtered_data_dir, "/", opt$analysis, "_snps_filtered.txt"), 
       sep="\t", quote=FALSE)
fwrite(gene_out, paste0(opt$filtered_data_dir, "/", opt$analysis, "_genes_filtered.txt"), 
       sep="\t", quote=FALSE)

cat("  Filtered data saved to:", opt$filtered_data_dir, "\n")

# ============================================================================
# 6. Prepare data for manhattantwin (genome-wide)
# ============================================================================
cat("\nPreparing genome-wide plot data...\n")

# For manhattantwin, prepare data with required columns
# SNP data: chr, pos, pvalue, gene (locus name for grouping)
snp_data <- data.frame(
  chr = snp_filtered$chr,
  pos = snp_filtered$bp,
  pvalue = snp_filtered$p,
  gene = snp_filtered$locus,  # Use locus name for grouping
  rsid = snp_filtered$snp
)

# Gene data: chr, pos, pvalue, gene (gene name)
gene_data_plot <- data.frame(
  chr = gene_filtered$chr,
  pos = gene_filtered$bp,
  pvalue = gene_filtered$p,
  gene = gene_filtered$gene,  # Use gene name
  rsid = gene_filtered$gene
)

cat("  SNP plot data:", nrow(snp_data), "rows\n")
cat("  Gene plot data:", nrow(gene_data_plot), "rows\n")

# ============================================================================
# 7. Generate genome-wide plots
# ============================================================================
cat("\nGenerating genome-wide Manhattan plot (SNPs)...\n")

# Single Manhattan plot for SNPs using plot_single_manhattan
manhattantwin::plot_single_manhattan(
  snp_data,
  plot_title_prefix = paste0("SNP-based GWAS - ", opt$analysis),
  chr_col = "chr",
  bp_col = "pos",
  p_col = "pvalue",
  n_cases = as.numeric(opt$n_cases),
  n_controls = as.numeric(opt$n_controls),
  file_name_prefix = paste0(opt$analysis, "_snp"),
  group_col = "gene",  # Group by locus
  gene_col = "gene",   # Label by locus
  output_folder = opt$genome_wide_dir,
  y_axis_squish_threshold = 30,
  label_threshold_colors = c("red" = 5e-8, "orange" = opt$snp_p_threshold),
  output_width = 15,
  output_height = 8
)

cat("  Saved: SNP Manhattan plot\n")

# Miami plot (SNPs vs Genes)
cat("\nGenerating genome-wide Miami plot (SNPs vs Genes)...\n")

manhattantwin::create_mirrored_manhattan_plot(
  snp_data,
  gene_data_plot,
  chr_col = "chr",
  bp_col = "pos",
  p_col = "pvalue",
  n_cases1 = as.numeric(opt$n_cases),
  n_controls1 = as.numeric(opt$n_controls),
  n_cases2 = as.numeric(opt$n_cases),
  n_controls2 = as.numeric(opt$n_controls),
  file_name_prefix = paste0(opt$analysis, "_miami"),
  group_col = "gene",
  gene_col = "gene",
  label_threshold_colors = c("red" = 5e-8, "orange" = opt$snp_p_threshold, "darkblue" = opt$gene_p_threshold),
  output_folder = opt$genome_wide_dir,
  y_axis_squish_threshold = 30,
  plot_title1 = paste0("SNP-based GWAS - ", opt$analysis),
  plot_title2 = paste0("Gene-based MAGMA - ", opt$analysis),
  label_alpha = 0.7,
  label_orientation = "vertical",
  output_width = 15,
  output_height = 10
)

cat("  Saved: Miami plot\n")

# ============================================================================
# 8. Generate per-chromosome plots
# ============================================================================
if (opt$plot_per_chromosome) {
  cat("\nGenerating per-chromosome plots...\n")
  
  # Determine which chromosomes to plot
  if (opt$plot_only_sig_chromosomes) {
    # Only chromosomes with data after filtering
    chr_with_snps <- unique(snp_filtered$chr)
    chr_with_genes <- unique(gene_filtered$chr)
    chromosomes <- sort(unique(c(chr_with_snps, chr_with_genes)))
  } else {
    # All autosomes
    chromosomes <- 1:22
  }
  
  cat("  Plotting", length(chromosomes), "chromosomes:", paste(chromosomes, collapse=", "), "\n")
  
  for (chr_i in chromosomes) {
    cat("  Processing chromosome", chr_i, "...\n")
    
    # Filter data for this chromosome
    snp_chr <- snp_filtered[chr == chr_i]
    gene_chr <- gene_filtered[chr == chr_i]
    
    # Skip if no data
    if (nrow(snp_chr) == 0 && nrow(gene_chr) == 0) {
      cat("    No data for chromosome", chr_i, ", skipping\n")
      next
    }
    
    # Prepare data frames for manhattantwin
    snp_chr_data <- data.frame(
      chr = snp_chr$chr,
      pos = snp_chr$bp,
      pvalue = snp_chr$p,
      gene = snp_chr$locus,
      rsid = snp_chr$snp
    )
    
    gene_chr_data <- data.frame(
      chr = gene_chr$chr,
      pos = gene_chr$bp,
      pvalue = gene_chr$p,
      gene = gene_chr$gene,
      rsid = gene_chr$gene
    )
    
    # Generate Miami plot for this chromosome using manhattantwin
    if (nrow(snp_chr_data) > 0 && nrow(gene_chr_data) > 0) {
      manhattantwin::create_mirrored_manhattan_plot(
        snp_chr_data,
        gene_chr_data,
        chr_col = "chr",
        bp_col = "pos",
        p_col = "pvalue",
        n_cases1 = as.numeric(opt$n_cases),
        n_controls1 = as.numeric(opt$n_controls),
        n_cases2 = as.numeric(opt$n_cases),
        n_controls2 = as.numeric(opt$n_controls),
        file_name_prefix = paste0(opt$analysis, "_chr", chr_i, "_miami"),
        group_col = "gene",
        gene_col = "gene",
        label_threshold_colors = c("red" = 5e-8, "orange" = opt$snp_p_threshold, "darkblue" = opt$gene_p_threshold),
        output_folder = opt$per_chromosome_dir,
        y_axis_squish_threshold = 20,
        plot_title1 = paste0("Chr", chr_i, " - SNP-based GWAS"),
        plot_title2 = paste0("Chr", chr_i, " - Gene-based MAGMA"),
        label_alpha = 0.7,
        label_orientation = "vertical",
        output_width = 12,
        output_height = 8
      )
      cat("    Saved: Chr", chr_i, " Miami plot\n")
    } else if (nrow(snp_chr_data) > 0) {
      # Only SNPs
      manhattantwin::plot_single_manhattan(
        snp_chr_data,
        plot_title_prefix = paste0("Chr", chr_i, " - SNP-based GWAS"),
        chr_col = "chr",
        bp_col = "pos",
        p_col = "pvalue",
        n_cases = as.numeric(opt$n_cases),
        n_controls = as.numeric(opt$n_controls),
        file_name_prefix = paste0(opt$analysis, "_chr", chr_i, "_snp"),
        group_col = "gene",
        gene_col = "gene",
        output_folder = opt$per_chromosome_dir,
        y_axis_squish_threshold = 20,
        label_threshold_colors = c("red" = 5e-8, "orange" = opt$snp_p_threshold),
        output_width = 12,
        output_height = 6
      )
      cat("    Saved: Chr", chr_i, " SNP plot\n")
    } else {
      # Only genes
      manhattantwin::plot_single_manhattan(
        gene_chr_data,
        plot_title_prefix = paste0("Chr", chr_i, " - Gene-based MAGMA"),
        chr_col = "chr",
        bp_col = "pos",
        p_col = "pvalue",
        n_cases = as.numeric(opt$n_cases),
        n_controls = as.numeric(opt$n_controls),
        file_name_prefix = paste0(opt$analysis, "_chr", chr_i, "_gene"),
        group_col = "gene",
        gene_col = "gene",
        output_folder = opt$per_chromosome_dir,
        y_axis_squish_threshold = 20,
        label_threshold_colors = c("red" = opt$gene_p_threshold),
        output_width = 12,
        output_height = 6
      )
      cat("    Saved: Chr", chr_i, " Gene plot\n")
    }
  }
}

# ============================================================================
# 9. Create done marker
# ============================================================================
cat("\nCreating done marker...\n")
writeLines(paste0("Plotting completed: ", Sys.time()), opt$output_done)

cat("\n========================================\n")
cat("Plot generation completed successfully!\n")
cat("========================================\n")
cat("Genome-wide plots:", opt$genome_wide_dir, "\n")
if (opt$plot_per_chromosome) {
  cat("Per-chromosome plots:", opt$per_chromosome_dir, "\n")
}
cat("Filtered data:", opt$filtered_data_dir, "\n")
