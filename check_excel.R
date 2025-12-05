#!/usr/bin/env Rscript
# Quick script to check Excel file contents

library(openxlsx2)
library(data.table)

excel_file <- "/s/nath-lab/alsamman/___Analysis___/Hispanic_2025/MAGMA_PipiLIne/results/example_analysis/magma/geneset_analysis/example_analysis_geneset_report.xlsx"

wb <- wb_load(excel_file)
sheet_names <- wb$get_sheet_names()
cat("Sheets in workbook:", paste(sheet_names, collapse=", "), "\n\n")

if ("Geneset_Loci_Matrix" %in% sheet_names) {
  data <- wb_read(wb, sheet="Geneset_Loci_Matrix")
  cat("Columns in Geneset_Loci_Matrix:\n")
  print(names(data))
  cat("\n")
  
  if ("N_Genes_Pass_Threshold" %in% names(data)) {
    cat("N_Genes_Pass_Threshold column found!\n")
    cat("Sample values (first 10 rows):\n")
    print(data[1:10, c("Geneset", "Total_Genes", "N_Genes_Pass_Threshold")])
    cat("\nSummary:\n")
    cat("Min:", min(data$N_Genes_Pass_Threshold, na.rm=TRUE), "\n")
    cat("Max:", max(data$N_Genes_Pass_Threshold, na.rm=TRUE), "\n")
    cat("Gene-sets with N > 0:", sum(data$N_Genes_Pass_Threshold > 0, na.rm=TRUE), "\n")
  } else {
    cat("N_Genes_Pass_Threshold column NOT found!\n")
    cat("Available columns:", paste(names(data), collapse=", "), "\n")
  }
}
