# ORA Analysis for Significant Modules (Skin & Liver) from WGCNA
# Author: V. Mantzana-Oikonomaki
# Description: Performs over-representation analysis (ORA) using the 'genekitr' package
# on selected WGCNA modules using Gene Ontology terms.
# Requirements: genekitr, geneset, and properly formatted WGCNA module files.
# Output: CSV files of enriched terms per module saved in 'enrichment_results/'.

# ================================
# Load Required Libraries
# ================================
library(genekitr)
library(geneset)

# ================================
# Set Working Directory
# ================================
setwd("D:/transcriptomic_analysis/WGCNA")

# ================================
# Define Module Colors
# ================================
skin_colors <- c("antiquewhite2", "antiquewhite4", "coral4", "firebrick2", "saddlebrown", "skyblue2")
liver_colors <- c("bisque1", "brown1", "cyan", "darkseagreen2", "forestgreen", "green3",
                  "khaki1", "lemonchiffon4", "lightcyan1", "lightpink", "paleturquoise",
                  "plum1", "royalblue2", "royalblue3", "salmon4", "sienna1", "thistle1", 
                  "violet", "wheat1")

# ================================
# Read Module Files Function
# ================================
read_selected_modules <- function(directory, colors) {
  files <- list.files(directory, pattern = "^Module_.*\\.csv$", full.names = TRUE)
  selected_files <- files[basename(files) %in% paste0("Module_", colors, ".csv")]
  module_list <- lapply(selected_files, read.csv)
  names(module_list) <- gsub("Module_|\\.csv", "", basename(selected_files))
  return(module_list)
}

# ================================
# Load Skin & Liver Modules
# ================================
skin_modules <- read_selected_modules("./skin_modules", skin_colors)
liver_modules <- read_selected_modules("./liver_modules", liver_colors)

# ================================
# Define Safe ORA Function
# ================================
safe_genORA <- function(id, geneset) {
  tryCatch({
    genORA(id, geneset)
  }, error = function(e) {
    message("Error in genORA: ", e$message)
    return(NULL)
  })
}

# ================================
# Prepare Gene Ontology Gene Sets
# ================================
hg_go_bp <- geneset::getGO(org = "human", ont = "bp")
hg_go_mf <- geneset::getGO(org = "human", ont = "mf")
hg_go_cc <- geneset::getGO(org = "human", ont = "cc")

# ================================
# Extract Gene IDs from Module
# ================================
get_module_genes <- function(module_data) {
  if (!"x" %in% colnames(module_data)) {
    stop("Column 'x' (gene IDs) not found in module file.")
  }
  return(module_data$x)
}

# ================================
# Run ORA for a Single Module
# ================================
run_enrichment <- function(module_data) {
  ids <- get_module_genes(module_data)
  
  if (length(ids) == 0) {
    message("Warning: No gene IDs found in module.")
    return(NULL)
  }
  
  list(
    GO_BP = safe_genORA(ids, hg_go_bp),
    GO_MF = safe_genORA(ids, hg_go_mf),
    GO_CC = safe_genORA(ids, hg_go_cc),
    ORA   = safe_genORA(ids, hg_go_bp)  # You can customize if ORA uses different genesets
  )
}

# ================================
# Run ORA for All Modules
# ================================
results <- list()

# Skin Modules
for (color in names(skin_modules)) {
  cat("Running ORA for skin module:", color, "\n")
  results[[paste0("skin_", color)]] <- run_enrichment(skin_modules[[color]])
}

# Liver Modules
for (color in names(liver_modules)) {
  cat("Running ORA for liver module:", color, "\n")
  results[[paste0("liver_", color)]] <- run_enrichment(liver_modules[[color]])
}

# ================================
# Save Results to CSV
# ================================
if (!dir.exists("enrichment_results")) dir.create("enrichment_results")

save_results_to_csv <- function(results_list) {
  for (module_name in names(results_list)) {
    module_result <- results_list[[module_name]]
    if (!is.null(module_result)) {
      for (category in names(module_result)) {
        result_data <- module_result[[category]]
        if (!is.null(result_data)) {
          out_path <- paste0("enrichment_results/", module_name, "_", category, ".csv")
          write.csv(result_data, out_path, row.names = FALSE)
          cat("Saved:", out_path, "\n")
        }
      }
    }
  }
}

save_results_to_csv(results)
