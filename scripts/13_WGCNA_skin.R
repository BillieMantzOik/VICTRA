# WGCNA Analysis for Skin and Liver Samples of O. vicentei
# Script by V. Mantzana-Oikonomaki
# Following WGCNA tutorial: https://www.dropbox.com/scl/fo/4vqfiysan6rlurfo2pbnk/h?rlkey=thqg8wlpdn4spu3ihjuc1kmlu&e=3&dl=0

# Load necessary libraries
library(WGCNA)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(knitr)
library(limma)
library(reshape2)
library(RColorBrewer)

# Set working directory
setwd("D:/transcriptomic_analysis")

# Import data
# Skin count data (already TMM normalized)
skin <- read.csv("vicentei_skin_logTMM_counts.csv", row.names = 1)
skin <- t(skin)  # Transpose data: genes as columns, samples as rows
expr_skin <- skin

# Import metadata for skin
skin_data <- read.csv("./metadata/vicentei_frogs_metadata_skin.csv", row.names = 1)
skin_data <- skin_data %>% filter(!row_number() %in% 7)  # Remove problematic row (cei5s)

# Select relevant traits (based on PCA results)
traits <- skin_data %>% dplyr::select(S1R, H4, CAA.Jc, melanophore)

# Enable multithreading for WGCNA (12 threads)
allowWGCNAThreads() 

# Choose a set of soft threshold parameters for network construction
powers = c(1:20, seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(expr_skin, powerVector = powers, verbose = 2)

# Plot to select the optimal soft threshold
sft_df <- data.frame(sft$fitIndices) %>% 
  mutate(model_fit = -sign(slope) * SFT.R.sq)
ggplot(sft_df, aes(x = Power, y = model_fit, label = Power)) +
  geom_point() + geom_text(nudge_y = 0.1) + 
  geom_hline(yintercept = 0.80, col = "red") + 
  ylim(c(min(sft_df$model_fit), 1.05)) +
  xlab("Soft Threshold (power)") + 
  ylab("Scale Free Topology Model Fit, signed R^2") +
  ggtitle("Scale independence") + theme_classic()

# Set the soft threshold power (chosen as 7 based on the plot)
power = 7

# Construct adjacency matrix and topological overlap matrix (TOM)
adjacency = adjacency(expr_skin, power = power, type = 'signed')
TOM = TOMsimilarity(adjacency, TOMType = 'signed')
dissTOM = 1 - TOM  # Dissimilarity from TOM

# Perform hierarchical clustering on the TOM dissimilarity
geneTree = hclust(as.dist(dissTOM), method = "average")

# Plot the gene dendrogram
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", 
     labels = FALSE, hang = 0.04)

# Identify modules using dynamic tree cut
minModuleSize = 10
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, 
                            pamRespectsDendro = FALSE, minClusterSize = minModuleSize)

# Convert numeric labels to colors
dynamicColors = labels2colors(dynamicMods)

# Plot the dendrogram with module colors
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, 
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05, 
                    main = "Gene dendrogram and module colors")

# Calculate module eigengenes
MEList = moduleEigengenes(expr_skin, colors = dynamicColors)
MEs = MEList$eigengenes
write.csv(MEs, "ModEigengenes_expr_skin.csv")

# Calculate module membership (KME)
KME = signedKME(expr_skin, MEs, corFnc = 'cor', corOptions = "use = 'p'")
write.csv(KME, "skin_KME.csv")

# Identify hub genes in each module
hub_genes <- chooseTopHubInEachModule(expr_skin, dynamicColors, omitColors = "grey", 
                                      power = power, type = "signed")
write.csv(hub_genes, "hub_genes_per_module_skin.csv")

# Correlate modules with traits
MTC = cor(MEs, traits, method = "pearson", use = "p")
write.csv(MTC, "./WGCNA/skin_modules/skin_module_trait_correlation.csv")

# Subset to strongly correlated modules (>0.75)
abs_MTC <- abs(MTC)
rows_to_keep <- apply(abs_MTC, 1, function(row) any(row > 0.75))
subset_MTC <- MTC[rows_to_keep, ]
write.csv(subset_MTC, "MTC_significant_skin.csv")

# Compute p-values for correlation significance
MTP = corPvalueStudent(subset_MTC, nSamples = nrow(expr_skin))
write.csv(MTP, "pvalue_skin.csv")

# Create a labeled heatmap of module-trait correlations
textMatrix <- matrix(NA, nrow = nrow(subset_MTC), ncol = ncol(subset_MTC))
for (i in 1:nrow(subset_MTC)) {
  for (j in 1:ncol(subset_MTC)) {
    textMatrix[i, j] <- paste0("cor = ", format(subset_MTC[i, j], scientific = FALSE, digits = 2), 
                               "\np = ", format(MTP[i, j], digits = 3))
  }
}

tiff("./skin_trait_cor.tif", width = 8, height = 8, units = "in", res = 600)
labeledHeatmap(Matrix = subset_MTC, xLabels = colnames(subset_MTC), yLabels = rownames(subset_MTC),
               ySymbols = rownames(subset_MTC), colorLabels = FALSE, colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 1, textAdj = c(0.5, 0.5), 
               main = "Module-trait relationships")
dev.off()

# Export module genes for Cytoscape
geneModuleDF <- data.frame(Gene = colnames(expr_skin), Module = dynamicColors)
moduleGenesList <- split(geneModuleDF$Gene, geneModuleDF$Module)

for(module in names(moduleGenesList)) {
  fileName <- paste0("./WGCNA/skin_modules_new/Module_", module, ".csv")
  write.csv(moduleGenesList[[module]], file = fileName, row.names = FALSE)
}

# Extract hub genes for significant modules and prepare for Cytoscape
modsi <- read.csv("MTC_significant_skin.csv")
hubGenes <- read.csv("hub_genes_per_module_skin.csv")
hubGenesSubset <- hubGenes[hubGenes$X %in% modsi, ]
write.csv(hubGenesSubset, "./WGCNA/skin_modules_new/hub_genes_skin_significant.csv")

# Export networks for Cytoscape (edges and nodes)
for (module in unique(sub("^ME", "", modsi$X))) {
  inModule <- (dynamicColors == module)
  genesInModule <- colnames(expr_skin)[inModule]
  modTOM <- TOM[inModule, inModule]
  dimnames(modTOM) <- list(genesInModule, genesInModule)
  
  hubGene <- hubGenesSubset$x[hubGenesSubset$X == module]
  isHub <- genesInModule %in% hubGene
  nodeAttr <- setNames(as.character(isHub), genesInModule)
  
  edgeFile <- paste0("./WGCNA/skin_modules_net/CytoscapeInput-edges-", module, ".txt")
  nodeFile <- paste0("./WGCNA/skin_modules_net/CytoscapeInput-nodes-", module, ".txt")
  
  exportNetworkToCytoscape(modTOM, edgeFile = edgeFile, nodeFile = nodeFile, weighted = TRUE, 
                           threshold = 0.02, nodeNames = genesInModule, altNodeNames = genesInModule, 
                           nodeAttr = nodeAttr)
}

# Final DEG and color-related gene analysis per module
module_dir <- "./WGCNA/skin_modules"
deg_genes <- read.csv("./limma_analysis_final/results/limma_all_contr_skin_significant.csv", row.names = 1)
deg_genes <- row.names(deg_genes)
col_genes <- read.csv("./genes list literature/literature_color_genes.csv")$gene
MTC_signif <- read.csv("./WGCNA/skin_modules/MTC_significant_skin.csv", row.names = 1)
significant_modules <- row.names(MTC_signif)

results <- data.frame(Module = character(), DEGs_Count = integer(), DEGs_Percentage = numeric(),
                      Color_Genes_Count = integer(), Color_Genes_Percentage = numeric(), stringsAsFactors = FALSE)

module_files <- list.files(module_dir, pattern = "*.csv", full.names = TRUE)

for (module_name in significant_modules) {
  file_pattern <- paste0("Module_", gsub("^ME", "", module_name), ".csv")
  module_file <- module_files[grepl(file_pattern, basename(module_files))]
  
  if (length(module_file) > 0) {
    module_genes <- read.csv(module_file)$x
    degs_in_module <- sum(module_genes %in% deg_genes)
    deg_percentage <- (degs_in_module / length(module_genes)) * 100
    color_genes_in_module <- sum(module_genes %in% col_genes)
    color_genes_percentage <- (color_genes_in_module / length(module_genes)) * 100
    
    results <- rbind(results, data.frame(Module = module_name, DEGs_Count = degs_in_module,
                                         DEGs_Percentage = deg_percentage, Color_Genes_Count = color_genes_in_module,
                                         Color_Genes_Percentage = color_genes_percentage))
  } else {
    warning(paste("Module file not found for:", module_name))
  }
}

# View and save the results
print(results)
write.csv(results, "./WGCNA/skin_modules/module_gene_analysis_results.csv", row.names = FALSE)
