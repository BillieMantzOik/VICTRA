# WGCNA analysis for liver samples of O. vicentei
# Script by V. Mantzana-Oikonomaki
# Based on WGCNA tutorial: https://www.dropbox.com/scl/fo/4vqfiysan6rlurfo2pbnk/h?rlkey=thqg8wlpdn4spu3ihjuc1kmlu&e=3&dl=0

# ===============================
# Load libraries
# ===============================
library(WGCNA)
library(dplyr)
library(tidyverse)
library(magrittr)
library(ggplot2)
library(knitr)
library(limma)
library(reshape2)
library(RColorBrewer)

# ===============================
# Set working directory
# ===============================
setwd("D:/transcriptomic_analysis")

# ===============================
# Load expression data
# ===============================
liver <- read.csv("./limma_analysis_final/results/vicentei_liver_logTMM_counts.csv", row.names = 1)
liver <- t(liver)
expr_liver <- liver

# ===============================
# Load and filter metadata
# ===============================
liver_data <- read.csv("./metadata/vicentei_frogs_metadata_liver.csv", row.names = 1)
liver_data <- liver_data %>% filter(!row_number() %in% 10)

# ===============================
# WGCNA Setup
# ===============================
allowWGCNAThreads()

powers <- c(1:40)
sft <- pickSoftThreshold(expr_liver, powerVector = powers, verbose = 5)

# Plot scale independence
sft_df <- data.frame(sft$fitIndices) %>%
  mutate(model_fit = -sign(slope) * SFT.R.sq)

ggplot(sft_df, aes(x = Power, y = model_fit, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.80, col = "red") +
  ylim(c(min(sft_df$model_fit), 1.05)) +
  xlab("Soft Threshold (power)") +
  ylab("Scale Free Topology Model Fit, signed R^2") +
  ggtitle("Scale independence") +
  theme_classic()

# ===============================
# Construct network
# ===============================
power <- sft$powerEstimate
adjacency <- adjacency(expr_liver, power = 11)
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM

# ===============================
# Module detection
# ===============================
geneTree <- hclust(as.dist(dissTOM), method = "average")
plot(geneTree, main = "Gene clustering on TOM-based dissimilarity_liver", labels = FALSE, hang = 0.04)

dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2,
                             pamRespectsDendro = FALSE, minClusterSize = 10)

dynamicColors <- labels2colors(dynamicMods)
pdf("4-module_tree.pdf", width = 8, height = 6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut Liver",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

# ===============================
# Eigengenes and hub genes
# ===============================
MEList <- moduleEigengenes(expr_liver, colors = dynamicColors)
MEs <- MEList$eigengenes
write.csv(MEs, "./WGCNA/liver_modules/ModEigengenes_liver.csv")

KME <- signedKME(expr_liver, MEs, corFnc = 'cor', corOptions = "use = 'p'")
write.csv(KME, "./WGCNA/liver_modules/KME_liver.csv")

hub_genes <- chooseTopHubInEachModule(expr_liver, dynamicColors,
                                      omitColors = "grey", power = 2, type = "signed")
write.csv(hub_genes, "./WGCNA/liver_modules/hub_genes_per_module_liver.csv")

# ===============================
# Module-trait correlation
# ===============================
traits <- liver_data %>% select(S1R, H4, CAA.Jc)
moduleTraitCor <- cor(MEs, traits, use = "p")
write.csv(moduleTraitCor, "./WGCNA/liver_modules/MTC_liver.csv")

nSamples <- nrow(expr_liver)
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
write.csv(moduleTraitPvalue, "./WGCNA/liver_modules/MTPv_liver.csv")

# Filter for significant module-trait correlations
abs_MTC <- abs(moduleTraitCor)
subset_MTC <- moduleTraitCor[apply(abs_MTC, 1, function(row) any(row > 0.8)), ]
write.csv(subset_MTC, "./WGCNA/liver_modules/MTC_signif_liver.csv")

# Heatmap
p_values_matrix <- as.matrix(moduleTraitPvalue)
correlations_matrix <- as.matrix(subset_MTC)

textMatrix <- matrix(NA, nrow = nrow(correlations_matrix), ncol = ncol(correlations_matrix))
for (i in 1:nrow(correlations_matrix)) {
  for (j in 1:ncol(correlations_matrix)) {
    textMatrix[i, j] <- paste0("cor = ", format(correlations_matrix[i, j], scientific = FALSE, digits = 2),
                               "\np = ", format(p_values_matrix[i, j], digits = 3))
  }
}

tiff("./liver_trait_cor.tif", width = 8, height = 8, units = "in", res = 600)
labeledHeatmap(Matrix = correlations_matrix,
               xLabels = colnames(correlations_matrix),
               yLabels = rownames(correlations_matrix),
               ySymbols = rownames(correlations_matrix),
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               main = "Module-trait relationships")
dev.off()

# ===============================
# Save genes per module
# ===============================
geneModuleDF <- data.frame(Gene = colnames(expr_liver), Module = dynamicColors)
moduleGenesList <- split(geneModuleDF$Gene, geneModuleDF$Module)
for (module in names(moduleGenesList)) {
  write.csv(moduleGenesList[[module]], paste0("Module_", module, ".csv"), row.names = FALSE)
}

# ===============================
# Cytoscape export
# ===============================
modsi <- read.csv("./WGCNA/liver_modules/MTC_signif_liver.csv")
hubGenes <- read.csv("./WGCNA/liver_modules/hub_genes_per_module_liver.csv")
hubGenesSubset <- hubGenes[hubGenes$X %in% modsi$X, ]
write.csv(hubGenesSubset, "./WGCNA/liver_modules/hub_genes_liver_signif.csv")

mods <- unique(sub("^ME", "", modsi$X))

for (module in mods) {
  inModule <- dynamicColors == module
  genesInModule <- colnames(expr_liver)[inModule]
  modTOM <- TOM[inModule, inModule]
  dimnames(modTOM) <- list(genesInModule, genesInModule)

  hubGene <- hubGenesSubset$x[hubGenesSubset$X == module]
  isHub <- genesInModule %in% hubGene
  nodeAttr <- setNames(as.character(isHub), genesInModule)

  exportNetworkToCytoscape(modTOM,
                           edgeFile = paste0("./WGCNA/liver_modules_net/CytoscapeInput-edges-", module, ".txt"),
                           nodeFile = paste0("./WGCNA/liver_modules_net/CytoscapeInput-nodes-", module, ".txt"),
                           weighted = TRUE, threshold = 0.02,
                           nodeNames = genesInModule, altNodeNames = genesInModule,
                           nodeAttr = nodeAttr)
}

# Add IsHub column
for (module in mods) {
  nodeFile <- paste0("./WGCNA/liver_modules_net/CytoscapeInput-nodes-", module, ".txt")
  nodes <- read.table(nodeFile, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  hubGenesInModule <- hubGenesSubset$x[hubGenesSubset$X == module]
  nodes$IsHub <- nodes$nodeName %in% hubGenesInModule
  write.table(nodes, file = nodeFile, sep = "\t", quote = FALSE, row.names = FALSE)
}

# ===============================
# DEG and color gene enrichment
# ===============================
deg_genes <- read.csv("./limma_analysis_final/results/limma_all_contr_liver_significant.csv", row.names = 1) %>% rownames()
col_genes <- read.csv("./genes list literature/literature_color_genes.csv")$gene

MTC_signif <- read.csv("./WGCNA/liver_modules/MTC_signif_liver.csv", row.names = 1)
significant_modules <- rownames(MTC_signif)
module_dir <- "./WGCNA/liver_modules"
module_files <- list.files(module_dir, pattern = "*.csv", full.names = TRUE)

results <- data.frame(Module = character(),
                      DEGs_Count = integer(),
                      DEGs_Percentage = numeric(),
                      Color_Genes_Count = integer(),
                      Color_Genes_Percentage = numeric(),
                      stringsAsFactors = FALSE)

for (module_name in significant_modules) {
  file_pattern <- paste0("Module_", gsub("^ME", "", module_name), ".csv")
  module_file <- module_files[grepl(file_pattern, basename(module_files))]

  if (length(module_file) > 0) {
    module_genes <- read.csv(module_file)$x
    degs_in_module <- sum(module_genes %in% deg_genes)
    deg_percentage <- (degs_in_module / length(module_genes)) * 100
    color_genes_in_module <- sum(module_genes %in% col_genes)
    color_genes_percentage <- (color_genes_in_module / length(module_genes)) * 100

    results <- rbind(results, data.frame(
      Module = module_name,
      DEGs_Count = degs_in_module,
      DEGs_Percentage = deg_percentage,
      Color_Genes_Count = color_genes_in_module,
      Color_Genes_Percentage = color_genes_percentage
    ))
  } else {
    warning(paste("Module file not found for:", module_name))
  }
}

write.csv(results, "./WGCNA/liver_modules/module_gene_analysis_results.csv", row.names = FALSE)