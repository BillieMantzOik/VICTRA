###############################################################################
# Differential gene expression analysis using limma-voom for RNA-seq data
# Organism: O. vicentei frogs | Tissues: Skin & Liver
# V. Mantzana-Oikonomaki, 2024
###############################################################################

### Load necessary libraries (Bioconductor + CRAN)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(ggrepel)
library(tidyverse)
library(tximport)
library(readr)
library(tximportData)
library(dplyr)
library(edgeR)
library(limma)
library(DESeq2)

### Set working directory
setwd("./transcriptomic_analysis")

###############################################################################
# Skin Data Analysis
###############################################################################

# Load skin metadata
skin <- read.csv("./metadata/vicentei_frogs_metadata_skin.csv") %>%
  filter(!row_number() %in% 7) %>%
  as.data.frame()
rownames(skin) <- skin$sample

# Load transcript-to-gene mapping and Kallisto abundance files
skin_tx2gene <- read.csv("./metadata/transcript_to_gene_map.csv", stringsAsFactors = FALSE)
files_skin <- paste(skin$path, "abundance.tsv", sep = "/")
names(files_skin) <- rownames(skin)

# Import counts
txi_skin <- tximport(files_skin, type = "kallisto", tx2gene = skin_tx2gene, ignoreAfterBar = TRUE)
dge_skin <- DGEList(txi_skin$counts)

# Filter low expression genes using design
group_skin <- interaction(skin$morph)
design_skin <- model.matrix(~0 + group_skin)
keep_skin <- filterByExpr(dge_skin, design_skin)
dge_skin <- dge_skin[keep_skin, ]
dge_skin <- calcNormFactors(dge_skin, method = "TMM")

# Voom transformation
v_skin <- voom(dge_skin, design_skin, plot = TRUE)
write.csv(v_skin$E, "./limma_analysis_final/results/voom_counts_skin_new.csv")

# Linear modeling
fit_skin <- lmFit(v_skin, design_skin)
contrasts_skin <- makeContrasts(
  aquamarine_vs_other = group_skinaquamarine - ((group_skinbrown + group_skingreen + group_skinred)/3),
  brown_vs_other      = group_skinbrown - ((group_skinaquamarine + group_skingreen + group_skinred)/3),
  green_vs_other      = group_skingreen - ((group_skinaquamarine + group_skinbrown + group_skinred)/3),
  red_vs_other        = group_skinred - ((group_skinaquamarine + group_skinbrown + group_skingreen)/3),
  aquamarine_vs_brown = group_skinaquamarine - group_skinbrown,
  aquamarine_vs_green = group_skinaquamarine - group_skingreen,
  aquamarine_vs_red   = group_skinaquamarine - group_skinred,
  brown_vs_green      = group_skinbrown - group_skingreen,
  brown_vs_red        = group_skinbrown - group_skinred,
  green_vs_red        = group_skingreen - group_skinred,
  levels = colnames(coef(fit_skin))
)

fit_skin <- contrasts.fit(fit_skin, contrasts_skin)
fit_skin <- eBayes(fit_skin)

# Extract DE genes
top_skin <- topTable(fit_skin, p.value = 0.05, lfc = 1, sort.by = "F", n = Inf)
all_skin <- topTable(fit_skin, sort.by = "F", n = Inf)

write.csv(top_skin, "./limma_analysis_final/results/limma_all_contr_skin_significant.csv")
write.csv(all_skin, "./limma_analysis_final/results/limma_all_contr_skin_allresults.csv")

logCPM_skin <- cpm(dge_skin, log = TRUE, prior.count = 3)
write.csv(logCPM_skin, "vicentei_skin_logTMM_counts.csv")

###############################################################################
# Liver Data Analysis
###############################################################################

# Load liver metadata
liver <- read.csv("./metadata/vicentei_frogs_metadata_liver.csv") %>%
  filter(!row_number() %in% 10) %>%
  as.data.frame()
rownames(liver) <- liver$sample

# Load counts
tx2gene <- read.csv("transcript_to_gene_map.csv", stringsAsFactors = FALSE)
files_liver <- paste(liver$path, "abundance.tsv", sep = "/")
names(files_liver) <- rownames(liver)

txi_liver <- tximport(files_liver, type = "kallisto", tx2gene = tx2gene)
dge_liver <- DGEList(txi_liver$counts)

# Filter and normalize
group_liver <- interaction(liver$morph)
design_liver <- model.matrix(~0 + group_liver)
keep_liver <- filterByExpr(dge_liver, design_liver)
dge_liver <- dge_liver[keep_liver, ]
dge_liver <- calcNormFactors(dge_liver, method = "TMM")

# Voom transformation
v_liver <- voom(dge_liver, design_liver, plot = TRUE)
write.csv(v_liver$E, "./limma_analysis_final/results/voom_counts_liver_new.csv")

# Linear modeling
fit_liver <- lmFit(v_liver, design_liver)
contrasts_liver <- makeContrasts(
  aquamarine_vs_other = group_liveraquamarine - ((group_liverbrown + group_livergreen + group_liverred)/3),
  brown_vs_other      = group_liverbrown - ((group_liveraquamarine + group_livergreen + group_liverred)/3),
  green_vs_other      = group_livergreen - ((group_liveraquamarine + group_liverbrown + group_liverred)/3),
  red_vs_other        = group_liverred - ((group_liveraquamarine + group_liverbrown + group_livergreen)/3),
  aquamarine_vs_brown = group_liveraquamarine - group_liverbrown,
  aquamarine_vs_green = group_liveraquamarine - group_livergreen,
  aquamarine_vs_red   = group_liveraquamarine - group_liverred,
  brown_vs_green      = group_liverbrown - group_livergreen,
  brown_vs_red        = group_liverbrown - group_liverred,
  green_vs_red        = group_livergreen - group_liverred,
  levels = colnames(coef(fit_liver))
)

fit_liver <- contrasts.fit(fit_liver, contrasts_liver)
fit_liver <- eBayes(fit_liver)

# Extract DE genes
top_liver <- topTable(fit_liver, p.value = 0.05, lfc = 1, sort.by = "F", n = Inf)
all_liver <- topTable(fit_liver, sort.by = "F", n = Inf)

write.csv(top_liver, "./limma_analysis_final/results/limma_all_contr_liver_significant.csv")
write.csv(all_liver, "./limma_analysis_final/results/limma_all_contr_liver_allresults.csv")

# PCA plot
logCPM_liver <- cpm(dge_liver, log = TRUE, prior.count = 3)
write.csv(logCPM_liver, "./limma_analysis_final/results/vicentei_liver_logTMM_counts.csv")