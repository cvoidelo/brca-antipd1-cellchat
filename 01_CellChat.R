# Script 1. Cell-cell communication inference per clinical-molecular subgroup (On-treatment)
# Dataset: BioKey Clinical Study (NCT03197389; Bassez et al., 2021)
# Author: Clayton Voidelo

# Adapted from the official CellChat v2 tutorial:
# https://github.com/jinworks/CellChat
# https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html

# Run this script for each of the four on-treatment subgroups by setting BC_TYPE and EXPANSION:
#   1) ER+/E   (BC_TYPE = "ER",   EXPANSION = "E")
#   2) ER+/NE  (BC_TYPE = "ER",   EXPANSION = "NE")
#   3) TNBC/E  (BC_TYPE = "TNBC", EXPANSION = "E")
#   4) TNBC/NE (BC_TYPE = "TNBC", EXPANSION = "NE")

#---

#packages
library(Seurat)
library(CellChat)
library(patchwork)
library(NMF)
library(ggalluvial)
options(stringsAsFactors = FALSE)

# Subgroup parameters
BC_TYPE   <- "TNBC"  # "ER" or "TNBC"
EXPANSION <- "NE"    # "E" or "NE"
SUB_DIR   <- ifelse(BC_TYPE == "ER", "ER+", "TNBC")

output_dir <- file.path("results", SUB_DIR, EXPANSION)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load raw counts and metadata (Bassez et al., 2021; BioKey cohort)
data.input.raw <- readRDS("data/1863-counts_cells_cohort1.rds")
meta <- read.csv("data/1872-BIOKEY_metaData_cohort1_web.csv", row.names = 1)

# Define cell type labels and subset on-treatment cells for the selected subgroup
meta$labels <- meta$cellType
cell.use <- rownames(meta)[meta$timepoint == "On" & meta$expansion == EXPANSION & meta$BC_type == BC_TYPE]

data.input.raw <- data.input.raw[, cell.use]
meta <- meta[cell.use, ]

# Create Seurat object and normalize (log-normalization, scale factor 10,000)
seurat.obj <- CreateSeuratObject(counts = data.input.raw, meta.data = meta)
seurat.obj <- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)

# Optional Seurat steps (used for UMAP visualization in Fig. 1B)
seurat.obj <- FindVariableFeatures(seurat.obj)
seurat.obj <- ScaleData(seurat.obj)
seurat.obj <- RunPCA(seurat.obj, verbose = FALSE)
seurat.obj <- RunUMAP(seurat.obj, dims = 1:20, verbose = FALSE)
DimPlot(seurat.obj, reduction = "umap", group.by = "labels") # Fig. 1B

# Extract log-normalized expression matrix and initialize CellChat object
data.input <- GetAssayData(seurat.obj, layer = "data")
meta$samples <- NULL
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

# Set ligand-receptor database (CellChatDB.human; use CellChatDB.mouse for mouse data)
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)

# Restrict to "Secreted Signaling" to evaluate paracrine and autocrine crosstalk
# Alternative: subsetDB(CellChatDB) includes all categories except "Non-protein Signaling"
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation")
cellchat@DB <- CellChatDB.use

# Preprocessing expression data for communication inference
cellchat <- subsetData(cellchat)
options(future.globals.maxSize = 10 * 1024^3)
future::plan("sequential") # Can be set to "multisession" with workers > 1 for parallel execution

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Optional: project gene expression onto human PPI network (requires raw.use = FALSE in computeCommunProb)
# cellchat <- projectData(cellchat, PPI.human)

#---
# Part II: Inference of cell-cell communication network
#---

# Compute communication probabilities using truncated mean (triMean) and filter groups < 10 cells
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Export inferred L-R interactions (use slot.name = "netP" for pathway-level table)
df.net <- subsetCommunication(cellchat)
write.csv(df.net, file = file.path(output_dir, paste0(SUB_DIR, "_", EXPANSION, "_cellchat_interactions.csv")), row.names = FALSE)

# Infer communication at the signaling pathway level and aggregate network
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# Circle plots of aggregated interaction count and strength
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction weights/strength")

#---
# Part III: Visualization of cell-cell communication networks
#---

# All significant signaling pathways in the subgroup: cellchat@netP$pathways
# Examples of pathways visualized in main and supplementary figures:
#   - Fig. 3B: exclusive pathways in ER+/E vs. ER+/NE
#   - Fig. 4B–D: "SPP1" (ER+/NE) and "CXCL" (ER+/E)
#   - Fig. 6A, 6C–F: "CXCL" and "GALECTIN" (TNBC/E); "CXCL", "GALECTIN", and "IGF" (TNBC/NE)
#   - Supp. Fig. 4B: exclusive pathways in TNBC/E vs. TNBC/NE
#   - Supp. Fig. 5: "MK" pathway (ER+/E and TNBC/E)

pathways.show <- c("CXCL")

# Circle plot, chord diagram, and heatmap for a selected pathway
par(mfrow = c(1, 1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord") # Fig. 3B, Fig. 6A, Supp. Fig. 4B, Supp. Fig. 5
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

# Relative contribution of each L-R pair to the selected pathway (Fig. 4B, Fig. 6C)
netAnalysis_contribution(cellchat, signaling = pathways.show)

# Bubble plots of significant L-R interactions from specific sender populations to target populations
# indices in sources.use / targets.use were adjusted to the cell populations evaluated in each figure
# (e.g., Supp. Fig. 9: cancer cells -> all populations in ER+/E)
netVisual_bubble(cellchat, sources.use = 1, targets.use = seq_along(levels(cellchat@idents)), remove.isolate = FALSE)

# Optional: restrict bubble plot to specific pathways or L-R pairs
# pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CXCL", "SPP1"))
# netVisual_bubble(cellchat, sources.use = 1, targets.use = seq_along(levels(cellchat@idents)), pairLR.use = pairLR.use, remove.isolate = TRUE)

#---
# Part IV: Systems analysis of cell-cell communication networks
#---

# Compute network centrality scores and export centrality heatmaps + chord diagrams for all active pathways
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
pathways.show.all <- cellchat@netP$pathways

pathway_plots_dir <- file.path(output_dir, "pathways_centrality_and_chord")
dir.create(pathway_plots_dir, recursive = TRUE, showWarnings = FALSE)

for (pathway_name in pathways.show.all) {
  # Centrality role heatmap: sender, receiver, mediator, influencer (Fig. 4C–D, Fig. 6D–F)
  pdf(file = file.path(pathway_plots_dir, paste0(pathway_name, "_centrality_heatmap.pdf")), width = 8, height = 2.5)
  netAnalysis_signalingRole_network(cellchat, signaling = pathway_name, width = 8, height = 2.5, font.size = 10)
  dev.off()
  
  # Pathway-level chord diagram (Fig. 3B, Fig. 6A, Supp. Fig. 4B, Supp. Fig. 5)
  cellchat_sub <- subsetCommunication(cellchat, signaling = pathway_name, slot.name = "netP")
  plot_height <- max(6, min(12, 0.25 * nrow(cellchat_sub)))
  
  pdf(file = file.path(pathway_plots_dir, paste0(pathway_name, "_chord_diagram.pdf")), width = 8, height = plot_height)
  netVisual_aggregate(cellchat, signaling = pathway_name, layout = "chord", slot.name = "netP")
  dev.off()
}

# Dominant senders and receivers in 2D space (single-dataset overview)
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("CXCL", "SPP1"))
gg1 + gg2

# Global communication patterns via Non-negative Matrix Factorization (NMF)
# 1) Outgoing patterns (secreting cells)
selectK(cellchat, pattern = "outgoing") # Supp. Fig. 10 (right panels)

# Set nPatterns (k) based on the Cophenetic and Silhouette drop observed in selectK for each subgroup
nPatterns_out <- 2
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns_out) # Supp. Fig. 11 & 12 (outgoing heatmaps)
netAnalysis_river(cellchat, pattern = "outgoing") # Fig. 7C (bottom panels for ER+/E and TNBC/E)

# 2) Incoming patterns (target cells)
selectK(cellchat, pattern = "incoming") # Supp. Fig. 10 (left panels)

nPatterns_in <- 2
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns_in) # Supp. Fig. 11 & 12 (incoming heatmaps)
netAnalysis_river(cellchat, pattern = "incoming") # Fig. 7C (top panels for ER+/E and TNBC/E)

# Save subgroup CellChat object for downstream comparative analysis (Script 02)
saveRDS(cellchat, file = file.path(output_dir, paste0("cellchat_BRCA_Pembro_On_", BC_TYPE, "_", EXPANSION, ".rds")))

sessionInfo()