# # Load required libraries
##############################################
# Analyze Bronchoalveolar lavage fluid (BALF), commonly gathered during the diagnostic workup of pulmonary sarcoidosis.
# It is thought to contain the immune cells found in lung alveoli, so can provide important information regarding the 
# immunological response that takes place during pulmonary disease (like with Covid-19). 

#This script analyzes single cell RNA-seq from BALF samples from a healthy control, and patients with mild and severe Covid 19. 
##############################################

library(Seurat)
library(dplyr)
library(patchwork)
library(hdf5r)
library(ggplot2)

# Parse input arguments
args <- commandArgs(trailingOnly = TRUE)
healthy_file <- args[1]
mild_file <- args[2]
severe_file <- args[3]

#function to save plots
save_plot <- function(plot, filename) {
  ggsave(filename, plot, width = 10, height = 7)
}

# Load the healthy data
healthy.data <- Read10X_h5(healthy_file)
healthy <- CreateSeuratObject(counts = healthy.data, project = "healthy", min.cells = 3, min.features = 200)
cells <- Cells(healthy)
sample.info <- data.frame(treatment=rep("healthy",length(cells)),row.names=cells)
healthy <- AddMetaData(object = healthy, metadata = sample.info)
# only look at a subset of the healthy sample
x <- subset(healthy, downsample=3000)
healthy <- x

# Load the mild data
mild.data <- Read10X_h5(mild_file)
mild <- CreateSeuratObject(counts = mild.data, project = "mild", min.cells = 3, min.features = 200)
cells <- Cells(mild)
sample.info <- data.frame(treatment=rep("mild",length(cells)),row.names=cells)
mild <- AddMetaData(object = mild, metadata = sample.info)


# Load the severe data
severe.data <- Read10X_h5(severe_file)
severe <- CreateSeuratObject(counts = severe.data, project = "severe", min.cells = 3, min.features = 200)
cells <- Cells(severe)
sample.info <- data.frame(treatment=rep("severe",length(cells)),row.names=cells)
severe <- AddMetaData(object = severe, metadata = sample.info)

# Visualize each dataset and save plots
plots <- list(
  FeatureScatter(healthy, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"),
  VlnPlot(healthy, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2),
  FeatureScatter(mild, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"),
  VlnPlot(mild, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2),
  FeatureScatter(severe, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"),
  VlnPlot(severe, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
)

plot_filenames <- c(
  "healthy_featurescatter.png",
  "healthy_violinplot.png",
  "mild_featurescatter.png",
  "mild_violinplot.png",
  "severe_featurescatter.png",
  "severe_violinplot.png"
)

for (i in seq_along(plots)) {
  save_plot(plots[[i]], plot_filenames[i])
}

# Filter cells from each sample
healthy <- subset(healthy, subset = nFeature_RNA > 300 & nFeature_RNA < 3500)
mild <- subset(mild, subset = nFeature_RNA > 300 & nFeature_RNA < 3500)
severe <- subset(severe, subset = nFeature_RNA > 300 & nFeature_RNA < 3500)

# Merge datasets
balf <- merge(x = healthy, y = c(mild, severe), add.cell.ids = c("healthy", "mild", "severe"), project = "balf")

############################
# Seurat Pipeline
############################

# Normalize, find variable features, scale, PCA, UMAP, neighbors, clusters
balf <- NormalizeData(balf) %>%
  # Identify most variable genes across cells
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  # Standardize gene expression value 
  ScaleData(verbose = FALSE) %>%
  # Dimension Reduction
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(dims = 1:30, reduction = "pca") %>%
  # Computes Nearest neighbor graph 
  FindNeighbors(dims = 1:30, reduction = "pca") %>%
  # Clusters similar cells
  FindClusters(resolution = 0.3)

# Save UMAP plots
umap_orig <- DimPlot(balf, reduction = "umap", group.by = "orig.ident", label = TRUE, repel = TRUE)
save_plot(umap_orig, "umap_origident.png")

umap_clusters <- DimPlot(balf, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)
save_plot(umap_clusters, "umap_clusters.png")

umap_split_by_treatment <- DimPlot(balf, reduction = "umap", split.by = "orig.ident")
save_plot(umap_split_by_treatment, "umap_split_by_treatment.png")

options(repr.plot.width=10, repr.plot.height=7)

cellcounts <- table(Idents(balf.combined),balf.combined$treatment)
cellcountsnorm <- t(cellcounts)/colSums(cellcounts)
barplotCellCount <- barplot(cellcountsnorm,beside = TRUE, ylab="Fraction of cells in each Cluster")
save_plot(barplotCellCount, "barplotCellCount.png")
options(repr.plot.width=7, repr.plot.height=7)

balf <- JoinLayers(balf)
# Find markers
balf.markers <- FindAllMarkers(balf, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Select top markers and save to CSV
bestMarkers <- balf.markers %>%
  group_by(cluster) %>%
  top_n(n = 3, wt = avg_log2FC) %>%
  arrange(desc(avg_log2FC))

dotPlotGenes <- DotPlot(balf.combined, features =x, cols = c("green", "blue", "red"), dot.scale = 5, split.by = "treatment") + RotatedAxis()
save_plot(dotPlotGenes, "dotPlotGenes.png")

# Save markers to a CSV
write.csv(bestMarkers, file = "best_markers.csv", row.names = FALSE)

topDEG_genes <- balf.markers %>%
  group_by(cluster) %>%
  top_n(n = 1, wt = avg_log2FC) %>%
  arrange(desc(avg_log2FC))
# Save plots as PDF
pdf("output.pdf")
topDEG_genes <- unique(topDEG_genes$gene)
for (gene in topDEG_genes) {
  p <- FeaturePlot(balf, features = gene, min.cutoff = "q9", split.by = "treatment") + 
    ggtitle(paste("Feature Plot for", gene))
  print(p)
}
dev.off()

