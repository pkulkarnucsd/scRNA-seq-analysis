# # Load required libraries
library(Seurat)
library(dplyr)
library(patchwork)
library(hdf5r)
library(ggplot2)
# library(tidyverse)


# args <- c("scRNA.balf.scRNA.healthy.feature_bc_matrix.h5", 
#           "scRNA.balf.scRNA.mild.feature_bc_matrix.h5", 
#           "scRNA.balf.scRNA.severe.feature_bc_matrix.h5")

# Parse input arguments
args <- commandArgs(trailingOnly = TRUE)
healthy_file <- args[1]
mild_file <- args[2]
severe_file <- args[3]
# 
# cat("Healthy file path: ", healthy_file, "\n")
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

# Normalize, find variable features, scale, PCA, UMAP, neighbors, clusters
balf <- NormalizeData(balf) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(dims = 1:30, reduction = "pca") %>%
  FindNeighbors(dims = 1:30, reduction = "pca") %>%
  FindClusters(resolution = 0.3)

# Save UMAP plots
umap_orig <- DimPlot(balf, reduction = "umap", group.by = "orig.ident", label = TRUE, repel = TRUE)
save_plot(umap_orig, "umap_origident.png")

umap_clusters <- DimPlot(balf, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)
save_plot(umap_clusters, "umap_clusters.png")

umap_split_by_treatment <- DimPlot(balf, reduction = "umap", split.by = "orig.ident")
save_plot(umap_split_by_treatment, "umap_split_by_treatment.png")

balf <- JoinLayers(balf)
# Find markers
balf.markers <- FindAllMarkers(balf, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Select top markers and save to CSV
bestMarkers <- balf.markers %>%
  group_by(cluster) %>%
  top_n(n = 1, wt = avg_log2FC) %>%
  arrange(desc(avg_log2FC))

# Save markers to a CSV
write.csv(bestMarkers, file = "best_markers.csv", row.names = FALSE)

# Save plots as PDF
pdf("output.pdf")
topDEG_genes <- unique(bestMarkers$gene)
for (gene in topDEG_genes) {
  p <- FeaturePlot(balf, features = gene, min.cutoff = "q9", split.by = "treatment") + 
    ggtitle(paste("Feature Plot for", gene))
  print(p)
}
dev.off()

