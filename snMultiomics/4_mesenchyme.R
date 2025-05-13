library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(SeuratDisk)


inhouse<-readRDS("/vol/RZRZK/data/somi/seurat_d0.rds")

inhouse_mesen <- subset(inhouse, subset = main_celltype == "Mesenchyme")

inhouse_mesen <- NormalizeData(inhouse_mesen)
inhouse_mesen <- FindVariableFeatures(inhouse_mesen, selection.method = "vst", nfeatures = 2000)
inhouse_mesen <- ScaleData(inhouse_mesen)
inhouse_mesen <- RunPCA(inhouse_mesen, npcs = 50)
inhouse_mesen <- RunUMAP(inhouse_mesen, dims = 1:15, return.model = TRUE)
inhouse_mesen <- FindNeighbors(inhouse_mesen, dims = 1:15)
inhouse_mesen <- FindClusters(inhouse_mesen,resolution = 0.6)

inhouse_mesen$celltype <- inhouse_mesen$seurat_clusters
levels(inhouse_mesen$celltype)[match('0',levels(inhouse_mesen$celltype))] <- "Pdgfra high"
levels(inhouse_mesen$celltype)[match('1',levels(inhouse_mesen$celltype))] <- "Pdgfra low"
levels(inhouse_mesen$celltype)[match('2',levels(inhouse_mesen$celltype))] <- "Pdgfra high"
levels(inhouse_mesen$celltype)[match('3',levels(inhouse_mesen$celltype))] <- "Myh11+"
levels(inhouse_mesen$celltype)[match('4',levels(inhouse_mesen$celltype))] <- "Pdgfra low"
levels(inhouse_mesen$celltype)[match('5',levels(inhouse_mesen$celltype))] <- "Pdgfra high"
levels(inhouse_mesen$celltype)[match('6',levels(inhouse_mesen$celltype))] <- "Pdgfra low"
levels(inhouse_mesen$celltype)[match('7',levels(inhouse_mesen$celltype))] <- "Pdgfra low"
levels(inhouse_mesen$celltype)[match('8',levels(inhouse_mesen$celltype))] <- "Pdgfra high"
levels(inhouse_mesen$celltype)[match('9',levels(inhouse_mesen$celltype))] <- NA
levels(inhouse_mesen$celltype)[match('10',levels(inhouse_mesen$celltype))] <- NA
levels(inhouse_mesen$celltype)[match('11',levels(inhouse_mesen$celltype))] <- NA

####fig S1B
options(repr.plot.width= 10, repr.plot.height=10)
p1 <- DimPlot(inhouse_mesen, reduction = "umap", label = TRUE,group.by = "celltype",label.size = 7)
p1
ggsave(filename = "/vol/RZRZK/output/inhouse_mesenchyme_dimplot_celltype.pdf",plot =p1, dpi = 300, width = 10, height = 10)


####fig 1B
options(repr.plot.width= 30, repr.plot.height=20)

p2 <- FeaturePlot(inhouse_mesen, features = c("Wnt2b","Wnt4","Wnt5a","Wnt5b","Pdgfra","Myh11"), reduction = "umap", order = TRUE)
p2
ggsave(filename = "/vol/RZRZK/output/inhouse_mesenchyme_feature_celltype.pdf",plot =p2, dpi = 300, width = 30, height = 20)
