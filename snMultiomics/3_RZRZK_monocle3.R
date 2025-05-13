library(Seurat)
library(SingleCellExperiment)
library(Signac)
library(ggplot2)
library(dittoSeq)
library(RColorBrewer)
library(paletteer)
library(patchwork)
library(DESeq2)
library(monocle3)
set.seed(1234)

image.name = "/hdds/hdd1/3_RZK_from_IMBA/3_RZK/0_FIN_scMultiome/RZ_RZK_object_v20220216.RData"
load(image.name)
rm(list = setdiff(ls(), "ready_to_visualize"))

my_col <- c("0_Wild_type"='#8936EF',"0_Mutant"='#8936EF',"1_Wild_type"='#F2CA19',"1_Mutant"='#F2CA19',"2_Wild_type"='#FF00BD',"2_Mutant"='#FF00BD',"3_Wild_type"='#E11845',"3_Mutant"='#E11845',"4_Wild_type"='#0057E9',"4_Mutant"='#0057E9',"5_Wild_type"='#87E911',"5_Mutant"='#87E911',"6_Wild_type"="#018300","6_Mutant"="#018300")
cluster_col <- c("Neck"='#8936EF',"Lgr5+"='#F2CA19',"SPEM"='#FF00BD',"Wnt7+"='#E11845',"Proliferating"='#0057E9',"Pre-Pit"='#87E911',"Pit"="#018300")

#clustering more to split several pit cell types
DefaultAssay(ready_to_visualize) <- "integrated"
options(repr.plot.width=6, repr.plot.height=6)
RZ.RZK <- FindClusters(ready_to_visualize, resolution = 0.6)
DimPlot(RZ.RZK, group.by = "seurat_clusters")

RZ.RZK@meta.data$cl_genotype[RZ.RZK@meta.data$seurat_clusters == "9" & RZ.RZK@meta.data$kras_genotype == "Wild_type"] <- "6_Wild_type"
RZ.RZK@meta.data$cl_genotype[RZ.RZK@meta.data$seurat_clusters == "9" & RZ.RZK@meta.data$kras_genotype == "Mutant"] <- "6_Mutant"

RZ.RZK@meta.data$clusters[RZ.RZK@meta.data$cl_genotype == "0_Wild_type"] <- "Neck"
RZ.RZK@meta.data$clusters[RZ.RZK@meta.data$cl_genotype == "0_Mutant"] <- "Neck"
RZ.RZK@meta.data$clusters[RZ.RZK@meta.data$cl_genotype == "1_Wild_type"] <- "Lgr5+"
RZ.RZK@meta.data$clusters[RZ.RZK@meta.data$cl_genotype == "1_Mutant"] <- "Lgr5+"
RZ.RZK@meta.data$clusters[RZ.RZK@meta.data$cl_genotype == "2_Wild_type"] <- "SPEM"
RZ.RZK@meta.data$clusters[RZ.RZK@meta.data$cl_genotype == "2_Mutant"] <- "SPEM"
RZ.RZK@meta.data$clusters[RZ.RZK@meta.data$cl_genotype == "3_Wild_type"] <- "Wnt7+"
RZ.RZK@meta.data$clusters[RZ.RZK@meta.data$cl_genotype == "3_Mutant"] <- "Wnt7+"
RZ.RZK@meta.data$clusters[RZ.RZK@meta.data$cl_genotype == "4_Wild_type"] <- "Proliferating"
RZ.RZK@meta.data$clusters[RZ.RZK@meta.data$cl_genotype == "4_Mutant"] <- "Proliferating"
RZ.RZK@meta.data$clusters[RZ.RZK@meta.data$cl_genotype == "5_Wild_type"] <- "Pre-Pit"
RZ.RZK@meta.data$clusters[RZ.RZK@meta.data$cl_genotype == "5_Mutant"] <- "Pre-Pit"
RZ.RZK@meta.data$clusters[RZ.RZK@meta.data$cl_genotype == "6_Wild_type"] <- "Pit"
RZ.RZK@meta.data$clusters[RZ.RZK@meta.data$cl_genotype == "6_Mutant"] <- "Pit"

rm("ready_to_visualize")

#Extract QC passed cells
rz.cell <- rownames(subset(RZ.RZK@meta.data, kras_genotype=='Wild_type'))
rz.prlf.cell <- rownames(subset(RZ.RZK@meta.data, cl_genotype=='4_Wild_type'))
rz.non.prlf.cell <- rz.cell[!rz.cell %in% rz.prlf.cell]
using.cell <- rz.cell
using.cell.non.prlf <- rz.non.prlf.cell
#
rzk.cell <- rownames(subset(RZ.RZK@meta.data, kras_genotype=='Mutant'))
rzk.prlf.cell <- rownames(subset(RZ.RZK@meta.data, cl_genotype=='4_Mutant'))
rzk.non.prlf.cell <- rzk.cell[!rzk.cell %in% rzk.prlf.cell]
using.cell <- append(using.cell, rzk.cell)
using.cell.non.prlf <- append(using.cell.non.prlf, rzk.non.prlf.cell)

RZRZK <- subset(RZ.RZK, subset = cl_genotype != "4_Mutant" & cl_genotype != "4_Wild_type")

##RZ, RZK counts in SCT monocle3

DefaultAssay(RZRZK) <- 'SCT'
RZRZK_m3 <- as.SingleCellExperiment(RZRZK, assay = "SCT")
rowData(RZRZK_m3)$gene_short_name <- rownames(RZRZK)
#Now create the Monocle object using the raw counts
cds_corrected_m3 <- new_cell_data_set(expression_data = counts(RZRZK_m3), 
                               cell_metadata = as.data.frame(colData(RZRZK_m3)),
                               gene_metadata = as.data.frame(rowData(RZRZK_m3)))
cds_corrected_m3 <- estimate_size_factors(cds_corrected_m3)

#Continue with the Monocle3 analysis pipeline
cds_corrected_m3 <- preprocess_cds(cds_corrected_m3, num_dim = 100, norm_method = "log")
cds_corrected_m3 <- align_cds(cds_corrected_m3, alignment_group='kras_genotype')
cds_corrected_m3 <- reduce_dimension(cds_corrected_m3, reduction_method = "UMAP")
cds_corrected_m3 <- cluster_cells(cds_corrected_m3, resolution = 1e-3,random_seed = 8888)
cds_corrected_m3 <- learn_graph(cds_corrected_m3)

#### fig 3f
options(repr.plot.width=12, repr.plot.height=12)
p1 <- plot_cells(cds_corrected_m3,
           color_cells_by = "clusters",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_principal_points=FALSE,
           graph_label_size=5,
           cell_size = .75,
           trajectory_graph_segment_size = 2,
           trajectory_graph_color = "black") + scale_color_manual(values = cluster_col, name = "clusters")
p1

options(repr.plot.width=12, repr.plot.height=12)
pData(cds_corrected_m3)$monocluster <- clusters(cds_corrected_m3)

p2 <- plot_cells(cds_corrected_m3,
           color_cells_by = "monocluster",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           label_principal_points=TRUE,
           graph_label_size=5,
           cell_size = .75)
p2

cds_corrected_m3$cluster_mono[cds_corrected_m3$monocluster == "1"] <- "Neck"
cds_corrected_m3$cluster_mono[cds_corrected_m3$monocluster == "2"] <- "Pre-Pit"
cds_corrected_m3$cluster_mono[cds_corrected_m3$monocluster == "3"] <- "Wnt7+"
cds_corrected_m3$cluster_mono[cds_corrected_m3$monocluster == "4"] <- "SPEM"
cds_corrected_m3$cluster_mono[cds_corrected_m3$monocluster == "5"] <- "Lgr5+"
cds_corrected_m3$cluster_mono[cds_corrected_m3$monocluster == "6"] <- "Wnt7+"
cds_corrected_m3$cluster_mono[cds_corrected_m3$monocluster == "7"] <- "Neck"
cds_corrected_m3$cluster_mono[cds_corrected_m3$monocluster == "8"] <- "Lgr5+"
cds_corrected_m3$cluster_mono[cds_corrected_m3$monocluster == "9"] <- "Neck"
cds_corrected_m3$cluster_mono[cds_corrected_m3$monocluster == "10"] <- "SPEM"
cds_corrected_m3$cluster_mono[cds_corrected_m3$monocluster == "11"] <- "Lgr5+"
cds_corrected_m3$cluster_mono[cds_corrected_m3$monocluster == "12"] <- "Neck"
cds_corrected_m3$cluster_mono[cds_corrected_m3$monocluster == "13"] <- "Pre-Pit"
cds_corrected_m3$cluster_mono[cds_corrected_m3$monocluster == "14"] <- "Pit"

options(repr.plot.width=12, repr.plot.height=12)

p1.total <- plot_cells(cds_corrected_m3,
           color_cells_by = "cluster_mono",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_principal_points=FALSE,
           graph_label_size=5,
           cell_size = .75,
           trajectory_graph_segment_size = 2,
           trajectory_graph_color = "black") + scale_color_manual(values = cluster_col, name = "clusters") 

ggsave(file="/hdds/hdd1/3_RZK_from_IMBA/youngchul/pdf/fig3f_monocle3_RZ_RZK_p1.pdf", width=10, height=10, dpi=300, plot=p1.total)

cds_corrected_m3 <- order_cells(cds_corrected_m3, root_pr_nodes=c('Y_32'))

spectral_colors <- colorRampPalette(rev(brewer.pal(11, "Spectral")))(100)

options(repr.plot.width=10, repr.plot.height=10)
p1.pseudo <- plot_cells(cds_corrected_m3,           
          color_cells_by = "pseudotime",           
           label_cell_groups=TRUE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=8,
           cell_size = 1,
           trajectory_graph_segment_size = 2,
           trajectory_graph_color = "black") + scale_color_gradientn(colors = spectral_colors) + labs(color = "Pseudotime")
ggsave(file="/hdds/hdd1/3_RZK_from_IMBA/youngchul/pdf/fig3f_monocle3_pseudotime.pdf", width=10, height=10, dpi=300, plot=p1.pseudo)

#### fig 3g
# log10(expression + min_expr)
log_expr <- normalized_counts(cds_corrected_m3)

expr_raw <- log_expr["Wnt7b", ]
expr_percent <- expr_raw / max(expr_raw, na.rm = TRUE) * 100

colData(cds_corrected_m3)$Wnt7b_percent <- expr_percent

p2.wnt <- plot_cells(cds_corrected_m3[, rz.non.prlf.cell],
           color_cells_by = "Wnt7b_percent",   
           label_cell_groups = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           label_principal_points = FALSE,
           graph_label_size = 5,
           cell_size = 3,
           trajectory_graph_segment_size = 2,
           cell_stroke = 0,01,
           trajectory_graph_color = "black") +
  scale_color_gradientn(colors = ht_custom_col, limits = c(0,100)) +
  labs(color = "Wnt7b (%)")+
theme(panel.border = element_blank())

p3.wnt <- plot_cells(cds_corrected_m3[, rzk.non.prlf.cell],
           color_cells_by = "Wnt7b_percent",   
           label_cell_groups = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           label_principal_points = FALSE,
           graph_label_size = 5,
           cell_size = 3,
           cell_stroke = 0,01,
           trajectory_graph_segment_size = 2,
           trajectory_graph_color = "black") +
  scale_color_gradientn(colors = ht_custom_col, limits = c(0,100)) +
  labs(color = "Wnt7b (%)")+
theme(panel.border = element_blank())

ggsave(file="/hdds/hdd1/3_RZK_from_IMBA/youngchul/pdf/fig3f_monocle3_RZ_RZK_p6_250414.pdf", width=10, height=10, dpi=300, plot=p2.wnt)
ggsave(file="/hdds/hdd1/3_RZK_from_IMBA/youngchul/pdf/fig3f_monocle3_RZ_RZK_p9_250414.pdf", width=10, height=10, dpi=300, plot=p3.wnt)