library(future)
library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(EnsDb.Mmusculus.v79)
library(EnsDb.Hsapiens.v86)
library(RColorBrewer)
library(clustree)
library(hash)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmusculus.UCSC.mm10)
library(dittoSeq)
library(motifmatchr)
library(JASPAR2020)
library(TFBSTools)
library(chromVAR)
library(httpgd)

##5-1 Customized Feature Plots
plot_featureplot <- function(obj, gene){
  DefaultAssay(obj) <- "RNA"
  p11 <- FeaturePlot(obj, features = gene, pt.size = 1, min.cutoff = 0.3, max.cutoff = 2)
  ggplot(data = p11$data, aes_string(x = 'UMAP_1', y = 'UMAP_2', fill = gene)) +    
    geom_point(shape = 21, stroke = 0.3, color = "black", alpha = 1, size = 2) +
    scale_fill_gradientn(colours = brewer.pal(n = 9, name = "YlOrRd")) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
    #labs(title = gene, subtitle = "Expression") +
} 
plot_featureplot_add_max <- function(obj, gene, max){
  DefaultAssay(obj) <- "RNA"
  #ht_custom_col <- colorRampPalette(colors = c("#FEF9E7", "#A93226"))(13)
  ht_custom_col <- brewer.pal(n = 9, name = "YlOrRd")
  p11 <- FeaturePlot(obj, features = gene, pt.size = 1, min.cutoff = 0.0, max.cutoff = max)
  ggplot(data = p11$data, aes_string(x = 'UMAP_1', y = 'UMAP_2', fill = gene)) +    
    geom_point(shape = 21, stroke = 0.3, color = "black", alpha = 1, size = 2) +
    scale_fill_gradientn(colours = ht_custom_col, limits=c(0.0, max)) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
    #labs(title = gene, subtitle = "Expression") +
} 
#if reverse: scale_fill_gradientn(colours = rev(brewer.pal(n = 9, name = "YlOrRd")))

enrichment_test_for_counts <- function(count_region, ident1, ident2){
    target_cells_1 <- colnames(subset(ready_to_visualize, subset=cl_genotype==ident1))
    target_cells_2 <- colnames(subset(ready_to_visualize, subset=cl_genotype==ident2))
    count_target_1 <- count_region[target_cells_1]
    count_target_2 <- count_region[target_cells_2]
    result <- t.test(count_target_1, count_target_2)
    print(c(result$p.value, result$estimate[1], result$estimate[2]))
    }
wilcoxon_for_expression <- function(ident1, ident2){
    target_cells_1 <- colnames(subset(RZ.RZK, subset=cl_genotype==ident1))
    target_cells_2 <- colnames(subset(RZ.RZK, subset=cl_genotype==ident2))
    count_target_1 <- RZ.RZK[["RNA"]]@data["Wnt7b",target_cells_1]
    count_target_2 <- RZ.RZK[["RNA"]]@data["Wnt7b",target_cells_2]
    result <- wilcox.test(count_target_1, count_target_2)
    print(c(result$p.value, result$estimate[1], result$estimate[2]))
    }

plot_featureplot_v3 <- function(obj, gene, min, max){
  DefaultAssay(obj) <- "RNA"
  ht_custom_col <- rev(paletteer_c("grDevices::Spectral", 30))
  p11 <- FeaturePlot(obj, features = gene, pt.size = 1, min.cutoff = min, max.cutoff = max, order = TRUE)
  ggplot(data = p11$data, aes_string(x = 'UMAP_1', y = 'UMAP_2', fill = gene)) +    
    geom_point(shape = 21, stroke = NA, color = "black", alpha = 1, size = 2) +
    scale_fill_gradientn(colours = ht_custom_col, name = "Expression\nLevel") +
    ggtitle(gene) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 16), plot.subtitle = element_text(hjust = 0.5))
}

plot_linkage_cl_genotype <- function(obj, assay, Target_gene){
    DefaultAssay(obj) <- assay
    # first compute the GC content for each peak
    obj <- RegionStats(obj, genome = linkage_genome)
    # link peaks to genes
    obj <- LinkPeaks(
    object = obj,
    peak.assay = assay,
    expression.assay = "RNA",
    genes.use = c(Target_gene)
    )
    lv.list <- c()
    for (i in c(1,2,0,4,5,6,3)){
        lv.list <- append(lv.list,sprintf("%s_Wild_type",i))
        lv.list <- append(lv.list,sprintf("%s_Mutant",i))    
    }
    obj$cl_genotype <- factor(obj$cl_genotype, levels=lv.list)
    CoveragePlot(
    object = obj,
    region = Target_gene,
    features = Target_gene,
    expression.assay = "RNA",  
    group.by="cl_genotype",
    extend.upstream = 0,   
    extend.downstream = 0
    )
} 

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

####fig 3A
options(repr.plot.width=8, repr.plot.height=9)
p1 <- DimPlot(RZ.RZK, group.by = "cl_genotype", cols = my_col, pt.size = 1.3) + NoLegend()
ggsave(p1, filename = "/hdds/hdd1/3_RZK_from_IMBA/youngchul/pdf/fig_3_A_UMAP.pdf",width = 8, height = 9)
p1
dev.off()

####fig 3B
options(repr.plot.width=6, repr.plot.height=7)
DefaultAssay(RZ.RZK) <- "RNA"
mylevel <- c("Lgr5+", "SPEM","Neck","Proliferating","Pre-Pit","Pit", "Wnt7+")
Idents(RZ.RZK) <- RZ.RZK@meta.data$clusters
Idents(RZ.RZK) <- factor(Idents(RZ.RZK), levels= mylevel)
myGenes = c("Lgr5", "Glipr1", "Cd44", "Cftr", "Muc6", "Mki67", "Foxm1","Hmgb2","Top2a","Smc2","Gkn1","Tff1","Gkn2","Muc5ac", "Wnt7b")
p2<-DotPlot(RZ.RZK, features=myGenes, dot.scale = 10) + scale_colour_gradient2(low="dodgerblue3", mid="ghostwhite", high="firebrick",limits = c(-2.5,2.5)) + xlab('Markers') +  ylab('Clusters') + coord_flip() + theme(axis.text=element_text(size=15),legend.position="top")
ggsave(p2, filename = "/hdds/hdd1/3_RZK_from_IMBA/youngchul/pdf/fig_3_B_dotplot_2.pdf",width = 6, height = 7)

####fig 3C
options(repr.plot.width=6, repr.plot.height=6)
mylevel <- c("Lgr5+", "SPEM","Neck","Proliferating","Pit","Pre-Pit", "Wnt7+")
Idents(RZ.RZK) <- RZ.RZK@meta.data$clusters
Idents(RZ.RZK) <- factor(Idents(RZ.RZK), levels= mylevel)
p3 <- dittoBarPlot(RZ.RZK, Idents(RZ.RZK), group.by = "kras_genotype", x.reorder=c(2,1), var.labels.reorder=c(1,6,2,5,4,3,7), color.panel = cluster_col)
ggsave(p3, filename = "/hdds/hdd1/3_RZK_from_IMBA/youngchul/pdf/fig_3_C_dittobarplot.pdf",width = 6, height = 6)

####fig 3D
options(repr.plot.height = 6, repr.plot.width = 12)
p1_0c <- plot_featureplot(subset(ready_to_visualize, subset=kras_genotype=="Wild_type"), "Wnt7b")
p2_0c <- plot_featureplot(subset(ready_to_visualize, subset=kras_genotype=="Mutant"), "Wnt7b")
p4 <- p1_0c+p2_0c
ggsave(file="/hdds/hdd1/3_RZK_from_IMBA/youngchul/pdf/fig_3_D_Wnt7b_feature.pdf", width=12, height=6, dpi=300, plot=p4)
p4
dev.off()

####fig 3E
cl_vec<- c('0_Wild_type','0_Mutant','1_Wild_type','1_Mutant','2_Wild_type','2_Mutant','3_Wild_type','3_Mutant','4_Wild_type','4_Mutant','5_Wild_type','5_Mutant','6_Wild_type','6_Mutant')

df <- data.frame(Clusters = cl_vec, Percentage = pp)
options(repr.plot.height = 8, repr.plot.width = 10)
p5 <- ggplot(df, aes(Clusters, Percentage, fill = Clusters)) + 
        geom_bar(stat="identity") + 
        scale_x_discrete(limits = c('1_Wild_type','1_Mutant','2_Wild_type','2_Mutant','0_Wild_type','0_Mutant','4_Wild_type','4_Mutant','5_Wild_type','5_Mutant','6_Wild_type','6_Mutant','3_Wild_type','3_Mutant'))+ 
        ylim(0,100) +
        theme_classic()+ 
        scale_fill_manual(values=my_col) + 
        ggtitle("Wnt7b-expressing cells (%)") + 
        theme(plot.title = element_text(hjust = 0.5))+
        NoLegend()

ggsave(file="/hdds/hdd1/3_RZK_from_IMBA/youngchul/pdf/fig_3_E_Wnt7b_expression cells.pdf", width=10, height=6, dpi=300, plot=p5)

options(repr.plot.width=15, repr.plot.height=3)

wnt7b.RZ <- RZ.RZK[["RNA"]]@data["Wnt7b",] * (RZ.RZK$kras_genotype == "Wild_type")
wnt7b.RZK <- RZ.RZK[["RNA"]]@data["Wnt7b",] * (RZ.RZK$kras_genotype == "Mutant")

RZ.RZK[["NEW"]] <- CreateAssayObject(data = rbind(wnt7b.RZ,wnt7b.RZK),)
DefaultAssay(RZ.RZK) <- "NEW"

Idents(RZ.RZK) <- RZ.RZK@meta.data$clusters
Idents(RZ.RZK) <- factor(Idents(RZ.RZK), levels= mylevel)
p2<-DotPlot(RZ.RZK, features=c("wnt7b.RZK","wnt7b.RZ"), dot.scale = 10) + scale_colour_gradient2(low="dodgerblue3", mid="ghostwhite", high="firebrick",limits = c(-2.5,2.5)) + xlab('Markers') +  ylab('Clusters') + coord_flip() + theme(axis.text=element_text(size=15),legend.position="top")
p2

ggsave(p2, filename = "/hdds/hdd1/3_RZK_from_IMBA/youngchul/pdf/fig3e_dotplot_wnt7b_RZRZK.pdf",width = 15, height = 3)


####fig 3E statistics
cl_list <- list(c('0_Wild_type','0_Mutant'),c('1_Wild_type','1_Mutant'),c('2_Wild_type','2_Mutant'),c('3_Wild_type','3_Mutant'),c('4_Wild_type','4_Mutant'),c('5_Wild_type','5_Mutant'),c('6_Wild_type','6_Mutant'))
pp <- c()
length_list <- c()
gc_list <- c()
for (i in cl_list){
    per <- c()
    wnt <- c()
    minwnt <- c()
    for (j in i){
        set <- subset(RZ.RZK, subset = cl_genotype == j)
        count <- set[["RNA"]]@data["Wnt7b",]
        length <- length(count)
        gc <- sum(count != 0)
        percent <- gc/length*100
        per <- c(per, percent)
        length_list <- c(length_list, length-gc)
        gc_list <- c(gc_list,gc)
        wnt <- c(wnt,gc)
        minwnt <- c(minwnt,length-gc)
    }
    fisher.df = matrix(c(wnt, minwnt), nrow=2)
    fisher_p = fisher.test(fisher.df, alternative="less")$p.value
    fisher_stat = fisher.test(fisher.df, alternative="less")$estimate
    pp <- c(pp, per)
    out_str <- paste(i, wnt[1], minwnt[1], wnt[2], minwnt[2], fisher_p, fisher_stat, sep=',')
    print(out_str)
}

####fig 3H
options(repr.plot.height = 8, repr.plot.width = 8)
p6 <- plot_linkage_cl_genotype(RZ.RZK, 'macs2', "Wnt7b")
p6
ggsave(file="/hdds/hdd1/3_RZK_from_IMBA/youngchul/pdf/fig_3_G_linkage_coverage_plot.pdf", width=8, height=8, dpi=300, plot=p6)

####fig 3H statistics
wilcoxon_for_expression("0_Wild_type","0_Mutant")
wilcoxon_for_expression("1_Wild_type","1_Mutant")
wilcoxon_for_expression("2_Wild_type","2_Mutant")
wilcoxon_for_expression("3_Wild_type","3_Mutant")
wilcoxon_for_expression("4_Wild_type","4_Mutant")
wilcoxon_for_expression("5_Wild_type","5_Mutant")
wilcoxon_for_expression("6_Wild_type","6_Mutant")

count_wnt7b_region1 <- CountsInRegion(
  object = RZ.RZK,
  assay = 'macs2',
  regions = StringToGRanges('chr15-85568400-85569900')
)
enrichment_test_for_counts(count_wnt7b_region1,'0_Mutant', '0_Wild_type')
enrichment_test_for_counts(count_wnt7b_region1,'1_Mutant', '1_Wild_type')
enrichment_test_for_counts(count_wnt7b_region1,'2_Mutant', '2_Wild_type')
enrichment_test_for_counts(count_wnt7b_region1,'3_Mutant', '3_Wild_type')
enrichment_test_for_counts(count_wnt7b_region1,'4_Mutant', '4_Wild_type')
enrichment_test_for_counts(count_wnt7b_region1,'5_Mutant', '5_Wild_type')
enrichment_test_for_counts(count_wnt7b_region1,'6_Mutant', '6_Wild_type')

count_wnt7b_region2 <- CountsInRegion(
  object = RZ.RZK,
  assay = 'macs2',
  regions = StringToGRanges('chr15-85574150-85575200')
)
enrichment_test_for_counts(count_wnt7b_region2,'0_Mutant', '0_Wild_type')
enrichment_test_for_counts(count_wnt7b_region2,'1_Mutant', '1_Wild_type')
enrichment_test_for_counts(count_wnt7b_region2,'2_Mutant', '2_Wild_type')
enrichment_test_for_counts(count_wnt7b_region2,'3_Mutant', '3_Wild_type')
enrichment_test_for_counts(count_wnt7b_region2,'4_Mutant', '4_Wild_type')
enrichment_test_for_counts(count_wnt7b_region2,'5_Mutant', '5_Wild_type')
enrichment_test_for_counts(count_wnt7b_region2,'6_Mutant', '6_Wild_type')

####fig 3I, S3C
DefaultAssay(ready_to_visualize) <- "ATAC"
plan("multicore", workers = 1)
da_peaks.31.m <- FindMarkers(ready_to_visualize, ident.1 = "3_Mutant", ident.2 = "1_Mutant", group.by="cl_genotype", min.pct = 0.05,  logfc.threshold = 0.1, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
da_peaks.31.w <- FindMarkers(ready_to_visualize, ident.1 = "3_Wild_type", ident.2 = "1_Wild_type", group.by="cl_genotype", min.pct = 0.05,  logfc.threshold = 0.1, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')

### Wnt7 vs Lgr5 in mutant
closest_genes <- ClosestFeature(ready_to_visualize, regions = rownames(da_peaks.31.m))
rownames(closest_genes) <- closest_genes$query_region
da_peaks.31.m$gene <- closest_genes$gene_name
da_peaks.31.m$gene_biotype <- closest_genes$gene_biotype
da_peaks.31.m$closest_region <- closest_genes$closest_region
da_peaks.31.m$distance <- closest_genes$distance
wnt.lgr.dapeak.table <- subset(da_peaks.31.m, subset=p_val_adj < .05)
wnt.lgr.dapeak.table.order <- wnt.lgr.dapeak.table[order(wnt.lgr.dapeak.table$avg_log2FC),]
wnt.lgr.dapeak.table.order$row.num <- seq.int(nrow(wnt.lgr.dapeak.table.order))
#--and Enriched motifs (and their associated TFs) of DAPs. --> The list of TFs are specifically activated on each of conditions.
#---Preparation Enriched motifs of 3's DAGs
pfm <- getMatrixSet( x = JASPAR2020, opts = list(collection = "CORE", tax_group='vertebrates', all_versions = FALSE))
ready_to_visualize <- AddMotifs( object = ready_to_visualize, genome = linkage_genome, pfm = pfm)
#---Wnt7 VS Lgr5
wnt.lgr.enriched.motifs <- FindMotifs(object = ready_to_visualize, features = rownames(da_peaks.31.m[da_peaks.31.m$p_val_adj < 0.005, ]))
wnt.lgr.enriched.motifs[wnt.lgr.enriched.motifs$p.adjust < 0.05,]

#Enriched motifs for Wnt7b population and Lgr5 stem cell population in mutant
options(repr.plot.width=6, repr.plot.height=6)
wnt.lgr.enriched.motifs$logP <- -log10(wnt.lgr.enriched.motifs$p.adjust)
head(wnt.lgr.enriched.motifs[wnt.lgr.enriched.motifs$p.adjust < 0.05,],15)
data <- wnt.lgr.enriched.motifs[1:15,]
y.orders <- rev(data$motif.name)
p7 <- ggplot(data, aes(x=logP, y=motif.name, fill=percent.observed, names.arg=logP)) + 
    geom_bar(stat="identity", width=0.2) + scale_fill_gradientn(colours = c("royalblue", "rosybrown2", "red"), limits=c(0,60)) + 
    scale_y_discrete(limits = y.orders) + xlim(0,30) +
    ggtitle("Wnt7b+ vs Lgr5+ \n (RZK)") + xlab("-logP") + ylab("TF motif") + theme_linedraw()

p7 <- p17 + theme(
       plot.title = element_text(hjust = 0.5, color = "black", size = 15, face="bold"),
       axis.title.x = element_text(color = "black", size = 12),
       axis.title.y = element_text(color = "black", size = 12)
    )

ggsave(file="/hdds/hdd1/3_RZK_from_IMBA/youngchul/pdf/fig_3_I_Enriched_motifs_and_TFs_Wnt_population_vs_Lgr5_on_RZK.pdf", width=6, height=6, dpi=300, plot=p7)

dev.off()

### Wnt7 vs Lgr5 in wlid type
closest_genes <- ClosestFeature(ready_to_visualize, regions = rownames(da_peaks.31.w))
rownames(closest_genes) <- closest_genes$query_region
da_peaks.31.w$gene <- closest_genes$gene_name
da_peaks.31.w$gene_biotype <- closest_genes$gene_biotype
da_peaks.31.w$closest_region <- closest_genes$closest_region
da_peaks.31.w$distance <- closest_genes$distance
wnt.lgr.dapeak.table <- subset(da_peaks.31.w, subset=p_val_adj < .05)
wnt.lgr.dapeak.table.order <- wnt.lgr.dapeak.table[order(wnt.lgr.dapeak.table$avg_log2FC),]
wnt.lgr.dapeak.table.order$row.num <- seq.int(nrow(wnt.lgr.dapeak.table.order))
#--and Enriched motifs (and their associated TFs) of DAPs. --> The list of TFs are specifically activated on each of conditions.
#---Preparation Enriched motifs of 3's DAGs
pfm <- getMatrixSet( x = JASPAR2020, opts = list(collection = "CORE", tax_group='vertebrates', all_versions = FALSE))
ready_to_visualize <- AddMotifs( object = ready_to_visualize, genome = linkage_genome, pfm = pfm)
#---Wnt7 VS Lgr5
wnt.lgr.enriched.motifs.mutant <- FindMotifs(object = ready_to_visualize, features = rownames(da_peaks.31.w[da_peaks.31.w$p_val_adj < 0.005, ]))
wnt.lgr.enriched.motifs.mutant[wnt.lgr.enriched.motifs.mutant$p.adjust < 0.05,]

#Enriched motifs for Wnt7b population and Lgr5 stem cell population in wild type
options(repr.plot.width=6, repr.plot.height=6)
wnt.lgr.enriched.motifs.mutant$logP <- -log10(wnt.lgr.enriched.motifs.mutant$p.adjust)
head(wnt.lgr.enriched.motifs.mutant[wnt.lgr.enriched.motifs.mutant$p.adjust < 0.05,],15)
data <- wnt.lgr.enriched.motifs.mutant[1:15,]
y.orders <- rev(data$motif.name)
p8 <- ggplot(data, aes(x=logP, y=motif.name, fill=percent.observed, names.arg=logP)) + 
    geom_bar(stat="identity", width=0.2) + scale_fill_gradientn(colours = c("royalblue", "rosybrown2", "red"), limits=c(0,60)) + 
    scale_y_discrete(limits = y.orders) + xlim(0,30) +
    ggtitle("Wnt7b+ vs Lgr5+ \n (RZ)") + xlab("-logP") + ylab("TF motif") + theme_linedraw()

p8 <- p8 + theme(
       plot.title = element_text(hjust = 0.5, color = "black", size = 15, face="bold"),
       axis.title.x = element_text(color = "black", size = 12),
       axis.title.y = element_text(color = "black", size = 12)
    )

ggsave(file="/hdds/hdd1/3_RZK_from_IMBA/youngchul/pdf/fig_3_H_Enriched_motifs_and_TFs_Wnt_population_vs_Lgr5_on_RZ.pdf", width=6, height=6, dpi=300, plot=p8)
dev.off()

#### fig S3A
options(repr.plot.width=10, repr.plot.height=10)

plot_featureplot_v3(RZ.RZK, "Muc5ac",0,3)
plot_featureplot_v3(RZ.RZK, "Gkn2",0,3)
plot_featureplot_v3(RZ.RZK, "Tff1",0.5,6)
plot_featureplot_v3(RZ.RZK, "Gkn1",0.5,5)
plot_featureplot_v3(RZ.RZK, "Smc2",0.5,3)
plot_featureplot_v3(RZ.RZK, "Top2a",0.5,4)
plot_featureplot_v3(RZ.RZK, "Hmgb2",1,4)
plot_featureplot_v3(RZ.RZK, "Foxm1",0,2)
plot_featureplot_v3(RZ.RZK, "Mki67",0.3,3)
plot_featureplot_v3(RZ.RZK, "Muc6",0.3,4.5)
plot_featureplot_v3(RZ.RZK, "Cftr",0.5,3.5)
plot_featureplot_v3(RZ.RZK, "Cd44",2,4.2)
plot_featureplot_v3(RZ.RZK, "Glipr1",1.5,4.5)
plot_featureplot_v3(RZ.RZK, "Lgr5",0.5,3.5)

#### fig S3B
DefaultAssay(RZ.RZK) <- "RNA"

genes <- c("Wnt1", "Wnt2", "Wnt2b", "Wnt3", "Wnt3a", "Wnt4", "Wnt5a", "Wnt5b", "Wnt6",
           "Wnt7a", "Wnt7b", "Wnt8a", "Wnt8b", "Wnt9a", "Wnt9b", "Wnt10a", "Wnt10b", "Wnt11", "Wnt16")

VlnPlot(
  object = RZ.RZK,
  features = genes,
  group.by = "kras_genotype" 
)
