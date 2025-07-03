#20250512 : This code is main code to make clusters for RZK project
#In v2 normalization == lognormalization (for DE analysis)
#Integration (SCT) ==> normalize --> vaiable features --> integration
#plot_linkage_cl_genotype is changed with gene info as argument
###Libararies
#Conda env
library(future)
library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(EnsDb.Mmusculus.v79)
library(RColorBrewer)
library(clustree)
library(hash)
library(BSgenome.Mmusculus.UCSC.mm10)
library(dittoSeq)
library(motifmatchr)
library(JASPAR2020)
library(TFBSTools)
library(chromVAR)
library(msigdbr)
library(gprofiler2)

#BiocManager (https://cran.r-project.org/web/packages/BiocManager/vignettes/BiocManager.html)
#R install.packages
#install.packages("BiocManager")
#conda install -c conda-forge r-biocmanager
#conda install -c conda-forge r-hdf5r
#BiocManager (install guide below)
#BiocManager::install("XVector")
#BiocManager::install("GO.db")
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38", force=TRUE)
#BiocManager::install("BSgenome.Mmusculus.UCSC.mm10", force=TRUE)
#BiocManager::install("dittoSeq", force=TRUE)
#BiocManager::install("motifmatchr", force=TRUE)
#BiocManager::install("JASPAR2020", force=TRUE)
#BiocManager::install("TFBSTools", force=TRUE)
#BiocManager::install("chromVAR", force=TRUE)
#BiocManager::install("biovizBase", force=TRUE)


###Setting Multiprocess
plan("multicore", workers = 20)
options(future.globals.maxSize = 80 * 1024 ^ 3) # for 80 Gb RAM
future.seed=TRUE
set.seed(1234)
###Files
#Wild_type
wild.h5.file='/data/snMultiomics/Mouse_RZ_WEN/filtered_feature_bc_matrix.h5'
wild.fragment.file='/data/snMultiomics/Mouse_RZ_WEN/atac_fragments.tsv.gz'
wild.fragment.tbi.file='/data/snMultiomics/Mouse_RZ_WEN/atac_fragments.tsv.gz.tbi'
wild.meta.file='/data/snMultiomics/Mouse_RZ_WEN/per_barcode_metrics.csv'
#Mutant
mutant.h5.file='/data/snMultiomics/Mouse_RZK_WEN/filtered_feature_bc_matrix.h5'
mutant.fragment.file='/data/snMultiomics/Mouse_RZK_WEN/atac_fragments.tsv.gz'
mutant.fragment.tbi.file='/data/snMultiomics/Mouse_RZK_WEN/atac_fragments.tsv.gz.tbi'
mutant.meta.file='/data/snMultiomics/Mouse_RZK_WEN/per_barcode_metrics.csv'
###Arguments
species <- "mouse"
if (species=='mouse'){
    linkage_genome = BSgenome.Mmusculus.UCSC.mm10   
} else {
    linkage_genome = BSgenome.Hsapiens.UCSC.hg38   
}
###Internal functions
set_annotation <- function(species){
    if (species == "mouse") {
        annotation <- GetGRangesFromEnsDb(ensdb=EnsDb.Mmusculus.v79)
        seqlevelsStyle(annotation) <- "UCSC"
        genome(annotation) <- "mm10"
    } else if (species == "human") {
        annotation <- GetGRangesFromEnsDb(ensdb=EnsDb.Hsapiens.v86)
        seqlevelsStyle(annotation) <- "UCSC"
        genome(annotation) <- "hg38"
    }
    return (annotation)    
}
###2. MAKE BASE OBJECT FILE for MULTIOME
generate_object <- function(h5_file, frag_file, meta_file, species){
    #Metadata
    metadata <- read.csv(file = meta_file, header = TRUE, row.names = 1)
    #RNA part
    input.10x <- Read10X_h5(h5_file)
    rna.counts <- input.10x$'Gene Expression'
    obj <- CreateSeuratObject(counts = rna.counts, assay = "RNA", meta.data = metadata)
    if (species=='human'){
        obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
    } else if (species=='mouse') {
        obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
    }    
    #Setting ATAC-seq with only standard chromosomes
    atac.counts <- input.10x$Peaks
    grange.counts <- StringToGRanges(rownames(atac.counts), sep = c(":", "-"))
    grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
    atac.counts <- atac.counts[as.vector(grange.use), ]
    #Add Chromatin Assay, with pre-set annotation
    obj[["ATAC"]] <- CreateChromatinAssay(
        counts = atac.counts,
        sep = c(":", "-"),
        genome = genome(annotations),
        fragments = frag_file,
        min.cells=10,
        annotation = annotations
    )
    return (obj)
}
add_genotype <- function(obj, geno){
    obj@meta.data$kras_genotype <- geno
    return (obj)
}
###3. Quality control based on ATAC-seq and RNA-seq data
qc_function <- function(obj){
    DefaultAssay(obj) <- "ATAC"
    obj <- NucleosomeSignal(obj)
    obj <- TSSEnrichment(obj)
    VlnPlot(
    object = obj,
    features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
    ncol = 4,
    pt.size = 0
    )
    return (obj)
}
subsetting_objects_QC <- function(obj) {
    obj <- qc_function(obj)
    obj <- subset(
        x = obj,
        subset = nCount_ATAC < 100000 &
            nCount_RNA < 25000 &
            nCount_ATAC > 1000 &
            nCount_RNA > 1000 &
            nucleosome_signal < 2 &
            TSS.enrichment > 1 &
            percent.mt < 15
    )   
    return (obj)
}
###Peak Calling by MACS2 tool_REQUIRED!!
call_macs2_counts <- function(obj, species){
    # call peaks using MACS2
    peaks <- CallPeaks(obj,macs2.path = "/opt/conda/envs/RZK_project/bin/macs2",verbose = TRUE)
    # remove peaks on nonstandard chromosomes and in genomic blacklist regions
    if (species=='mouse'){
        blacklist = blacklist_mm10
    } else if (species=='human'){
        blacklist = blacklist_hg38_unified
    } else {print('check_species')}  
    peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")    
    peaks <- subsetByOverlaps(x = peaks, ranges = blacklist, invert = TRUE)

    # quantify counts in each peak
    macs2_counts <- FeatureMatrix(
        fragments = Fragments(obj),
        features = peaks,
        cells = colnames(obj)
    )
    return (macs2_counts)
}
add_peak_to_object <- function(obj, species, frag_file){
    obj[["macs2"]] <- CreateChromatinAssay(
        counts = call_macs2_counts(obj, species),
        fragments = frag_file,
        annotation = annotations
    )
    return (obj)
}
### 4. Dimensional reduction, nomralization, and clustering
##4-1. Gene expression data processing 
process_expression <- function(obj){
    obj <- SCTransform(obj) #Normalize, Scale data, FVF are included and normalized data is stored in [["SCT"]]@scale.data
    DefaultAssay(obj) <- "RNA"
    obj <- NormalizeData(obj)
    obj <- FindVariableFeatures(obj, selection.method="vst", nfeatures=2000)    
    return (obj)
}
##4-2. DNA accessibility data processing
process_accessibility <- function(obj){
    DefaultAssay(obj) <- "ATAC" 
    obj <- FindTopFeatures(obj, min.cutoff = 'q0')      
    obj <- RunTFIDF(obj)    
    obj <- RunSVD(obj)  
    DefaultAssay(obj) <- "macs2" 
    obj <- FindTopFeatures(obj, min.cutoff = 'q0')    
    obj <- RunTFIDF(obj)    
    obj <- RunSVD(obj)    
    return (obj)
}
##4-1. Integration bsaed on scRNA-seq results
#https://github.com/satijalab/seurat/issues/1836
integrate_based_on_rna <- function(obj1, obj2) {
    DefaultAssay(obj1) <- "SCT"
    DefaultAssay(obj2) <- "SCT"
    obj.list <- list(obj1, obj2)
    features <- SelectIntegrationFeatures(obj.list)
    obj.list <- PrepSCTIntegration(obj.list, anchor.features = features)
    obj.anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features, normalization.method = "SCT")
    obj.combined <- IntegrateData(anchorset = obj.anchors, normalization.method = "SCT")
    return (obj.combined)
}
###5. UMAP visualization
generate_umap_integrated <- function(obj){    
    DefaultAssay(obj) <- "integrated"    
    obj <- ScaleData(obj, verbose = FALSE)
    obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
    obj <- RunUMAP(obj, reduction = "pca", dims = 1:30)    
    #obj <- RunTSNE(obj, reduction = "pca",dims = 1:19)
    obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30)
    obj <- FindClusters(obj, resolution = 0.3, algorithm=1)    
    return (obj)
}
visualize_umap_integrated <- function(obj){
    obj$kras_genotype <- factor(obj$kras_genotype, levels = c("Wild_type", "Mutant"))
    cols <- c('#8936EF','#F2CA19','#FF00BD','#E11845','#0057E9','#87E911')
    p1 <- DimPlot(obj, reduction = "umap", group.by = "kras_genotype")
    p2 <- DimPlot(obj, reduction = "umap", label = TRUE, repel = TRUE, cols=cols)
    p1 + p2
}

#############################========
#+++++++++++++++++++++++++++++++++++#############################========
#+++++++++++++++++++++++++++++++++++#############################========
#+++++++++++++++++++++++++++++++++++#############################========
#+++++++++++++++++++++++++++++++++++#############################========
#+++++++++++++++++++++++++++++++++++
##Calling and generating data for figures
#Genome setting
annotations <- set_annotation(species)
#Generate Seurat object
wild <- generate_object(wild.h5.file, wild.fragment.file, wild.meta.file, 'mouse')
mutant <- generate_object(mutant.h5.file, mutant.fragment.file, mutant.meta.file, 'mouse')
#Add metadata
wild <- add_genotype(wild,'Wild_type')
mutant <- add_genotype(mutant,'Mutant')
#VlnPlot(wild, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#VlnPlot(mutant, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#QC and subsetting
wild <- subsetting_objects_QC(wild)
mutant <- subsetting_objects_QC(mutant)

Fragments(wild)
Fragments(mutant)
#MACS2 Peak calling
wild <- add_peak_to_object(wild, species, wild.fragment.file)
mutant <- add_peak_to_object(mutant, species, mutant.fragment.file)
#RNA data processing
wild <- process_expression(wild)
mutant <- process_expression(mutant)
#ATAC-seq data processing
wild <- process_accessibility(wild)
mutant <- process_accessibility(mutant)
#Integration objects based on RNA assay
wild_mutant.combined <- integrate_based_on_rna(wild, mutant)
#Clustering
ready_to_visualize <- generate_umap_integrated(wild_mutant.combined)
ready_to_visualize@meta.data$cl_genotype <- paste0(ready_to_visualize@meta.data$integrated_snn_res.0.3, "_", ready_to_visualize@meta.data$kras_genotype)

##Analyzing and Visualizing datasets from DE, DA analyses
#0. Basic UMAP plots
#ready_to_visualize$kras_genotype <- factor(ready_to_visualize$kras_genotype, levels = c("Wild_type", "Mutant"))
#0-A. Kras_genotype and cluster
DefaultAssay(ready_to_visualize) <- "RNA"
cols <- c('#8936EF','#F2CA19','#FF00BD','#E11845','#0057E9','#87E911')
##Cell type Identification of RZ_vs_RZK
#Clu 0 – Neck cells (Muc6, Cftr)
#Clu 1 – Lgr5+ (Lgr5)
#Clu 2 – SPEM (Glipr1) 
#Clu 3 – Wnt7b+ 
#Clu 4 – Proliferating (Ki67, Stmn1)
#Clu 5 – Pit (Tff1, Gkn1)
mylevel <- c("Lgr5+", "SPEM","Neck","Proliferating","Pit", "Wnt7+")
Idents(ready_to_visualize) <- ready_to_visualize@meta.data$integrated_snn_res.0.3
Idents(ready_to_visualize) <- factor(Idents(ready_to_visualize), levels= mylevel)
#my_cols <- c('Neck'='#8936EF','Lgr5+'='#F2CA19','SPEM'='#FF00BD','Wnt7+'='#E11845','Proliferating'='#0057E9','Pit'='#87E911')
my_cols <- c('0'='#8936EF','1'='#F2CA19','2'='#FF00BD','3'='#E11845','4'='#0057E9','5'='#87E911')
p_0a <- visualize_umap_integrated(ready_to_visualize)
p_0a1 <- visualize_umap_integrated(ready_to_visualize)

#### DEG for fig 7, S7
DEG <- FindMarkers(ready_to_visualize, group.by = "kras_genotype", ident.1 = "Mutant", test.use = "MAST")
write.csv(DEG,"/data/snMultiomics/RZRZK_DEG_scMAST.csv")
DEG <- read.csv("/data/snMultiomics/RZRZK_DEG_scMAST.csv")

temp.deg<- DEG[DEG$p_val_adj < 0.01 & DEG$avg_log2FC > 1,]

hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)

kras_signaling <- hallmark[hallmark$gs_name == "HALLMARK_KRAS_SIGNALING_UP", ]
kras_signaling <- kras_signaling$gene_symbol

length(kras_signaling) # 220 genes in Hallmark KRAS_SIGNALING_UP 

temp_sc_2 <- gconvert(query = temp.deg$X, organism = "hsapiens", target = "ENSG")

length(temp_sc_2$name) # 512 genes in RZ vs RZK DEGs
intersect(temp_sc_2$name, kras_signaling) # 33 intersect gene list

save.image(file="/data/snMultiomics/RZ_RZK_object_v20250701.RData")
