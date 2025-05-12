library(magrittr)
library(ggplot2)
library(pheatmap)
library(biomaRt)
library(dplyr)
library(readr)
library(tibble)
library(edgeR)
library(GSEABase)
library(GSVAdata)
library(tximport)

tx2gene_data = as.data.frame(read.csv("./0_result/star_salmon/salmon_tx2gene.tsv",sep="\t"))
files<- list.files(path = "./0_result/star_salmon", pattern = "quant.sf", full.names = TRUE, 
                   recursive = TRUE)
names(files)<- stringr::str_split(files, pattern = "/", simplify = TRUE)[,12] %>% # 12 should be changed based on path
  stringr::str_replace("_quant", "")
final_files = c(final_files,files)

txi.salmon <- tximport(final_files, type = "salmon", tx2gene = tx2gene_data)

cts <- txi.salmon$counts
normMat <- txi.salmon$length

# Obtaining per-observation scaling factors for length, adjusted to avoid
# changing the magnitude of the counts.
normMat <- normMat/exp(rowMeans(log(normMat)))
normCts <- cts/normMat

# Computing effective library sizes from scaled counts, to account for
# composition biases between samples.

eff.lib <- calcNormFactors(normCts) * colSums(normCts)

# Combining effective library sizes with the length factors, and calculating
# offsets for a log-link GLM.
normMat <- sweep(normMat, 2, eff.lib, "*")
normMat <- log(normMat)

# Creating a DGEList object for use in edgeR.
y <- DGEList(cts)
y <- scaleOffset(y, normMat)

cpms <- edgeR::cpm(y, offset = y$offset, log = FALSE)

write.csv(cpms, "00250508_tmm_normalized_counts_all_sample_final.csv")