library(tidyverse)
library(readxl)
library(ggnewscale)
library(dplyr)
library(readr)
library(dendextend)

expr_df   <- read_csv("/vol/RZRZK/data/OmicsExpressionProteinCodingGenesTPMLogp1BatchCorrected.csv")

colnames(expr_df) <- sapply(colnames(expr_df), function(col) {
  if (grepl("\\(", col)) {
    strsplit(col, " ")[[1]][1]  # example: "EGFR (some_info)" → "EGFR"
  } else {
    col
  }
})

expr_df <- expr_df %>%
  rename(depmap_id = '...1') %>%
  column_to_rownames("depmap_id")

expr_df <- tibble::rownames_to_column(expr_df, "depmap_id")

library(tidyverse)

# Load correlation data
correlation_data <- read_csv("/vol/RZRZK/output/correlation/Tissue_WNT7B_corr_hallmark_spearman.csv") %>%
  rename(Gene = `...1`) %>%
  column_to_rownames("Gene") %>%
  tibble::rownames_to_column("Gene")
metadata <- read_csv("/vol/RZRZK/data/CCLE_metadata.csv")

# Extract list of genes with non-NA and positive values in the 'Stomach' column
genes_of_interest <- correlation_data %>%
  filter(!is.na(Stomach), Stomach > 0) %>%
  pull(Gene)

# Define a manually selected gene list and filter by the above list 33 intersect genes
genes_of_interest2 <- c('EREG','CPE','EMP1','CCND2','IGFBP3','TSPAN1','KIF5C','ITGA2','KCNN4','ANXA10',
                        'PTPRR','ANGPTL4','PTGS2','GPRC5B','TMEM158','CCL20','ETV4','PLAU','WNT7A',
                        'CIDEA','INHBA','ADAMDEC1','RGS16','PPBP','LIF','TRIB2','ETV1','MALL',
                        'CFB','GNG11','ALDH1A2','CSF2','PRDM1')
genes_of_interest2 <- intersect(genes_of_interest2, genes_of_interest)

# Convert expression data to long format
expression_long <- expr_df %>%
  pivot_longer(-depmap_id, names_to = "Gene", values_to = "Expression")

# Add tissue information based on depmap_id
expression_long <- expression_long %>%
  left_join(metadata, by = c("depmap_id" = "depmap_id"))

# Initialize result storage list
result_list <- list()

# Iterate over each tissue
for (tissue in unique(metadata$Tissue)) {
  # Get sample IDs belonging to the current tissue
  samples_in_tissue <- metadata %>%
    filter(Tissue == tissue) %>%
    pull(depmap_id)
  
  # Filter expression data for the current tissue samples
  tissue_expression <- expression_long %>%
    filter(depmap_id %in% samples_in_tissue)
  
  # Iterate over each gene of interest
  for (gene in genes_of_interest2) {
    # Divide samples into quartiles based on expression of the current gene
    gene_expression <- tissue_expression %>%
      filter(Gene == gene) %>%
      select(depmap_id, Expression) %>%
      mutate(quartile = ntile(Expression, 4))
    
    # Define high (4th quartile) and low (1st quartile) groups
    high_samples <- gene_expression %>% filter(quartile == 4) %>% pull(depmap_id)
    low_samples  <- gene_expression %>% filter(quartile == 1) %>% pull(depmap_id)
    
    # Filter WNT7B expression values
    wnt7b_expression <- tissue_expression %>%
      filter(Gene == "WNT7B") %>%
      select(depmap_id, Expression)
    
    # Calculate average WNT7B expression in high group
    high_mean_wnt7b <- wnt7b_expression %>%
      filter(depmap_id %in% high_samples) %>%
      summarise(mean_WNT7B = mean(Expression, na.rm = TRUE)) %>%
      pull(mean_WNT7B)
    
    # Calculate average WNT7B expression in low group
    low_mean_wnt7b <- wnt7b_expression %>%
      filter(depmap_id %in% low_samples) %>%
      summarise(mean_WNT7B = mean(Expression, na.rm = TRUE)) %>%
      pull(mean_WNT7B)
    
    # Retrieve correlation value between WNT7B and the current gene in this tissue (if available)
    if (tissue %in% colnames(correlation_data)) {
      correlation_value <- correlation_data %>%
        filter(Gene == gene) %>%
        select(all_of(tissue)) %>%
        pull()
    } else {
      correlation_value <- NA
    }
    
    # Extract WNT7B expression values for high and low groups (for Wilcoxon test)
    high_wnt_values <- wnt7b_expression %>% filter(depmap_id %in% high_samples) %>% pull(Expression)
    low_wnt_values  <- wnt7b_expression %>% filter(depmap_id %in% low_samples) %>% pull(Expression)
    
    # Perform one-tailed Wilcoxon test (High > Low)
    if(length(high_wnt_values) >= 2 & length(low_wnt_values) >= 2) {
      ttest_result <- wilcox.test(high_wnt_values, low_wnt_values, alternative = "greater")
      ttest_p <- ttest_result$p.value
    } else {
      ttest_p <- NA
    }
    
    # Save results in a standardized format (same t-test p-value for High and Low)
    result_list[[paste(tissue, gene, "high", sep = "_")]] <- data.frame(
      Tissue = tissue,
      Gene = gene,
      Group = "High",
      WNT7B = high_mean_wnt7b,
      Correlation = correlation_value,
      Ttest_P = ttest_p
    )
    
    result_list[[paste(tissue, gene, "low", sep = "_")]] <- data.frame(
      Tissue = tissue,
      Gene = gene,
      Group = "Low",
      WNT7B = low_mean_wnt7b,
      Correlation = correlation_value,
      Ttest_P = ttest_p
    )
  }
}

# Combine all results into a single data frame
final_result <- bind_rows(result_list)

# Print result
print(final_result)

# If the Tissue_clean column already exists, use it; otherwise, extract the part before the parenthesis from Tissue
final_result <- final_result %>%
  mutate(Tissue_clean = sub(" \\(.*", "", Tissue))

# For each gene, count how many unique Tissue_clean values have Correlation > 0
gene_positive_counts <- final_result %>%
  group_by(Gene, Tissue_clean) %>%
  summarise(is_positive = any(Correlation > 0, na.rm = TRUE)) %>% 
  ungroup() %>%
  group_by(Gene) %>%
  summarise(num_positive = sum(is_positive)) %>%
  ungroup()

# Keep only genes with positive correlation in more than 5 tissues
genes_to_keep <- gene_positive_counts %>%
  filter(num_positive > 5) %>%
  pull(Gene)

# Filter final result to include only the selected genes
final_result_filtered <- final_result %>%
  filter(Gene %in% genes_to_keep)

# Create a new variable combining Gene and Group (High/Low) for Y-axis sorting
final_result_filtered <- final_result_filtered %>%
  mutate(Gene_Group = paste(Gene, Group, sep = "_"))

# Count number of samples per Tissue
tissue_sample_count <- metadata %>%
  group_by(Tissue) %>%
  summarise(sample_count = n()) %>%
  ungroup()

# Append sample count to Tissue name and merge with the final data
final_result_filtered <- final_result_filtered %>%
  left_join(tissue_sample_count, by = "Tissue")

# Remove tissues with 20 or fewer samples
final_result_filtered <- final_result_filtered %>%
  filter(sample_count > 20)

# Append sample count to Tissue names for display
final_result_filtered <- final_result_filtered %>%
  mutate(Tissue = paste0(Tissue, " (", sample_count, ")"))

# Check final result
print(final_result_filtered)


############## fig 7B upper heatmap
# 1. Extract unique Gene, Tissue, and Correlation values from final_result_filtered
heatmap_df <- final_result_filtered %>%
  distinct(Tissue_clean, Gene, Correlation)

# 2. Convert to wide format: rows = Gene, columns = Tissue
heatmap_wide <- heatmap_df %>%
  pivot_wider(names_from = Tissue_clean, values_from = Correlation)

# 3. Convert to matrix: set rownames to Gene
heatmap_mat <- as.data.frame(heatmap_wide)
rownames(heatmap_mat) <- heatmap_mat$Gene
heatmap_mat <- heatmap_mat %>% select(-Gene)
heatmap_mat <- as.matrix(heatmap_mat)

# Create binary matrix (1 if correlation > 0, 0 otherwise, NA = 0)
binary_mat <- ifelse(is.na(heatmap_mat), 0, ifelse(heatmap_mat > 0, 1, 0))

# Set desired gene and tissue order for dendrogram rotation
new_order_gene_rotate <- rev(c("ANGPTL4","EMP1","KCNN4","PTPRR","CSF2","EREG","LIF","INHBA","PTGS2","IGFBP3","TSPAN1","ITGA2"))
target_tissues_new_order <- c('Stomach','Lung','Bone','Peripheral Nervous System','Ovary/Fallopian Tube','Skin',
                              'Soft Tissue','Kidney','Lymphoid','Pancreas','Bowel','Liver','Uterus','Breast',
                              'Bladder/Urinary Tract','Biliary Tract','Head and Neck','Myeloid','Esophagus',
                              'CNS/Brain',"Pleura",'Thyroid')

# Column (Tissue) clustering: compute distance on transposed binary matrix
col_dist <- dist(t(binary_mat), method = "binary")
col_hclust <- hclust(col_dist, method = "average")
# Convert to dendrogram for dendextend
col_dend <- as.dendrogram(col_hclust)

# Reorder dendrogram using rotate()
col_dend <- rotate(col_dend, order = target_tissues_new_order)

col_hclust2 <- as.hclust(col_dend)

# Row (gene) clustering
row_dist <- dist(binary_mat, method = "binary")
row_hclust <- hclust(row_dist, method = "average")
# Convert to dendrogram for dendextend
row_dend <- as.dendrogram(row_hclust)

# Reorder dendrogram using rotate()
row_dend <- rotate(row_dend, order = new_order_gene_rotate)

row_hclust2 <- as.hclust(row_dend)

# 6. Count sample numbers for each Tissue_clean using metadata
tissue_counts <- metadata %>%
  mutate(Tissue_clean = sub(" \\(.*", "", Tissue)) %>%
  count(Tissue_clean, name = "sample_count")

# 7. Append sample count to heatmap column names
new_colnames <- sapply(colnames(heatmap_mat), function(x) {
  cnt <- tissue_counts %>% filter(Tissue_clean == x) %>% pull(sample_count)
  if(length(cnt) == 0) x else paste0(x, " (", cnt, ")")
})
colnames(heatmap_mat) <- new_colnames

# Plot heatmap with clustering
p1 <- pheatmap(heatmap_mat,
               cluster_rows = row_hclust2,
               cluster_cols = col_hclust2,
               color = colorRampPalette(c("#195696","#2F79B5","#4E9AC6","#86BDDA","#B6D7E8","#DBEAF2",
                                          "#F7F7F7","#FBE3D4","#F9C3A9","#EF9B7A","#DA6954","#C13639", "#9C1127"))(50),
               breaks = seq(-0.6, 0.6, length.out = 51),
               angle_col = 45,
               main = "Gene-Tissue Correlation Heatmap"
)
print(p1)

# Save the plot as a PDF
ggsave(plot = p1, width = 10, height = 10, filename = "/vol/RZRZK/output/pdf/tissue_wnt_corr_filter_dendro_deg_250408.pdf")


############## fig 7B middle heatmap
# 1. Prepare p-value data: extract Tissue_clean and compute -log10(p-value), set p > 0.05 to NA
pval_df <- final_result_filtered %>%
  mutate(Tissue_clean = sub(" \\(.*", "", Tissue)) %>%
  distinct(Tissue_clean, Gene, Ttest_P) %>% 
  mutate(logP = -log10(Ttest_P))

# 2. Convert to wide format: rows = Gene, columns = Tissue_clean, values = logP
pval_wide <- pval_df %>%
  select(Tissue_clean, Gene, logP) %>%
  pivot_wider(names_from = Tissue_clean, values_from = logP)

# 3. Convert to matrix and set rownames to Gene
pval_mat <- as.data.frame(pval_wide)
rownames(pval_mat) <- pval_mat$Gene
pval_mat <- pval_mat %>% select(-Gene)
pval_mat <- as.matrix(pval_mat)

# 4. Manual ordering: reorder genes and tissues according to user-defined order
manual_gene_order <- rev(c("ANGPTL4","EMP1","KCNN4","PTPRR","CSF2","EREG","LIF","INHBA","PTGS2","IGFBP3","TSPAN1","ITGA2"))
manual_tissue_order <- c('Stomach','Lung','Bone','Peripheral Nervous System','Ovary/Fallopian Tube','Skin',
                         'Soft Tissue','Kidney','Lymphoid','Pancreas','Bowel','Liver','Uterus','Breast',
                         'Bladder/Urinary Tract','Biliary Tract','Head and Neck','Myeloid','Esophagus',
                         'CNS/Brain',"Pleura",'Thyroid')

# Subset and reorder matrix based on the defined gene and tissue order
pval_mat <- pval_mat[intersect(manual_gene_order, rownames(pval_mat)),
                     intersect(manual_tissue_order, colnames(pval_mat)), drop = FALSE]

options(repr.plot.width = 10, repr.plot.height = 10)

# Count number of samples per Tissue_clean using metadata
tissue_counts <- metadata %>%
  mutate(Tissue_clean = sub(" \\(.*", "", Tissue)) %>%
  count(Tissue_clean, name = "sample_count")

# 7. Append sample count to column names of the pval matrix
new_colnames <- sapply(colnames(pval_mat), function(x) {
  cnt <- tissue_counts %>% filter(Tissue_clean == x) %>% pull(sample_count)
  if(length(cnt) == 0) x else paste0(x, " (", cnt, ")")
})
colnames(pval_mat) <- new_colnames

# 5. Draw heatmap using pheatmap: no clustering, preserve manual order, rotate x-axis labels
p_heatmap <- pheatmap(pval_mat,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("#FFFFFF", "#FFC6FF", "#DD7EFF", "#A020F0", "#5C1DB5"))(50),
         breaks = seq(0, 3, length.out = 51),
         main = "-log10(t-test p-value) Heatmap\n(p > 0.05 set to NA)",
         angle_col = 45,
         fontsize = 10,
         na_col = "grey90")

# Save the heatmap as PDF
ggsave(filename = "/vol/RZRZK/output/pdf/heatmap_ttest_purple_spearman_deg_250408.pdf",
       plot = p_heatmap, width = 10, height = 10)

############## fig 7B under heatmap

# 1. Extract only WNT7B expression (assumes expression_long and metadata are preloaded)
wnt7b_df <- expression_long %>% 
  filter(Gene == "WNT7B") %>%
  left_join(metadata, by = "depmap_id")  # Add tissue information

# 2. Compute average WNT7B expression for each tissue
wnt7b_mean <- wnt7b_df %>%
  group_by(Tissue.y) %>%
  summarise(WNT7B_mean = mean(Expression, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Tissue_clean = sub(" \\(.*", "", Tissue.y))

# 3. Convert to wide format: row = "WNT7B", columns = Tissue_clean, values = WNT7B_mean  
#    (only one gene, so single-row matrix)
heatmap_df <- wnt7b_mean %>%
  select(Tissue_clean, WNT7B_mean) %>%
  pivot_wider(names_from = Tissue_clean, values_from = WNT7B_mean)

# 4. Convert to matrix
heatmap_mat <- as.data.frame(heatmap_df)
rownames(heatmap_mat) <- "WNT7B"  # Only one row
heatmap_mat <- heatmap_mat %>% select(-1)  # Remove the first column
heatmap_mat <- as.matrix(heatmap_mat)

# 5. Manual ordering:
#    manual_gene_order contains only "WNT7B", manual_tissue_order defines custom tissue order
manual_gene_order <- "WNT7B"
manual_tissue_order <- c('Stomach','Lung','Bone','Peripheral Nervous System','Ovary/Fallopian Tube','Skin',
                         'Soft Tissue','Kidney','Lymphoid','Pancreas','Bowel','Liver','Uterus','Breast',
                         'Bladder/Urinary Tract','Biliary Tract','Head and Neck','Myeloid','Esophagus',
                         'CNS/Brain',"Pleura",'Thyroid')

# 6. Reorder columns using manual_tissue_order (keep only tissues present in data)
available_tissues <- intersect(manual_tissue_order, colnames(heatmap_mat))
heatmap_mat <- heatmap_mat[, available_tissues, drop = FALSE]
options(repr.plot.width = 10, repr.plot.height = 2)

# 7. Plot heatmap using pheatmap (no clustering, custom order, rotate x-axis labels)
p <- pheatmap(heatmap_mat,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("white", "black"))(50),
         breaks = seq(0, 4, length.out = 51),
         main = "Average WNT7B Expression (per tissue)",
         angle_col = 45)
p

# 8. Save result
ggsave(plot = p, width = 10, height = 2, filename = "/vol/RZRZK/output/pdf/tissue_wnt7b_expression_spearman_deg_dendro_250403.pdf")



