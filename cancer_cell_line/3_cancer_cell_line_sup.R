library(tidyverse)
library(readxl)
library(ggnewscale)
library(dplyr)
library(readr)
library(cowplot)

# Load correlation data
correlation_data <- read_csv("/vol/RZRZK/output/correlation/Tissue_WNT7B_corr_hallmark_spearman.csv")
correlation_data <- correlation_data %>%
  rename(Gene = '...1') %>%
  column_to_rownames("Gene")
correlation_data <- tibble::rownames_to_column(correlation_data, "Gene")

metadata <- read_csv("/vol/RZRZK/data/CCLE_metadata.csv")

# Extract genes with non-NA and positive correlation values in the "Stomach" column
genes_of_interest <- correlation_data %>%
  filter(!is.na(Stomach), Stomach > 0) %>%
  pull(Gene)

genes_of_interest2 <- c('EREG','CPE','EMP1','CCND2','IGFBP3','TSPAN1','KIF5C','ITGA2','KCNN4','ANXA10','PTPRR','ANGPTL4','PTGS2','GPRC5B','TMEM158','CCL20','ETV4','PLAU','WNT7A','CIDEA','INHBA','ADAMDEC1','RGS16','PPBP','LIF','TRIB2','ETV1','MALL','CFB','GNG11','ALDH1A2','CSF2','PRDM1')

# Convert expression data to long format
expression_long <- expr_df %>%
  pivot_longer(-depmap_id, names_to = "Gene", values_to = "Expression")

# Add tissue information using depmap_id
expression_long <- expression_long %>%
  left_join(metadata, by = "depmap_id")

# Initialize results list
result_list <- list()

# Loop through each tissue
for (tissue in unique(metadata$Tissue)) {
  # Get sample IDs for the current tissue
  samples_in_tissue <- metadata %>%
    filter(Tissue == tissue) %>%
    pull(depmap_id)

  # Filter expression data for the selected tissue
  tissue_expression <- expression_long %>%
    filter(depmap_id %in% samples_in_tissue)

  # Loop through genes of interest
  for (gene in genes_of_interest2) {
    # Divide gene expression into quartiles
    gene_expression <- tissue_expression %>%
      filter(Gene == gene) %>%
      select(depmap_id, Expression) %>%
      mutate(quartile = ntile(Expression, 4))

    # Define high (Q4) and low (Q1) groups
    high_samples <- gene_expression %>% filter(quartile == 4) %>% pull(depmap_id)
    low_samples  <- gene_expression %>% filter(quartile == 1) %>% pull(depmap_id)

    # Extract WNT7B expression
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

    # Perform Wilcoxon test if both groups have ≥ 2 samples
    high_vals <- wnt7b_expression %>% filter(depmap_id %in% high_samples) %>% pull(Expression)
    low_vals  <- wnt7b_expression %>% filter(depmap_id %in% low_samples) %>% pull(Expression)
    if(length(high_vals) >= 2 & length(low_vals) >= 2){
      t_test_pvalue <- wilcox.test(high_vals, low_vals, alternative = "greater")$p.value
    } else {
      t_test_pvalue <- NA
    }

    # Retrieve correlation value if tissue exists in correlation_data
    if(tissue %in% colnames(correlation_data)){
      correlation_value <- correlation_data %>%
        filter(Gene == gene) %>%
        pull(!!sym(tissue))
      if(length(correlation_value) == 0) { correlation_value <- NA }
    } else {
      correlation_value <- NA
    }

    # Save both high and low group results (same p-value)
    result_list[[paste(tissue, gene, "high", sep = "_")]] <- data.frame(
      Tissue = tissue,
      Gene = gene,
      Group = "High",
      WNT7B = high_mean_wnt7b,
      Correlation = correlation_value,
      t_test_pvalue = t_test_pvalue
    )

    result_list[[paste(tissue, gene, "low", sep = "_")]] <- data.frame(
      Tissue = tissue,
      Gene = gene,
      Group = "Low",
      WNT7B = low_mean_wnt7b,
      Correlation = correlation_value,
      t_test_pvalue = t_test_pvalue
    )
  }
}

# Combine all results
final_result <- bind_rows(result_list)

# Create Gene_Group (e.g., "Gene_High", "Gene_Low") for Y-axis
final_result <- final_result %>%
  mutate(Gene_Group = paste(Gene, Group, sep = "_"))

# Add sample count per tissue
tissue_sample_count <- metadata %>%
  group_by(Tissue) %>%
  summarise(sample_count = n()) %>%
  ungroup()

# Append sample count to tissue names and merge
final_result <- final_result %>%
  left_join(tissue_sample_count, by = "Tissue") %>%
  filter(sample_count > 20) %>%
  mutate(Tissue = paste0(Tissue, " (", sample_count, ")"))

# Sort X-axis (Tissue): prioritize Stomach and sort by # of positive correlations
tissue_order <- final_result %>%
  group_by(Tissue) %>%
  summarise(positive_corr_count = sum(Correlation > 0, na.rm = TRUE)) %>%
  arrange(desc(positive_corr_count)) %>%
  pull(Tissue)
tissue_order <- c(tissue_order[grepl("^Stomach", tissue_order)], setdiff(tissue_order, tissue_order[grepl("^Stomach", tissue_order)]))

final_result <- final_result %>%
  mutate(Tissue = factor(Tissue, levels = tissue_order, ordered = TRUE))

# Sort Y-axis: Gene_Group (Low before High)
gene_order <- final_result %>%
  group_by(Gene) %>%
  summarise(na_count = sum(is.na(Correlation))) %>%
  arrange(desc(na_count)) %>%
  pull(Gene)

final_result <- final_result %>%
  mutate(Gene = factor(Gene, levels = gene_order, ordered = TRUE)) %>%
  arrange(Gene, factor(Group, levels = c("Low", "High"))) %>%
  mutate(Gene_Group = factor(Gene_Group, levels = unique(Gene_Group), ordered = TRUE))

# Optional CSV export
# write_csv(final_result, "/vol/RZRZK/output/division4/high_low_ttest_results.csv")

# Manual gene and tissue order for figure layout
gene_rotate <- c('CCND2','KIF5C','TMEM158','ETV4','ADAMDEC1','RGS16','TRIB2','ETV1','GNG11','CIDEA','ALDH1A2','CPE','GPRC5B','PPBP','CCL20','PLAU','PRDM1','ANXA10','CFB','WNT7A','MALL',"ANGPTL4","EMP1","KCNN4","PTPRR","CSF2","EREG","LIF","INHBA","PTGS2","IGFBP3","TSPAN1","ITGA2")
tissue_rotate <- c('Stomach','Lung','Bone','Peripheral Nervous System','Ovary/Fallopian Tube','Skin','Soft Tissue','Kidney','Lymphoid','Pancreas','Bowel','Liver','Uterus','Breast','Bladder/Urinary Tract','Biliary Tract','Head and Neck','Myeloid','Esophagus','CNS/Brain',"Pleura",'Thyroid')

# Extract Tissue_clean from "Tissue (sample_count)"
final_result <- final_result %>%
  mutate(Tissue_clean = sub(" \\(.*", "", Tissue))

# Reorder tissue based on manual list
tissue_ordered <- final_result %>% 
  distinct(Tissue, Tissue_clean) %>% 
  arrange(match(Tissue_clean, tissue_rotate)) %>% 
  pull(Tissue)

final_result <- final_result %>%
  mutate(Tissue = factor(Tissue, levels = tissue_ordered, ordered = TRUE))

# Reorder gene and gene_group
final_result <- final_result %>%
  mutate(Gene = factor(Gene, levels = gene_rotate, ordered = TRUE))

all_gene_groups <- c()
for(g in gene_rotate) {
  all_gene_groups <- c(all_gene_groups, paste0(g, "_Low"), paste0(g, "_High"))
}

final_result <- final_result %>%
  mutate(Gene_Group = paste0(Gene, "_", Group)) %>%
  mutate(Gene_Group = factor(Gene_Group, levels = all_gene_groups, ordered = TRUE))

final_result <- final_result[!is.na(final_result$Gene), ]

# ----------------------------
# Plot: Figure S7B - Dot heatmap
# ----------------------------
grid_data <- expand.grid(Tissue = levels(final_result$Tissue), Gene_Group = levels(final_result$Gene_Group))
grid_data <- grid_data %>% filter(row_number() %% 2 == 0)
y_grid_lines <- final_result %>% group_by(Gene) %>% summarise(y_pos = max(as.numeric(Gene_Group))) %>% pull(y_pos)
x_grid_lines <- seq(1.5, length(tissue_ordered) - 0.5, by = 1)
xmin <- 0.5
xmax <- length(tissue_ordered) + 0.5
ymin <- 0.5
ymax <- length(levels(final_result$Gene_Group)) + 0.5

options(repr.plot.width = 10, repr.plot.height = 14)

p <- ggplot(final_result, aes(x = Tissue, y = Gene_Group)) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA, color = "black", size = 1) +
  geom_tile(aes(fill = Correlation), color = NA) +
  geom_point(data = subset(final_result, t_test_pvalue <= 0.05),
             aes(size = WNT7B),
             color = "black", stroke = 0.5, shape = 21, fill = "white") +
  geom_hline(yintercept = y_grid_lines + 0.5, color = "black", size = 0.5) +
  geom_vline(xintercept = x_grid_lines, color = "black", size = 0.5) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-0.6, 0.6), oob = scales::squish, na.value = "grey90") +
  scale_size(range = c(0, 6)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    title = "Tissue-Gene High/Low WNT7B Expression & Correlation\n(with Wilcoxon test p-values)",
    x = "Tissue (Sample Count > 20)",
    y = "Kras Signature geneset (High/Low)",
    fill = "Correlation",
    size = "WNT7B \nMean Expression"
  )

ggsave(plot = p, width = 10, height = 14, filename = "/vol/RZRZK/output/pdf/tissue_wnt_corr_spearman_deg_250408.pdf")

# ----------------------------
# Plot: Figure S7C - Spearman rank scatterplots in Stomach
# ----------------------------
stomach_expr <- expression_long %>%
  mutate(Tissue_clean = sub(" \\(.*", "", Tissue)) %>%
  filter(Tissue_clean == "Stomach")

genes_to_analyze <- c("ANGPTL4","EMP1","PTPRR","EREG","CSF2","LIF","INHBA","KCNN4","PTGS2","TSPAN1","IGFBP3","ITGA2")

corr_results_list <- list()

for(g in genes_to_analyze) {
  df_gene <- stomach_expr %>% filter(Gene == g)
  df_wnt7b <- stomach_expr %>%
    filter(Gene == "WNT7B") %>%
    select(depmap_id, Expression) %>%
    rename(WNT7B_Expression = Expression)

  df_join <- left_join(df_gene, df_wnt7b, by = "depmap_id") %>%
    filter(!is.na(Expression), !is.na(WNT7B_Expression))

  if(nrow(df_join) < 5) next

  cor_test <- cor.test(df_join$WNT7B_Expression, df_join$Expression, method = "spearman")
  r_val <- cor_test$estimate
  p_val <- cor_test$p.value

  corr_results_list[[g]] <- data.frame(
    Gene = g,
    spearman_R = r_val,
    spearman_P = p_val
  )

  df_join <- df_join %>%
    mutate(
      WNT7B_Expression_Rank = rank(WNT7B_Expression),
      Gene_Expression_Rank = rank(Expression)
    )

  p_scatter <- ggplot(df_join, aes(x = WNT7B_Expression_Rank, y = Gene_Expression_Rank)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(
      title = paste(g, "vs WNT7B (Rank-based)"),
      subtitle = paste("Spearman R =", round(r_val, 3), ", p =", signif(p_val, 3)),
      x = "WNT7B Expression (Rank)",
      y = paste(g, "Expression (Rank)")
    ) +
    theme_minimal()

  ggsave(filename = paste0("/vol/RZRZK/output/stomach_spearman/Correlation_ScatterPlots_rank/", g, "_ScatterPlot.pdf"),
         plot = p_scatter, width = 6, height = 5)
}

corr_results_df <- bind_rows(corr_results_list)
write_csv(corr_results_df, "/vol/RZRZK/output/stomach/stomach_WNT7B_correlation_result_rank.csv")
