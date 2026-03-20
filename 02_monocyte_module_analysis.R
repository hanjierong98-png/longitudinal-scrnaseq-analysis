suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
  library(cluster)
})

# --------------------------------------------------
# 1. Load DEG results
# --------------------------------------------------
deg_all <- readRDS("results/deg_all_celltypes.rds")

padj_candidates  <- c("p_val_adj", "p_adj", "FDR")
logfc_candidates <- c("avg_log2FC", "avg_logFC", "log2FoldChange", "logFC")

padj_col  <- intersect(padj_candidates, colnames(deg_all))[1]
logfc_col <- intersect(logfc_candidates, colnames(deg_all))[1]

# --------------------------------------------------
# 2. Focus on monocytes
# --------------------------------------------------
deg_mono <- deg_all %>%
  filter(celltype == "Monocytes")

# --------------------------------------------------
# 3. Build effect-size matrix
# strong genes: padj <= 0.05 and |logFC| >= 0.25 in any comparison
# display: show raw logFC where padj <= 0.05, otherwise NA
# --------------------------------------------------
strong_genes <- deg_mono %>%
  filter(.data[[padj_col]] <= 0.05, abs(.data[[logfc_col]]) >= 0.25) %>%
  pull(gene) %>%
  unique()

mat_show_df <- deg_mono %>%
  filter(gene %in% strong_genes) %>%
  mutate(value_for_plot = ifelse(.data[[padj_col]] <= 0.05, .data[[logfc_col]], NA_real_)) %>%
  select(gene, comparison, value_for_plot) %>%
  distinct() %>%
  pivot_wider(names_from = comparison, values_from = value_for_plot)

rownames(mat_show_df) <- mat_show_df$gene
mat_show <- as.matrix(mat_show_df[, -1])

# matrix for clustering: fill NA as 0
mat_clust <- mat_show
mat_clust[is.na(mat_clust)] <- 0

# --------------------------------------------------
# 4. Choose number of clusters using silhouette
# --------------------------------------------------
choose_k_silhouette <- function(mat, k_min = 2, k_max = 8) {
  d <- dist(mat)
  hc <- hclust(d, method = "ward.D2")

  sil_df <- lapply(k_min:k_max, function(k) {
    cl <- cutree(hc, k = k)
    sil <- silhouette(cl, d)
    data.frame(k = k, avg_silhouette = mean(sil[, 3]))
  }) %>% bind_rows()

  sil_df
}

sil_df <- choose_k_silhouette(mat_clust, k_min = 2, k_max = 8)
best_k <- sil_df$k[which.max(sil_df$avg_silhouette)]

print(sil_df)
message("Selected k = ", best_k)

# --------------------------------------------------
# 5. Cluster genes into modules
# --------------------------------------------------
hc <- hclust(dist(mat_clust), method = "ward.D2")
clusters <- cutree(hc, k = best_k)

module_df <- data.frame(
  gene = rownames(mat_clust),
  module = paste0("Module_", clusters),
  stringsAsFactors = FALSE
)

row_split <- factor(module_df$module, levels = unique(module_df$module))

# --------------------------------------------------
# 6. Heatmap visualization
# --------------------------------------------------
ht <- Heatmap(
  mat_show,
  name = "logFC",
  col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
  na_col = "white",
  show_row_names = FALSE,
  cluster_rows = as.dendrogram(hc),
  cluster_columns = FALSE,
  row_split = row_split,
  column_title = "Monocyte differential expression modules"
)

draw(ht)

# --------------------------------------------------
# 7. Save outputs
# --------------------------------------------------
write.csv(module_df, "results/monocyte_modules.csv", row.names = FALSE)
saveRDS(mat_show, "results/monocyte_effect_matrix.rds")
write.csv(sil_df, "results/monocyte_silhouette_scores.csv", row.names = FALSE)
