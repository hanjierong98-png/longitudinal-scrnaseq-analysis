suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
})

# --------------------------------------------------
# 1. Load single-cell object
# --------------------------------------------------
dir.create("results", showWarnings = FALSE)
obj <- readRDS("data/single_cell_object.rds")
DefaultAssay(obj) <- "RNA"

# --------------------------------------------------
# 2. Prepare metadata
# --------------------------------------------------
meta <- obj@meta.data

# adapt these column names to your object
stopifnot(all(c("celltype", "timepoint", "patient_id") %in% colnames(obj@meta.data)))

obj$celltype   <- as.character(obj$celltype)
obj$timepoint  <- factor(obj$timepoint, levels = c("T0", "T1", "T2", "T3"))
obj$patient_id <- as.character(obj$patient_id)

obj <- subset(
  obj,
  subset = !is.na(celltype) &
           !is.na(timepoint) &
           !is.na(patient_id)
)

# --------------------------------------------------
# 3. Paired MAST differential expression by cell type
# --------------------------------------------------
run_mast_paired <- function(seu, celltype_name, group1, group2,
                            min_cells = 20, logfc.threshold = 0) {

  seu_ct <- subset(seu, subset = celltype == celltype_name & timepoint %in% c(group1, group2))
  if (ncol(seu_ct) < min_cells) return(NULL)

  # keep only paired patients
  paired_ids <- seu_ct@meta.data %>%
    distinct(patient_id, timepoint) %>%
    count(patient_id) %>%
    filter(n == 2) %>%
    pull(patient_id)

  seu_ct <- subset(seu_ct, subset = patient_id %in% paired_ids)
  if (ncol(seu_ct) < min_cells) return(NULL)

  Idents(seu_ct) <- "timepoint"

  res <- FindMarkers(
    object = seu_ct,
    ident.1 = group1,
    ident.2 = group2,
    test.use = "MAST",
    latent.vars = "patient_id",
    logfc.threshold = logfc.threshold,
    min.pct = 0.05
  )

  if (nrow(res) == 0) return(NULL)

  res$gene <- rownames(res)
  res$celltype <- celltype_name
  res$comparison <- paste0(group1, "_vs_", group2)
  res
}

celltypes <- sort(unique(obj$celltype))
comparisons <- list(
  c("T0", "T3"),
  c("T1", "T3"),
  c("T2", "T3")
)

deg_list <- list()

for (ct in celltypes) {
  for (cmp in comparisons) {
    g1 <- cmp[1]
    g2 <- cmp[2]
    key <- paste(ct, g1, g2, sep = "__")
    message("Running: ", key)

    deg_list[[key]] <- tryCatch(
      run_mast_paired(obj, ct, g1, g2),
      error = function(e) {
        message("Skipped ", key, ": ", e$message)
        NULL
      }
    )
  }
}

deg_all <- bind_rows(deg_list)

# --------------------------------------------------
# 4. Summarize DEG counts across cell types
# --------------------------------------------------
padj_col  <- "p_val_adj"
logfc_col <- "avg_log2FC"

deg_count <- deg_all %>%
  filter(.data[[padj_col]] <= 0.05, abs(.data[[logfc_col]]) >= 0.25) %>%
  count(celltype, comparison, name = "n_deg") %>%
  arrange(desc(n_deg))

print(deg_count)

# optional visualization
ggplot(deg_count, aes(x = comparison, y = n_deg, fill = celltype)) +
  geom_col(position = "dodge") +
  theme_bw() +
  labs(
    title = "Differentially expressed genes across cell types",
    x = "Comparison",
    y = "Number of DEGs"
  )

# --------------------------------------------------
# 5. Save results for downstream analysis
# --------------------------------------------------
saveRDS(deg_all, "results/deg_all_celltypes.rds")
write.csv(deg_count, "results/deg_count_by_celltype.csv", row.names = FALSE)
