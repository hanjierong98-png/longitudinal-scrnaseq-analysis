suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
})

# --------------------------------------------------
# 1. Load object and module assignments
# --------------------------------------------------
obj <- readRDS("single_cell_object.rds")
DefaultAssay(obj) <- "RNA"

module_df <- read.csv("results/monocyte_modules.csv", stringsAsFactors = FALSE)

stopifnot(all(c("gene", "module") %in% colnames(module_df)))
stopifnot(all(c("celltype", "timepoint") %in% colnames(obj@meta.data)))

obj$timepoint <- factor(obj$timepoint, levels = c("T0", "T1", "T2", "T3"))

# --------------------------------------------------
# 2. Focus on monocytes
# --------------------------------------------------
obj_mono <- subset(obj, subset = celltype == "Monocytes")

expr_mat <- GetAssayData(obj_mono, slot = "data")

# --------------------------------------------------
# 3. Calculate average expression per module and time point
# --------------------------------------------------
module_list <- split(module_df$gene, module_df$module)

avg_expr_list <- lapply(names(module_list), function(mod) {
  genes_use <- intersect(module_list[[mod]], rownames(expr_mat))
  if (length(genes_use) == 0) return(NULL)

  df <- data.frame(
    timepoint = obj_mono$timepoint,
    expr = Matrix::colMeans(expr_mat[genes_use, , drop = FALSE])
  ) %>%
    group_by(timepoint) %>%
    summarise(mean_expression = mean(expr), .groups = "drop") %>%
    mutate(module = mod)

  df
})

module_trend_df <- bind_rows(avg_expr_list)

# --------------------------------------------------
# 4. Plot module trends
# --------------------------------------------------
p <- ggplot(module_trend_df, aes(x = timepoint, y = mean_expression, group = module)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~ module, scales = "free_y") +
  theme_bw(base_size = 14) +
  labs(
    title = "Average module expression across time points",
    x = "Time point",
    y = "Mean expression"
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    text = element_text(colour = "black"),
    axis.text = element_text(colour = "black")
  )

ggsave(
  filename = "results/module_trends.png",
  plot = p,
  width = 10,
  height = 6,
  dpi = 300
)
