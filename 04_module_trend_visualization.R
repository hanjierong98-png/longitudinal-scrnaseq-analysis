suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(UCell)
})

# --------------------------------------------------
# 1. Parameters
# --------------------------------------------------
ora_csv <- "results/module_ORA_all_sig.csv"

celltype_col   <- "celltype"
celltype_label <- "Monocytes"
timepoint_col  <- "timepoint"

padj_cut    <- 0.05
top_n       <- 15
maxRank_use <- 5000
min_gs_size <- 5

out_dir <- "results/pathway_auc"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------
# 2. Helper: parse geneID_list
# --------------------------------------------------
split_genes <- function(x) {
  if (is.na(x) || !nzchar(x)) return(character(0))

  x <- stringr::str_replace_all(x, "^c\\(|\\)$", "")
  x <- stringr::str_replace_all(x, "[\"'`]", "")

  parts <- unlist(
    stringr::str_split(x, "\\s*[/;,|\\s]\\s*")
  )

  parts <- parts[nzchar(parts)]
  unique(parts)
}

# --------------------------------------------------
# 3. Load ORA results
# --------------------------------------------------
stopifnot(file.exists(ora_csv))
ora_raw <- readr::read_csv(ora_csv, show_col_types = FALSE)

stopifnot(all(
  c("Module", "Set", "Description", "p.adjust", "Count", "geneID_list") %in% colnames(ora_raw)
))

ora <- ora_raw %>%
  mutate(
    geneID_list = purrr::map(geneID_list, split_genes)
  )

ora_base <- ora %>%
  transmute(
    Module      = as.character(Module),
    Set         = as.character(Set),
    Description = as.character(Description),
    padj        = as.numeric(`p.adjust`),
    Count_ORA   = as.integer(Count),
    geneID_list
  ) %>%
  filter(
    Set %in% c("GO_BP", "Reactome"),
    padj <= padj_cut
  )

readr::write_csv(
  ora_base %>% select(-geneID_list),
  file.path(out_dir, "ORA_sig_filtered.csv")
)

# --------------------------------------------------
# 4. Select top pathways per Module × Set
# --------------------------------------------------
ora_top <- ora_base %>%
  group_by(Module, Set) %>%
  arrange(padj, desc(Count_ORA), .by_group = TRUE) %>%
  slice_head(n = top_n) %>%
  ungroup()

readr::write_csv(
  ora_top %>% select(-geneID_list),
  file.path(out_dir, "ORA_top15_perModule_perSet.csv")
)

# --------------------------------------------------
# 5. Load single-cell object and focus on monocytes
# --------------------------------------------------
obj <- readRDS("data/single_cell_object.rds")
DefaultAssay(obj) <- "RNA"

stopifnot(all(c(celltype_col, timepoint_col) %in% colnames(obj@meta.data)))

Seurat::Idents(obj) <- celltype_col
obj_mono <- subset(obj, idents = celltype_label)

mono_genes <- rownames(obj_mono)

# --------------------------------------------------
# 6. Build pathway gene sets from geneID_list
# --------------------------------------------------
ora_sets <- ora_top %>%
  mutate(
    FeatureID = sprintf("PW%04d", row_number()),
    genes_ORA = purrr::map(geneID_list, ~ unique(.x)),
    genes_use = purrr::map(genes_ORA, ~ intersect(.x, mono_genes)),
    n_genes   = purrr::map_int(genes_use, length)
  ) %>%
  filter(n_genes >= min_gs_size)

stopifnot(nrow(ora_sets) > 0)

features_list <- setNames(
  ora_sets$genes_use,
  ora_sets$FeatureID
)

readr::write_csv(
  ora_sets %>%
    select(
      FeatureID,
      Module,
      Set,
      Description,
      padj,
      Count_ORA,
      n_genes
    ),
  file.path(out_dir, "UCell_feature_table.csv")
)

# --------------------------------------------------
# 7. Run UCell
# --------------------------------------------------
meta_before <- colnames(obj_mono@meta.data)

obj_mono <- UCell::AddModuleScore_UCell(
  obj_mono,
  features = features_list,
  name = "PW_",
  maxRank = maxRank_use
)

meta_after <- colnames(obj_mono@meta.data)
new_cols <- setdiff(meta_after, meta_before)

score_cols <- new_cols[
  sapply(obj_mono@meta.data[, new_cols, drop = FALSE], is.numeric)
]

pattern <- paste0("(", paste(ora_sets$FeatureID, collapse = "|"), ")")

col_map <- tibble::tibble(
  ScoreCol = score_cols
) %>%
  mutate(
    FeatureID = stringr::str_extract(ScoreCol, pattern)
  ) %>%
  left_join(
    ora_sets %>%
      select(FeatureID, Module, Set, Description),
    by = "FeatureID"
  )

readr::write_csv(
  col_map,
  file.path(out_dir, "UCell_score_map.csv")
)

# --------------------------------------------------
# 8. Summarize mean UCell score by time point
# --------------------------------------------------
plot_df <- obj_mono@meta.data %>%
  mutate(
    Timepoint = factor(
      .data[[timepoint_col]],
      levels = c("T0", "T1", "T2", "T3")
    )
  ) %>%
  select(Timepoint, all_of(col_map$ScoreCol)) %>%
  pivot_longer(
    cols = all_of(col_map$ScoreCol),
    names_to = "ScoreCol",
    values_to = "Score"
  ) %>%
  left_join(col_map, by = "ScoreCol") %>%
  group_by(Module, Set, Description, Timepoint) %>%
  summarise(
    mean_score = mean(Score, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(
  plot_df,
  file.path(out_dir, "UCell_mean_by_timepoint.csv")
)

# --------------------------------------------------
# 9. Plot pathway dynamics for each module
# --------------------------------------------------
plot_df$Set <- factor(
  plot_df$Set,
  levels = c("GO_BP", "Reactome")
)

modules <- unique(plot_df$Module)

plot_one <- function(mod) {
  df <- plot_df %>%
    filter(Module == mod)

  ggplot(
    df,
    aes(
      x = Timepoint,
      y = mean_score,
      group = Description,
      color = Description,
      linetype = Set
    )
  ) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.7) +
    scale_linetype_manual(
      values = c(
        GO_BP = "solid",
        Reactome = "dashed"
      )
    ) +
    theme_classic(base_size = 13) +
    theme(
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 9),
      text = element_text(colour = "black"),
      axis.text = element_text(colour = "black")
    ) +
    labs(
      title = paste("Monocyte pathway activity -", mod),
      y = "Mean UCell score",
      x = NULL,
      color = "Pathway",
      linetype = "Set"
    )
}

plots <- lapply(modules, plot_one)
names(plots) <- modules

for (mod in modules) {
  p <- plots[[mod]]

  ggsave(
    file.path(
      out_dir,
      paste0(
        gsub("[^A-Za-z0-9_]+", "_", mod),
        "_UCell_lineplot.png"
      )
    ),
    p,
    width = 9,
    height = 6,
    dpi = 300
  )
}
