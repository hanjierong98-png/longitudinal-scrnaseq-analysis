suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
})

# --------------------------------------------------
# 1. Parameters
# --------------------------------------------------
dir.create("results", showWarnings = FALSE)
PADJ_SHOW <- 0.05
TOP_N <- 5

# --------------------------------------------------
# 2. Load module assignments
# --------------------------------------------------
module_df <- read.csv("results/monocyte_modules.csv", stringsAsFactors = FALSE)

if (!all(c("gene", "module") %in% colnames(module_df))) {
  stop("module_df must contain columns: gene, module")
}

# --------------------------------------------------
# 3. Helper functions
# --------------------------------------------------
map_symbols_to_entrez <- function(symbols) {
  symbols <- unique(symbols[!is.na(symbols) & symbols != ""])
  if (length(symbols) == 0) return(character(0))

  mp <- suppressWarnings(
    bitr(
      symbols,
      fromType = "SYMBOL",
      toType   = "ENTREZID",
      OrgDb    = org.Hs.eg.db
    )
  )

  if (is.null(mp) || nrow(mp) == 0) return(character(0))
  unique(mp$ENTREZID)
}

normalize_geneid <- function(x) {
  x <- as.character(x)
  empty <- is.na(x) | x == ""
  x[empty] <- NA_character_
  x <- stringr::str_replace_all(x, "[,;\\s]+", "/")
  x <- stringr::str_replace_all(x, "/+", "/")
  x <- stringr::str_replace_all(x, "^/|/$", "")
  x
}

format_ora_df <- function(enrich_obj, module_id, set_name, padj_filter = NULL) {
  if (is.null(enrich_obj)) return(NULL)

  df <- as.data.frame(enrich_obj)
  if (is.null(df) || nrow(df) == 0) return(NULL)

  if (!"p.adjust" %in% colnames(df)) df$p.adjust <- NA_real_
  if (!"Count" %in% colnames(df)) df$Count <- NA_real_
  if (!"BgRatio" %in% colnames(df)) df$BgRatio <- NA_character_
  if (!"geneID" %in% colnames(df)) df$geneID <- NA_character_

  suppressWarnings(df$p.adjust <- as.numeric(df$p.adjust))
  suppressWarnings(df$Count <- as.numeric(df$Count))

  df$TermSize <- suppressWarnings(
    as.numeric(stringr::str_extract(df$BgRatio, "^[0-9]+"))
  )

  df$Gene_ratio <- NA_real_
  ok <- !is.na(df$Count) & !is.na(df$TermSize) & df$TermSize > 0
  df$Gene_ratio[ok] <- (df$Count[ok] / df$TermSize[ok]) * 100

  df_out <- df %>%
    mutate(
      Module = module_id,
      Set = set_name,
      geneID = as.character(.data$geneID),
      geneID_list = normalize_geneid(.data$geneID)
    ) %>%
    select(Module, Set, Description, p.adjust, Count, geneID, geneID_list, TermSize, Gene_ratio) %>%
    filter(!is.na(p.adjust), !is.na(Gene_ratio), Gene_ratio > 0)

  if (!is.null(padj_filter)) {
    df_out <- df_out %>% filter(p.adjust <= padj_filter)
  }

  if (nrow(df_out) == 0) return(NULL)
  df_out
}

run_module_ora <- function(genes) {
  entrez <- map_symbols_to_entrez(genes)

  ego <- try(
    enrichGO(
      gene          = entrez,
      OrgDb         = org.Hs.eg.db,
      keyType       = "ENTREZID",
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 1,
      qvalueCutoff  = 1,
      readable      = TRUE
    ),
    silent = TRUE
  )
  ego <- if (inherits(ego, "try-error")) NULL else ego

  react <- try(
    enrichPathway(
      gene          = entrez,
      organism      = "human",
      pAdjustMethod = "BH",
      pvalueCutoff  = 1,
      qvalueCutoff  = 1,
      readable      = TRUE
    ),
    silent = TRUE
  )
  react <- if (inherits(react, "try-error")) NULL else react

  list(GO_BP = ego, Reactome = react, ENTREZ = entrez)
}

choose_topN_sig <- function(df, top_n = 2, padj_show = 0.05) {
  df %>%
    filter(!is.na(p.adjust), p.adjust <= padj_show, !is.na(Gene_ratio), Gene_ratio > 0) %>%
    group_by(Module, Set) %>%
    arrange(p.adjust, desc(Gene_ratio), desc(Count)) %>%
    slice_head(n = top_n) %>%
    ungroup()
}

# --------------------------------------------------
# 4. Run ORA by module
# --------------------------------------------------
dir.create("results", showWarnings = FALSE)

module_genes <- split(module_df$gene, module_df$module)

module_gene_counts <- list()
all_full <- list()
all_sig <- list()

for (mod in names(module_genes)) {
  genes_sym <- unique(module_genes[[mod]])
  ora_res <- run_module_ora(genes_sym)

  module_gene_counts[[length(module_gene_counts) + 1]] <- tibble(
    Module = mod,
    n_symbol = length(genes_sym),
    n_entrez = length(ora_res$ENTREZ)
  )

  df_go_full <- format_ora_df(ora_res$GO_BP, mod, "GO_BP")
  df_re_full <- format_ora_df(ora_res$Reactome, mod, "Reactome")
  df_full <- bind_rows(df_go_full, df_re_full)

  if (!is.null(df_full) && nrow(df_full) > 0) {
    write_csv(df_full, file.path("results", paste0(mod, "_ORA_full.csv")))
    all_full[[length(all_full) + 1]] <- df_full
  }

  df_go_sig <- format_ora_df(ora_res$GO_BP, mod, "GO_BP", padj_filter = PADJ_SHOW)
  df_re_sig <- format_ora_df(ora_res$Reactome, mod, "Reactome", padj_filter = PADJ_SHOW)
  df_sig <- bind_rows(df_go_sig, df_re_sig)

  if (!is.null(df_sig) && nrow(df_sig) > 0) {
    write_csv(df_sig, file.path("results", paste0(mod, "_ORA_sig.csv")))
    all_sig[[length(all_sig) + 1]] <- df_sig
  }
}

# --------------------------------------------------
# 5. Save combined tables
# --------------------------------------------------
write_csv(bind_rows(module_gene_counts), "results/module_gene_counts.csv")

ora_all_full <- bind_rows(all_full)
ora_all_sig <- bind_rows(all_sig)

write_csv(ora_all_full, "results/module_ORA_all_full.csv")
write_csv(ora_all_sig, "results/module_ORA_all_sig.csv")

# --------------------------------------------------
# 6. Visualization: Top N significant terms per Module × Set
# --------------------------------------------------
ora_top <- choose_topN_sig(ora_all_sig, top_n = TOP_N, padj_show = PADJ_SHOW) %>%
  mutate(
    Set = factor(Set, levels = c("GO_BP", "Reactome")),
    Module = factor(Module, levels = unique(module_df$module))
  )

p <- ggplot(
  ora_top,
  aes(x = Set, y = forcats::fct_rev(stringr::str_wrap(Description, 46)))
) +
  geom_point(aes(size = Gene_ratio, colour = p.adjust), alpha = 0.9) +
  facet_grid(Module ~ ., scales = "free_y", space = "free") +
  scale_colour_gradient(
    low = "#67000D",
    high = "#FDE0DD"
  ) +
  theme_bw(base_size = 14) +
  theme(
    text = element_text(colour = "black"),
    axis.text = element_text(colour = "black"),
    strip.text = element_text(face = "bold", colour = "black"),
    axis.text.y = element_text(size = 11)
  ) +
  labs(
    title = paste0(
      "Module enrichment analysis (padj ≤ ", PADJ_SHOW,
      ", Top ", TOP_N, " terms per Module × Set)"
    ),
    x = NULL,
    y = NULL,
    colour = "Adjusted p-value",
    size = "Gene ratio (%)"
  )

ggsave(
  filename = "results/module_ORA_top_terms.png",
  plot = p,
  width = 8,
  height = 10,
  dpi = 300
)
