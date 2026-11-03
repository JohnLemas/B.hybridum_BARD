################################################################################
# Robustness of the drought / pinoxaden results to the counting strategy
#   fractional multimapper-inclusive (primary)  vs  unique-only (sensitivity)
#
# Needs both pipelines already run:
#   03_outputs/05_overlap        (fractional)   <- DE_overlap_upset_v2.R
#   03_outputs/05_overlap_uniq   (unique-only)  <- same script, pointed at *_uniq
#
# Outputs (03_outputs/07_robustness):
#   class_crosstab.csv            fractional class x unique class, gene counts
#   class_concordance.csv         per class: retained / Jaccard between methods
#   contrast_concordance.csv      per contrast: LFC correlation + DE-call overlap
#   DE_membership_annotated_supported.csv   fractional table + support flags
#   support_by_class.csv, multimap_by_class.csv (if unique counts supplied)
#   uniq_only_genes.csv           genes called DE only with unique counts
################################################################################

suppressPackageStartupMessages({
  library(DESeq2)      # BEFORE tidyverse (masking); also needed to read DESeqResults RDS
  library(tidyverse)
})

# ============================== CONFIG ========================================
setwd("/Users/jmlemas/Library/CloudStorage/OneDrive-MichiganStateUniversity/Dissertation/Chapter3/B.hybridum_BARD")

frac_dir <- "03_outputs/05_overlap"
uniq_dir <- "03_outputs/05_overlap_uniq"
d_f <- "03_outputs/04_DESeq2/Drought_results"
d_u <- "03_outputs/04_DESeq2/Drought_results_uniq"
p_f <- "03_outputs/04_DESeq2/Pinoxaden_results"
p_u <- "03_outputs/04_DESeq2/Pinoxaden_results_uniq"
out_dir <- "03_outputs/07_robustness"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Raw featureCounts tables, used only for the per-gene multimapping fraction.
# Set uniq_counts to your unique-only table (NA = skip that part).
frac_counts <- "03_outputs/03_featureCounts/counts_multi.txt"
uniq_counts <- "03_outputs/03_featureCounts/counts_unique.txt"
high_multi  <- 0.3    # genes with > 30% of reads multimapped are flagged

# Match primary_tier in DE_overlap_upset_v2.R
tier <- list(padj = 0.05, lfc = 1)
# ==============================================================================


# ---- 1. Gene-level class comparison ------------------------------------------
frac <- read_csv(file.path(frac_dir, "DE_membership_annotated.csv"), show_col_types = FALSE)
uniq <- read_csv(file.path(uniq_dir, "DE_membership_annotated.csv"), show_col_types = FALSE)

fu <- full_join(
  frac %>% dplyr::select(gene_id, class_frac = class),
  uniq %>% dplyr::select(gene_id, class_uniq = class),
  by = "gene_id"
) %>%
  mutate(class_frac = replace_na(class_frac, "not_DE"),
         class_uniq = replace_na(class_uniq, "not_DE"))

xt <- fu %>%
  count(class_frac, class_uniq) %>%
  pivot_wider(names_from = class_uniq, values_from = n, values_fill = 0)
write_csv(xt, file.path(out_dir, "class_crosstab.csv"))

classes <- setdiff(union(fu$class_frac, fu$class_uniq), "not_DE")
conc <- map_dfr(classes, function(cl) {
  a <- fu$gene_id[fu$class_frac == cl]
  b <- fu$gene_id[fu$class_uniq == cl]
  tibble(
    class = cl,
    n_frac = length(a),
    n_uniq = length(b),
    n_both = length(intersect(a, b)),
    pct_frac_retained = round(100 * n_both / n_frac, 1),
    jaccard = round(n_both / length(union(a, b)), 3),
    n_frac_not_DE_in_uniq = sum(fu$class_frac == cl & fu$class_uniq == "not_DE")
  )
}) %>% arrange(desc(n_frac))
print(conc)
write_csv(conc, file.path(out_dir, "class_concordance.csv"))


# ---- 2. Per-gene multimapping fraction (optional) -----------------------------
read_total <- function(file) {
  # standard featureCounts layout: Geneid Chr Start End Strand Length <samples...>
  x <- read.table(file, header = TRUE, sep = "\t", comment.char = "#",
                  row.names = 1, check.names = FALSE, quote = "")
  rowSums(x[, -(1:5), drop = FALSE])
}

if (!is.na(uniq_counts) && file.exists(uniq_counts)) {
  tf <- read_total(frac_counts)
  tu <- read_total(uniq_counts)
  multi <- tibble(gene_id = names(tf),
                  total_frac = as.numeric(tf),
                  total_uniq = as.numeric(tu[names(tf)])) %>%
    mutate(multimap_frac = if_else(total_frac > 0,
                                   pmax(0, 1 - total_uniq / total_frac), NA_real_))
} else {
  message("uniq_counts not set / not found: skipping multimapping fraction")
  multi <- tibble(gene_id = character(), total_frac = numeric(),
                  total_uniq = numeric(), multimap_frac = numeric())
}


# ---- 3. Support flags on the fractional table ---------------------------------
annotated <- frac %>%
  left_join(
    uniq %>% dplyr::select(gene_id, class_uniq = class,
                           drought_uniq = drought_status,
                           pinoxaden_uniq = pinoxaden_status,
                           n_sets_strict_uniq = n_sets_strict),
    by = "gene_id") %>%
  mutate(class_uniq = replace_na(class_uniq, "not_DE")) %>%
  left_join(multi, by = "gene_id") %>%
  mutate(
    high_multimap = multimap_frac > high_multi,
    support = case_when(
      class == class_uniq                              ~ "supported_both",
      class_uniq == "not_DE" & high_multimap %in% TRUE ~ "fractional_only_multimapper_heavy",
      class_uniq == "not_DE"                           ~ "fractional_only",
      TRUE                                             ~ "class_differs"
    )
  ) %>%
  relocate(support, class_uniq, multimap_frac, .after = class)

write_csv(annotated, file.path(out_dir, "DE_membership_annotated_supported.csv"))

support_by_class <- annotated %>%
  count(class, support) %>%
  pivot_wider(names_from = support, values_from = n, values_fill = 0)
print(support_by_class)
write_csv(support_by_class, file.path(out_dir, "support_by_class.csv"))

if (any(!is.na(annotated$multimap_frac))) {
  annotated %>%
    group_by(class) %>%
    summarise(n = n(),
              median_multimap_frac = median(multimap_frac, na.rm = TRUE),
              pct_high_multimap = round(100 * mean(high_multimap, na.rm = TRUE), 1),
              .groups = "drop") %>%
    write_csv(file.path(out_dir, "multimap_by_class.csv"))
}

# Genes DE only with unique counts (possibly diluted by fractional counting)
uniq %>%
  anti_join(frac, by = "gene_id") %>%
  write_csv(file.path(out_dir, "uniq_only_genes.csv"))


# ---- 4. Contrast-level concordance --------------------------------------------
sets <- tribble(
  ~set,        ~path_f,                                          ~path_u,
  "Drt_all",   file.path(frac_dir, "cache_Drt_all.rds"),          file.path(uniq_dir, "cache_Drt_all.rds"),
  "Drt_D",     file.path(d_f, "resLFC_trt_drt.rds"),              file.path(d_u, "resLFC_trt_drt.rds"),
  "Drt_H",     file.path(d_f, "resLFC_trt_hrb.rds"),              file.path(d_u, "resLFC_trt_hrb.rds"),
  "Pin_all",   file.path(frac_dir, "cache_Pin_all.rds"),          file.path(uniq_dir, "cache_Pin_all.rds"),
  "Pin_D",     file.path(p_f, "resLFC_trt_drt.rds"),              file.path(p_u, "resLFC_trt_drt.rds"),
  "Pin_H",     file.path(p_f, "resLFC_trt_hrb.rds"),              file.path(p_u, "resLFC_trt_hrb.rds"),
  "DrtInt_D",  file.path(d_f, "resLFC_interaction_drt.rds"),      file.path(d_u, "resLFC_interaction_drt.rds"),
  "DrtInt_H",  file.path(d_f, "resLFC_interaction_hrb.rds"),      file.path(d_u, "resLFC_interaction_hrb.rds"),
  "PinInt_D",  file.path(p_f, "resLFC_interaction_drt.rds"),      file.path(p_u, "resLFC_interaction_drt.rds"),
  "PinInt_H",  file.path(p_f, "resLFC_interaction_hrb.rds"),      file.path(p_u, "resLFC_interaction_hrb.rds"),
  "PinBase_D", file.path(p_f, "resLFC_baseline_drt.rds"),         file.path(p_u, "resLFC_baseline_drt.rds"),
  "PinBase_H", file.path(p_f, "resLFC_baseline_hrb.rds"),         file.path(p_u, "resLFC_baseline_hrb.rds")
)
sets <- sets %>% filter(file.exists(path_f), file.exists(path_u))

read_res <- function(file) {
  as.data.frame(readRDS(file)) %>%
    rownames_to_column("gene_id") %>%
    dplyr::select(gene_id, lfc = log2FoldChange, padj)
}

de_call <- function(lfc, padj) {
  case_when(padj < tier$padj & lfc >  tier$lfc ~ "up",
            padj < tier$padj & lfc < -tier$lfc ~ "down",
            TRUE ~ "none")
}

contrast_conc <- pmap_dfr(sets, function(set, path_f, path_u) {
  d <- inner_join(read_res(path_f), read_res(path_u), by = "gene_id", suffix = c("_f", "_u")) %>%
    filter(!is.na(padj_f), !is.na(padj_u)) %>%
    mutate(call_f = de_call(lfc_f, padj_f), call_u = de_call(lfc_u, padj_u))
  de_f <- d$call_f != "none"
  de_u <- d$call_u != "none"
  tibble(
    set = set,
    n_tested_both = nrow(d),
    rho_lfc = cor(d$lfc_f, d$lfc_u, method = "spearman"),
    n_DE_frac = sum(de_f),
    n_DE_uniq = sum(de_u),
    n_DE_both_same_dir = sum(d$call_f == d$call_u & de_f),
    frac_DE_retained = n_DE_both_same_dir / n_DE_frac,
    jaccard = n_DE_both_same_dir / (n_DE_frac + n_DE_uniq - n_DE_both_same_dir)
  )
})
print(contrast_conc)
write_csv(contrast_conc, file.path(out_dir, "contrast_concordance.csv"))

message("Done. Outputs in ", out_dir)