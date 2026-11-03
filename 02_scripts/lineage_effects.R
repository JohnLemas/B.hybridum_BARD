################################################################################
# Did recurrent selection alter the transcriptome?  Lineage-effect checks.
#
#  1. Threshold grid: how many genes are called in every interaction/baseline
#     contrast across padj x |LFC| cutoffs (+ p-value histograms). Shows whether
#     lineage effects exist but are small (hidden by a strict |LFC| filter).
#  2. Priming test: do genes shifted at BASELINE in a selected lineage overlap /
#     correlate with genes that respond to stress?
#  3. Gene-set shifts: baseline LFC of shared_up / shared_down / drought_specific
#     genes vs all other genes.
#  4. PCA: lineage vs treatment variance in each experiment, plus PCA after
#     removing the treatment effect.
#
# Set run <- "uniq" to repeat everything on the unique-only results.
#
# CAUTION (priming test): contrasts from the SAME experiment share samples (e.g.
# well-watered plants contribute to both the baseline and the drought-response
# estimate), which induces correlated noise. Cross-experiment pairs (e.g.
# PinBase_H vs Drt_all) use independent plants and are the cleaner test; the
# output flags them with same_experiment = FALSE.
# CAUTION (interpretation): a baseline difference between lines can also reflect
# genetic divergence among the lines, not only selection.
################################################################################

suppressPackageStartupMessages({
  library(DESeq2)      # BEFORE tidyverse (masking)
  library(tidyverse)
})

# ============================== CONFIG ========================================
setwd("/Users/jmlemas/Library/CloudStorage/OneDrive-MichiganStateUniversity/Dissertation/Chapter3/B.hybridum_BARD")

run <- "frac"                                   # "frac" or "uniq"
sfx <- if (run == "uniq") "_uniq" else ""

d_dir       <- paste0("03_outputs/04_DESeq2/Drought_results", sfx)
p_dir       <- paste0("03_outputs/04_DESeq2/Pinoxaden_results", sfx)
overlap_dir <- paste0("03_outputs/05_overlap", sfx)     # holds cache_Drt_all.rds / cache_Pin_all.rds
robust_dir  <- "03_outputs/07_robustness"               # supported_both classes (fractional run)
out_dir     <- paste0("03_outputs/09_lineage_effects", sfx)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tier <- list(padj = 0.05, lfc = 1)       # call used for baseline / response DE in sections 2-3

grid_padj <- c(0.1, 0.05, 0.01)          # section 1 grid
grid_lfc  <- c(0, 0.5, 1, 2)

ntop <- 1000                             # most-variable genes used for PCA

pca_cfg <- list(
  drought   = list(dds = file.path(d_dir, "dds_full.rds"), lineage = "lineage", trt = "trt"),
  pinoxaden = list(dds = file.path(p_dir, "dds_full.rds"), lineage = "line",    trt = "trt")
)

# Same labels/colours as shared_heatmap.R
lineage_map   <- c(C = "Control", Control = "Control", D = "Drought", Drought = "Drought",
                   H = "Herbicide", Herbicide = "Herbicide")
treatment_map <- c(WW = "Untreated / well-watered", Control = "Untreated / well-watered",
                   DW = "Drought", Pinoxaden = "Pinoxaden")
lineage_cols  <- c(Control = "#1B9E77", Drought = "#7570B3", Herbicide = "#E7298A")

sets <- tribble(
  ~set,        ~experiment, ~type,         ~path,
  "Drt_all",   "drought",   "response",    file.path(overlap_dir, "cache_Drt_all.rds"),
  "Pin_all",   "pinoxaden", "response",    file.path(overlap_dir, "cache_Pin_all.rds"),
  "DrtInt_D",  "drought",   "interaction", file.path(d_dir, "resLFC_interaction_drt.rds"),
  "DrtInt_H",  "drought",   "interaction", file.path(d_dir, "resLFC_interaction_hrb.rds"),
  "PinInt_D",  "pinoxaden", "interaction", file.path(p_dir, "resLFC_interaction_drt.rds"),
  "PinInt_H",  "pinoxaden", "interaction", file.path(p_dir, "resLFC_interaction_hrb.rds"),
  "DrtBase_D", "drought",   "baseline",    file.path(d_dir, "resLFC_baseline_drt.rds"),
  "DrtBase_H", "drought",   "baseline",    file.path(d_dir, "resLFC_baseline_hrb.rds"),
  "PinBase_D", "pinoxaden", "baseline",    file.path(p_dir, "resLFC_baseline_drt.rds"),
  "PinBase_H", "pinoxaden", "baseline",    file.path(p_dir, "resLFC_baseline_hrb.rds")
)
# ==============================================================================


# ---- 0. Load results ----------------------------------------------------------
missing <- !file.exists(sets$path)
if (any(missing)) message("Skipping (file not found): ", paste(sets$set[missing], collapse = ", "))
sets <- sets[!missing, ]

read_res <- function(f) {
  as.data.frame(readRDS(f)) %>%
    rownames_to_column("gene_id") %>%
    dplyr::select(gene_id, lfc = log2FoldChange, padj, any_of("pvalue"))
}

res_all <- sets %>%
  mutate(data = map(path, read_res)) %>%
  dplyr::select(set, experiment, type, data) %>%
  unnest(data)

res_lin <- res_all %>% filter(type %in% c("interaction", "baseline"))

de_call <- function(lfc, padj) {
  case_when(padj < tier$padj & lfc >  tier$lfc ~ "up",
            padj < tier$padj & lfc < -tier$lfc ~ "down",
            TRUE ~ "none")
}


# ---- 1. Threshold grid + p-value histograms ----------------------------------------
thr <- expand_grid(padj_cut = grid_padj, lfc_cut = grid_lfc)

lin_counts <- res_lin %>%
  group_by(set) %>%
  group_split() %>%
  map_dfr(function(d) {
    thr %>% mutate(
      set      = d$set[1],
      n_tested = sum(!is.na(d$padj)),
      n_up     = map2_int(padj_cut, lfc_cut, function(p, l) sum(d$padj < p & d$lfc >  l, na.rm = TRUE)),
      n_down   = map2_int(padj_cut, lfc_cut, function(p, l) sum(d$padj < p & d$lfc < -l, na.rm = TRUE)),
      n_total  = n_up + n_down
    )
  }) %>%
  relocate(set)
write_csv(lin_counts, file.path(out_dir, "lineage_threshold_grid.csv"))

p1 <- ggplot(lin_counts, aes(factor(lfc_cut), factor(padj_cut), fill = log10(n_total + 1))) +
  geom_tile(colour = "white") +
  geom_text(aes(label = n_total), size = 3) +
  facet_wrap(~ set, ncol = 4) +
  scale_fill_viridis_c(name = "log10(n+1)", direction = -1) +
  labs(x = "|shrunk log2FC| cutoff", y = "padj cutoff",
       title = "Genes called in each lineage-effect contrast") +
  theme_bw()
ggsave(file.path(out_dir, "lineage_threshold_grid.pdf"), p1, width = 12, height = 6)

if ("pvalue" %in% names(res_lin)) {
  p2 <- res_lin %>% filter(!is.na(pvalue)) %>%
    ggplot(aes(pvalue)) +
    geom_histogram(bins = 40, boundary = 0, fill = "grey30") +
    facet_wrap(~ set, scales = "free_y", ncol = 4) +
    labs(title = "p-value distributions (a pile-up near 0 = real signal; flat = none)") +
    theme_bw()
  ggsave(file.path(out_dir, "lineage_pvalue_histograms.pdf"), p2, width = 12, height = 6)
}


# ---- 2. Priming: baseline shifts vs stress response -----------------------------------
base_sets <- sets$set[sets$type == "baseline"]
resp_sets <- sets$set[sets$type == "response"]
exp_of    <- setNames(sets$experiment, sets$set)

pairs <- expand_grid(base = base_sets, resp = resp_sets)

pair_data <- pmap(pairs, function(base, resp) {
  inner_join(
    res_all %>% filter(set == base, !is.na(padj)) %>% dplyr::select(gene_id, lfc_b = lfc, padj_b = padj),
    res_all %>% filter(set == resp, !is.na(padj)) %>% dplyr::select(gene_id, lfc_r = lfc, padj_r = padj),
    by = "gene_id") %>%
    mutate(call_b = de_call(lfc_b, padj_b), call_r = de_call(lfc_r, padj_r),
           base = base, resp = resp)
})

overlap_tab <- map2_dfr(pair_data, seq_along(pair_data), function(d, i) {
  expand_grid(bdir = c("up", "down"), rdir = c("up", "down")) %>%
    pmap_dfr(function(bdir, rdir) {
      b <- d$call_b == bdir
      r <- d$call_r == rdir
      ft <- tryCatch(fisher.test(table(factor(b, c(TRUE, FALSE)), factor(r, c(TRUE, FALSE)))),
                     error = function(e) NULL)
      tibble(base = d$base[1], resp = d$resp[1],
             same_experiment = exp_of[[d$base[1]]] == exp_of[[d$resp[1]]],
             n_tested = nrow(d), baseline_dir = bdir, response_dir = rdir,
             n_base = sum(b), n_resp = sum(r), n_overlap = sum(b & r),
             expected = sum(b) * sum(r) / nrow(d),
             odds_ratio = if (is.null(ft)) NA_real_ else unname(ft$estimate),
             p = if (is.null(ft)) NA_real_ else ft$p.value)
    })
}) %>%
  mutate(p_adj = p.adjust(p, method = "BH"),
         pattern = if_else(baseline_dir == response_dir,
                           "priming-like (baseline shifted toward response)",
                           "opposite"))
write_csv(overlap_tab, file.path(out_dir, "priming_overlap.csv"))

safe_cor <- function(x, y) if (length(x) >= 3) suppressWarnings(cor(x, y, method = "spearman")) else NA_real_
cor_tab <- map_dfr(pair_data, function(d) {
  de_r <- d$call_r != "none"
  tibble(base = d$base[1], resp = d$resp[1],
         same_experiment = exp_of[[d$base[1]]] == exp_of[[d$resp[1]]],
         n_tested = nrow(d),
         rho_all_genes = safe_cor(d$lfc_b, d$lfc_r),
         n_response_DE = sum(de_r),
         rho_response_DE_genes = safe_cor(d$lfc_b[de_r], d$lfc_r[de_r]))
})
print(cor_tab)
write_csv(cor_tab, file.path(out_dir, "priming_correlation.csv"))

scatter_df <- bind_rows(pair_data) %>%
  mutate(pair = paste(base, "vs", resp), resp_DE = call_r != "none")
p3 <- ggplot(scatter_df, aes(lfc_b, lfc_r)) +
  geom_hline(yintercept = 0, colour = "grey70") + geom_vline(xintercept = 0, colour = "grey70") +
  geom_point(data = filter(scatter_df, !resp_DE), colour = "grey70", size = 0.4, alpha = 0.4) +
  geom_point(data = filter(scatter_df, resp_DE), colour = "firebrick", size = 0.9, alpha = 0.7) +
  facet_wrap(~ pair, scales = "free", ncol = 4) +
  labs(x = "Baseline log2FC (selected lineage vs Control, untreated)",
       y = "Stress-response log2FC", title = "Baseline shift vs stress response (red = response DE)") +
  theme_bw()
ggsave(file.path(out_dir, "priming_scatter.pdf"), p3,
       width = 14, height = 3.6 * ceiling(nrow(pairs) / 4))


# ---- 3. Gene-set shifts at baseline -----------------------------------------------------
sup_file <- file.path(robust_dir, "DE_membership_annotated_supported.csv")
if (file.exists(sup_file)) {
  cls <- read_csv(sup_file, show_col_types = FALSE, guess_max = Inf) %>%
    filter(support == "supported_both", class %in% c("shared_up", "shared_down", "drought_specific")) %>%
    dplyr::select(gene_id, class)

  shift_tab <- expand_grid(base = base_sets, class = unique(cls$class)) %>%
    pmap_dfr(function(base, class) {
      d <- res_all %>% filter(set == base, !is.na(padj))
      g <- cls$gene_id[cls$class == class]
      inn <- d$gene_id %in% g
      if (sum(inn) < 5) return(NULL)
      tibble(base = base, class = class, n_genes = sum(inn),
             median_lfc_in = median(d$lfc[inn]),
             median_lfc_other = median(d$lfc[!(d$gene_id %in% cls$gene_id)]),
             p_wilcox = wilcox.test(d$lfc[inn], d$lfc[!inn])$p.value)
    }) %>%
    mutate(p_adj = p.adjust(p_wilcox, method = "BH"))
  print(shift_tab)
  write_csv(shift_tab, file.path(out_dir, "geneset_baseline_shift.csv"))

  box_df <- map_dfr(base_sets, function(b) {
    d <- res_all %>% filter(set == b, !is.na(padj))
    bind_rows(
      map_dfr(unique(cls$class), function(cl)
        tibble(set = b, group = cl, lfc = d$lfc[d$gene_id %in% cls$gene_id[cls$class == cl]])),
      tibble(set = b, group = "all other genes", lfc = d$lfc[!(d$gene_id %in% cls$gene_id)])
    )
  })
  p4 <- ggplot(box_df, aes(group, lfc, fill = group)) +
    geom_hline(yintercept = 0, colour = "grey60") +
    geom_boxplot(outlier.size = 0.3, show.legend = FALSE) +
    coord_cartesian(ylim = c(-1.5, 1.5)) +
    facet_wrap(~ set, ncol = 4) +
    labs(x = NULL, y = "Baseline log2FC (selected lineage vs Control)",
         title = "Do shared / specific genes shift at baseline?") +
    theme_bw() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(file.path(out_dir, "geneset_baseline_boxplot.pdf"), p4, width = 12, height = 6)
} else {
  message("Skipping section 3: ", sup_file, " not found")
}


# ---- 4. PCA: lineage vs treatment variance ------------------------------------------------
pca_res <- imap(pca_cfg, function(cfg, nm) {
  dds <- readRDS(cfg$dds)
  cd  <- as.data.frame(colData(dds))
  stopifnot(all(c(cfg$lineage, cfg$trt) %in% names(cd)))
  lin <- unname(lineage_map[as.character(cd[[cfg$lineage]])])
  trt <- unname(treatment_map[as.character(cd[[cfg$trt]])])
  if (anyNA(lin) || anyNA(trt)) stop(nm, ": unmapped lineage/treatment values; edit lineage_map / treatment_map")

  v   <- assay(vst(dds, blind = FALSE))
  top <- order(matrixStats::rowVars(v), decreasing = TRUE)[seq_len(min(ntop, nrow(v)))]
  mat <- t(v[top, ])

  # remove the treatment effect: centre each gene within each treatment group
  mat_c <- mat
  for (lv in unique(trt)) {
    idx <- trt == lv
    mat_c[idx, ] <- sweep(mat[idx, , drop = FALSE], 2, colMeans(mat[idx, , drop = FALSE]))
  }

  do_pca <- function(x, label) {
    pr  <- prcomp(x, center = TRUE, scale. = FALSE)
    pct <- 100 * pr$sdev^2 / sum(pr$sdev^2)
    list(
      scores = tibble(sample = rownames(x), PC1 = pr$x[, 1], PC2 = pr$x[, 2],
                      lineage = lin, treatment = trt,
                      panel = sprintf("%s, %s (PC1 %.0f%%, PC2 %.0f%%)",
                                      str_to_title(nm), label, pct[1], pct[2])),
      vp = map_dfr(1:5, function(i) tibble(
        experiment = nm, data = label, pc = i, pct_var = pct[i],
        r2_lineage   = summary(lm(pr$x[, i] ~ lin))$r.squared,
        r2_treatment = summary(lm(pr$x[, i] ~ trt))$r.squared))
    )
  }
  list(raw = do_pca(mat, "all samples"), centred = do_pca(mat_c, "treatment effect removed"))
})

scores <- map_dfr(pca_res, function(x) bind_rows(x$raw$scores, x$centred$scores))
vp     <- map_dfr(pca_res, function(x) bind_rows(x$raw$vp, x$centred$vp))
print(vp, n = Inf)
write_csv(vp, file.path(out_dir, "pca_variance_partition.csv"))

p5 <- ggplot(scores, aes(PC1, PC2, colour = lineage, shape = treatment)) +
  geom_point(size = 2, alpha = 0.8) +
  facet_wrap(~ panel, scales = "free", ncol = 2) +
  scale_colour_manual(values = lineage_cols, name = "Lineage") +
  labs(shape = "Treatment") +
  theme_bw()
ggsave(file.path(out_dir, "pca_lineage_treatment.pdf"), p5, width = 12, height = 9)

message("Done. Outputs in ", out_dir)