################################################################################
# DE overlap (UpSet) + drought-specific vs shared-stress classification
# + LFC concordance + InterProScan annotation      (v2)
#
# Differences from v1:
#  - Reads the FULL shrunk DESeq2 results (resLFC_*.rds / dds_full.rds), not the
#    pre-filtered sig_*.txt files, so every threshold lives in this CONFIG block.
#  - Main-effect sets (Drt_all, Pin_all) are recomputed from dds_full.rds with the
#    same recipe as every other set (test vs 0 + ashr shrinkage), removing the
#    lfcThreshold = 2 mismatch. Results are cached in out_dir.
#  - Two threshold tiers (lenient / strict). UpSets are drawn for both; the
#    classification uses `primary_tier`; strict hits are flagged in the output.
#  - "Specific" now requires the gene to be CLEARLY not DE in the other
#    experiment (not just non-significant). In-between genes are "ambiguous".
#  - Adds a threshold-free LFC concordance check (drought vs pinoxaden).
################################################################################

suppressPackageStartupMessages({
  # DESeq2 must load BEFORE tidyverse: its dependencies (matrixStats, S4Vectors,
  # IRanges) mask dplyr::count, rename, desc, etc. if attached afterwards.
  library(DESeq2)
  library(tidyverse)
  library(UpSetR)   # install.packages("UpSetR")
})

# ============================== CONFIG ========================================
setwd("/Users/jmlemas/Library/CloudStorage/OneDrive-MichiganStateUniversity/Dissertation/Chapter3/B.hybridum_BARD")

drt_dir <- "03_outputs/04_DESeq2/Drought_results"
pin_dir <- "03_outputs/04_DESeq2/Pinoxaden_results"
out_dir <- "03_outputs/05_overlap"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

iprscan_file <- "01_inputs/BraHy_463_v1.2.tsv"

# iprscan: Brahy.S03G0161600.1.v1.1  ->  DESeq2 id: Brahy.S03G0161600.v1.1
to_gene <- function(x) sub("\\.\\d+(\\.v\\d+\\.\\d+)$", "\\1", x)

# --- Thresholds (applied to shrunk LFC; padj is the Wald test vs 0) ----------
tiers <- list(
  lenient = list(padj = 0.05, lfc = 1),
  strict  = list(padj = 0.01, lfc = 2)
)
primary_tier <- "lenient"   # used for classification, membership table, scatter

# "Clearly not DE" (used to call a gene SPECIFIC to one experiment):
# tested (padj not NA), padj above this, AND |shrunk LFC| below this.
# Genes that are neither DE nor clearly null are labelled "ambiguous".
not_de_padj <- 0.1
not_de_lfc  <- 0.5

# --- Sets ---------------------------------------------------------------------
# type:  response    = treatment effect (DW vs WW, Pinoxaden vs Control)
#        interaction = lineage difference in treatment response, vs Control lineage
#        baseline    = untreated lineage vs untreated Control lineage
# kind:  "rds" = readRDS(path) of a shrunk DESeqResults object
#        "dds" = compute `contrast` ("factor,numerator,denominator") from a dds
#                with results() + lfcShrink(type = "ashr"); cached in out_dir
# Sets whose file is missing are skipped.
sets_cfg <- tribble(
  ~set,        ~experiment, ~type,         ~kind, ~path,                                        ~contrast,
  "Drt_all",   "drought",   "response",    "dds", file.path(drt_dir, "dds_full.rds"),            "trt,DW,WW",
  "Drt_D",     "drought",   "response",    "rds", file.path(drt_dir, "resLFC_trt_drt.rds"),      NA_character_,
  "Drt_H",     "drought",   "response",    "rds", file.path(drt_dir, "resLFC_trt_hrb.rds"),      NA_character_,
  "Pin_all",   "pinoxaden", "response",    "dds", file.path(pin_dir, "dds_full.rds"),            "trt,Pinoxaden,Control",
  "Pin_D",     "pinoxaden", "response",    "rds", file.path(pin_dir, "resLFC_trt_drt.rds"),      NA_character_,
  "Pin_H",     "pinoxaden", "response",    "rds", file.path(pin_dir, "resLFC_trt_hrb.rds"),      NA_character_,
  "DrtInt_D",  "drought",   "interaction", "rds", file.path(drt_dir, "resLFC_interaction_drt.rds"), NA_character_,
  "DrtInt_H",  "drought",   "interaction", "rds", file.path(drt_dir, "resLFC_interaction_hrb.rds"), NA_character_,
  "PinInt_D",  "pinoxaden", "interaction", "rds", file.path(pin_dir, "resLFC_interaction_drt.rds"), NA_character_,
  "PinInt_H",  "pinoxaden", "interaction", "rds", file.path(pin_dir, "resLFC_interaction_hrb.rds"), NA_character_,
  "PinBase_D", "pinoxaden", "baseline",    "rds", file.path(pin_dir, "resLFC_baseline_drt.rds"),    NA_character_,
  "PinBase_H", "pinoxaden", "baseline",    "rds", file.path(pin_dir, "resLFC_baseline_hrb.rds"),    NA_character_
)

# Which sets go into which UpSet figure (drawn for each tier x any/up/down)
panels <- list(
  response        = c("Drt_all", "Drt_D", "Drt_H", "Pin_all", "Pin_D", "Pin_H"),
  lineage_effects = c("DrtInt_D", "DrtInt_H", "PinInt_D", "PinInt_H", "PinBase_D", "PinBase_H")
)

# Drought-vs-pinoxaden pairs for the LFC concordance check
pairs_cfg <- tribble(
  ~drt,       ~pin,
  "Drt_D",    "Pin_D",
  "Drt_H",    "Pin_H",
  "Drt_all",  "Pin_all"
)
# ==============================================================================


# ---- 1. Load full results ----------------------------------------------------
missing <- !file.exists(sets_cfg$path)
if (any(missing)) {
  message("Skipping (file not found): ",
          paste(sets_cfg$set[missing], collapse = ", "))
}
sets_cfg <- sets_cfg[!missing, ]

load_res <- function(set, kind, path, contrast) {
  if (kind == "rds") return(readRDS(path))

  cache <- file.path(out_dir, paste0("cache_", set, ".rds"))   # delete to recompute
  if (file.exists(cache)) return(readRDS(cache))

  message("Computing ", set, " from ", path, " (", contrast, ")")
  dds <- readRDS(path)
  res <- results(dds, contrast = strsplit(contrast, ",")[[1]], alpha = 0.05)
  res <- lfcShrink(dds, res = res, type = "ashr")
  saveRDS(res, cache)
  res
}

res_all <- sets_cfg %>%
  mutate(data = pmap(list(set, kind, path, contrast), function(set, kind, path, contrast) {
    as.data.frame(load_res(set, kind, path, contrast)) %>%
      rownames_to_column("gene_id") %>%
      select(gene_id, log2FoldChange, padj)
  })) %>%
  select(set, experiment, type, data) %>%
  unnest(data)

# Per-set call at the primary tier (drives classification)
tp <- tiers[[primary_tier]]
res_all <- res_all %>%
  mutate(call = case_when(
    is.na(padj)                                              ~ "untested",
    padj < tp$padj & log2FoldChange >  tp$lfc                ~ "up",
    padj < tp$padj & log2FoldChange < -tp$lfc                ~ "down",
    padj > not_de_padj & abs(log2FoldChange) < not_de_lfc    ~ "null",
    TRUE                                                     ~ "ambiguous"
  ))

# DE genes per tier
pass_tier <- function(d, t) {
  d %>%
    filter(!is.na(padj), padj < t$padj, abs(log2FoldChange) > t$lfc) %>%
    mutate(direction = if_else(log2FoldChange > 0, "up", "down"))
}
de_by_tier <- map(tiers, ~ pass_tier(res_all, .x))
de_long    <- de_by_tier[[primary_tier]]

set_sizes <- imap_dfr(de_by_tier, ~ count(.x, set, direction) %>% mutate(tier = .y)) %>%
  pivot_wider(names_from = direction, values_from = n, values_fill = 0) %>%
  mutate(total = rowSums(across(any_of(c("up", "down")))))
print(set_sizes, n = Inf)
write_csv(set_sizes, file.path(out_dir, "set_sizes_by_tier.csv"))


# ---- 2. UpSet plots ----------------------------------------------------------
# up/down panels only count a gene if it moves in that direction, so
# opposite-direction calls are never counted as "shared".
make_upset <- function(d_long, tier, panel_name, keys, dir_filter = c("any", "up", "down")) {
  dir_filter <- match.arg(dir_filter)
  d <- d_long %>% filter(set %in% keys)
  if (dir_filter != "any") d <- filter(d, direction == dir_filter)

  lst <- d %>% group_by(set) %>% summarise(g = list(unique(gene_id)), .groups = "drop")
  lst <- setNames(lst$g, lst$set)
  lst <- lst[intersect(keys, names(lst))]
  if (length(lst) < 2) {
    message("Skipping upset ", tier, "/", panel_name, "/", dir_filter, ": <2 non-empty sets")
    return(invisible(NULL))
  }

  pdf(file.path(out_dir, sprintf("upset_%s_%s_%s.pdf", tier, panel_name, dir_filter)),
      width = 11, height = 6)
  print(upset(fromList(lst), sets = names(lst), keep.order = TRUE,
              nintersects = 40, order.by = "freq",
              mb.ratio = c(0.6, 0.4), text.scale = 1.3,
              mainbar.y.label = sprintf("Intersection size (%s, %s)", tier, dir_filter),
              sets.x.label = "DE genes"))
  dev.off()
}

for (tr in names(de_by_tier)) {
  for (p in names(panels)) {
    for (dr in c("any", "up", "down")) {
      make_upset(de_by_tier[[tr]], tr, p, panels[[p]], dr)
    }
  }
}


# ---- 3. Classification (primary tier) ---------------------------------------
# Per experiment, summarise each gene across that experiment's RESPONSE sets:
#   up / down   = DE in that direction only
#   mixed       = DE up in one response set and down in another
#   null        = clearly not DE (and no DE call anywhere)
#   ambiguous   = neither DE nor clearly null
#   untested    = not in the results, or padj NA (e.g. low counts)
exp_status <- res_all %>%
  filter(type == "response") %>%
  group_by(experiment, gene_id) %>%
  summarise(
    status = case_when(
      any(call == "up") & any(call == "down")                   ~ "mixed",
      any(call == "up")                                         ~ "up",
      any(call == "down")                                       ~ "down",
      any(call == "null") & all(call %in% c("null", "untested")) ~ "null",
      all(call == "untested")                                   ~ "untested",
      TRUE                                                      ~ "ambiguous"
    ),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = experiment, values_from = status)

for (e in c("drought", "pinoxaden")) {
  if (!e %in% names(exp_status)) exp_status[[e]] <- NA_character_
}

# Wide table: log2FC where the gene is DE (primary tier) in each set
set_cols <- intersect(sets_cfg$set, unique(de_long$set))

wide <- de_long %>%
  select(gene_id, set, log2FoldChange) %>%
  pivot_wider(names_from = set, values_from = log2FoldChange) %>%
  select(gene_id, all_of(set_cols))

mat <- as.matrix(wide[set_cols])

# How many sets does each gene also pass at the strict tier?
n_strict <- de_by_tier$strict %>%
  count(gene_id, name = "n_sets_strict")

is_de <- function(s) s %in% c("up", "down", "mixed")

membership <- wide %>%
  mutate(
    n_sets = rowSums(!is.na(mat)),
    sets   = apply(!is.na(mat), 1, function(r) paste(set_cols[r], collapse = "+"))
  ) %>%
  left_join(n_strict, by = "gene_id") %>%
  mutate(n_sets_strict = replace_na(n_sets_strict, 0L)) %>%
  left_join(exp_status, by = "gene_id") %>%
  mutate(
    drought   = replace_na(drought,   "untested"),
    pinoxaden = replace_na(pinoxaden, "untested"),
    class = case_when(
      drought == "up"   & pinoxaden == "up"   ~ "shared_up",
      drought == "down" & pinoxaden == "down" ~ "shared_down",
      is_de(drought) & is_de(pinoxaden)       ~ "discordant",
      is_de(drought) & pinoxaden == "null"      ~ "drought_specific",
      is_de(drought) & pinoxaden == "ambiguous" ~ "drought_ambiguous_in_pinoxaden",
      is_de(drought) & pinoxaden == "untested"  ~ "drought_only_not_tested_in_pinoxaden",
      is_de(pinoxaden) & drought == "null"      ~ "pinoxaden_specific",
      is_de(pinoxaden) & drought == "ambiguous" ~ "pinoxaden_ambiguous_in_drought",
      is_de(pinoxaden) & drought == "untested"  ~ "pinoxaden_only_not_tested_in_drought",
      TRUE                                      ~ "lineage_effects_only"
    )
  ) %>%
  rename(drought_status = drought, pinoxaden_status = pinoxaden)

class_levels <- c("shared_up", "shared_down", "discordant",
                  "drought_specific", "drought_ambiguous_in_pinoxaden",
                  "drought_only_not_tested_in_pinoxaden",
                  "pinoxaden_specific", "pinoxaden_ambiguous_in_drought",
                  "pinoxaden_only_not_tested_in_drought",
                  "lineage_effects_only")


# ---- 4. Threshold-free LFC concordance --------------------------------------
# If drought and pinoxaden share a response, drought-DE genes should move the
# same way in pinoxaden even when they miss significance there.
pairs_cfg <- pairs_cfg %>% filter(drt %in% res_all$set, pin %in% res_all$set)

if (nrow(pairs_cfg) > 0) {
  pair_df <- pmap_dfr(pairs_cfg, function(drt, pin) {
    a <- res_all %>% filter(set == drt, !is.na(padj)) %>%
      select(gene_id, lfc_drt = log2FoldChange, call_drt = call)
    b <- res_all %>% filter(set == pin, !is.na(padj)) %>%
      select(gene_id, lfc_pin = log2FoldChange)
    inner_join(a, b, by = "gene_id") %>% mutate(pair = paste(drt, "vs", pin))
  })

  safe_cor <- function(x, y) {
    if (length(x) >= 3) suppressWarnings(cor(x, y, method = "spearman")) else NA_real_
  }

  concordance <- pair_df %>%
    group_by(pair) %>%
    summarise(
      n_tested_both       = n(),
      rho_all_genes       = safe_cor(lfc_drt, lfc_pin),
      n_drought_DE        = sum(call_drt %in% c("up", "down")),
      rho_drought_DE      = safe_cor(lfc_drt[call_drt %in% c("up", "down")],
                                     lfc_pin[call_drt %in% c("up", "down")]),
      frac_same_sign_DE   = mean(sign(lfc_drt[call_drt %in% c("up", "down")]) ==
                                 sign(lfc_pin[call_drt %in% c("up", "down")])),
      .groups = "drop"
    )
  print(concordance)
  write_csv(concordance, file.path(out_dir, "pair_concordance.csv"))

  p <- ggplot(pair_df, aes(lfc_drt, lfc_pin)) +
    geom_hline(yintercept = 0, colour = "grey70") +
    geom_vline(xintercept = 0, colour = "grey70") +
    geom_point(data = filter(pair_df, !call_drt %in% c("up", "down")),
               colour = "grey70", size = 0.6, alpha = 0.4) +
    geom_point(data = filter(pair_df, call_drt %in% c("up", "down")),
               colour = "firebrick", size = 1.2, alpha = 0.7) +
    facet_wrap(~ pair, scales = "free") +
    labs(x = "Drought log2FC (shrunk)", y = "Pinoxaden log2FC (shrunk)",
         title = "LFC concordance; red = DE in the drought set",
         subtitle = sprintf("Drought DE call: padj < %s & |LFC| > %s", tp$padj, tp$lfc)) +
    theme_bw()
  ggsave(file.path(out_dir, "lfc_concordance_scatter.pdf"), p,
         width = 4 * nrow(pairs_cfg) + 1, height = 4.5)
}


# ---- 5. InterProScan annotation ---------------------------------------------
ipr_names <- c("protein", "md5", "len", "analysis", "sig_acc", "sig_desc", "start", "stop",
               "score", "status", "date", "ipr_acc", "ipr_desc", "go", "pathways")

ipr <- read_tsv(iprscan_file, col_names = ipr_names, col_types = cols(.default = "c"),
                na = c("", "-"), quote = "", progress = FALSE) %>%
  mutate(gene_id = to_gene(protein))

collapse_u <- function(x, sep = "; ") {
  x <- unique(na.omit(x))
  if (length(x) == 0) NA_character_ else paste(sort(x), collapse = sep)
}

annot <- ipr %>%
  group_by(gene_id) %>%      # union across isoforms
  summarise(
    interpro = collapse_u(if_else(is.na(ipr_acc), NA_character_, paste(ipr_acc, ipr_desc))),
    pfam     = collapse_u(sig_desc[analysis == "Pfam"]),
    go_terms = collapse_u(unlist(str_extract_all(go, "GO:\\d{7}"))),
    pathways = collapse_u(pathways),
    .groups = "drop"
  )

# Long GO table for topGO / clusterProfiler (enricher TERM2GENE)
gene2go <- ipr %>%
  filter(!is.na(go)) %>%
  mutate(go = str_extract_all(go, "GO:\\d{7}")) %>%
  unnest(go) %>%
  distinct(go, gene_id)
write_tsv(gene2go, file.path(out_dir, "gene2go.tsv"))

# Proteins with no hits are absent from iprscan output, so <100% is normal;
# ~0% means to_gene() doesn't match your gene ids.
hit <- mean(membership$gene_id %in% annot$gene_id)
message(sprintf("%.1f%% of DE genes have InterProScan annotation", 100 * hit))
if (hit < 0.3) {
  warning("Low annotation match. DE ids look like: ",
          paste(head(membership$gene_id, 3), collapse = ", "),
          " | iprscan ids (after to_gene) look like: ",
          paste(head(unique(ipr$gene_id), 3), collapse = ", "))
}

annotated <- membership %>%
  left_join(annot, by = "gene_id") %>%
  arrange(factor(class, levels = class_levels), desc(n_sets_strict), desc(n_sets))

write_csv(annotated, file.path(out_dir, "DE_membership_annotated.csv"))

annotated %>%
  count(class, name = "n_genes") %>%
  arrange(factor(class, levels = class_levels)) %>%
  write_csv(file.path(out_dir, "class_summary.csv"))

annotated %>%
  group_by(class) %>%
  group_walk(~ write_csv(.x, file.path(out_dir, paste0("genes_", .y$class, ".csv"))))

print(count(annotated, class))
message("Done. Outputs in ", out_dir)