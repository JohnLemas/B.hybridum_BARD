################################################################################
# DE overlap (UpSet) + drought-specific vs shared-stress classification
# + InterProScan annotation
#
# Reads the sig_*.txt files written by Drought_DESeq2.R and Pinox_DESeq2.R,
# builds UpSet plots, classifies genes, and joins InterProScan annotations.
# Edit the CONFIG block only.
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(UpSetR)   # install.packages("UpSetR")
})

# ============================== CONFIG ========================================
setwd("/Users/jmlemas/Library/CloudStorage/OneDrive-MichiganStateUniversity/Dissertation/Chapter3/B.hybridum_BARD")

drt_dir <- "03_outputs/04_DESeq2/Drought_results"
pin_dir <- "03_outputs/04_DESeq2/Pinoxaden_results"
out_dir <- "03_outputs/05_overlap"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Full InterProScan TSV (the 15-column format, no header)
iprscan_file <- "01_inputs/BraHy_463_v1.2.tsv"   # <-- EDIT

# iprscan reports PROTEIN/TRANSCRIPT ids; DESeq2 rows are GENE ids.
# Map one to the other here. Default strips a trailing ".1", ".2", ...
# Change if your isoform naming differs (e.g. "-T1", ".t1", "_mRNA").
to_gene <- function(x) sub("\\.\\d+(\\.v\\d+\\.\\d+)$", "\\1", x)     # <-- EDIT if needed

# type:
#   response    = a treatment (DW vs WW or Pinoxaden vs Control) effect.
#                 Only these are used to classify genes as drought-specific/shared.
#   interaction = lineage-specific DIFFERENCE in treatment response (vs Control lineage)
#   baseline    = untreated lineage vs Control lineage (constitutive change from selection)
# Files that don't exist are skipped automatically.
sets_cfg <- tribble(
  ~set,         ~experiment,  ~type,          ~file,
  "Drt_all",    "drought",    "response",     file.path(drt_dir, "sig_trt_all.txt"),
  "Drt_D",      "drought",    "response",     file.path(drt_dir, "sig_trt_drt.txt"),
  "Drt_H",      "drought",    "response",     file.path(drt_dir, "sig_trt_hrb.txt"),
  "Pin_all",    "pinoxaden",  "response",     file.path(pin_dir, "sig_trt_genes.txt"),
#  "Pin_D",      "pinoxaden",  "response",     file.path(pin_dir, "sig_trt_drt.txt"),  # optional, see notes
#  "Pin_H",      "pinoxaden",  "response",     file.path(pin_dir, "sig_trt_hrb.txt"),  # optional, see notes
  "DrtInt_D",   "drought",    "interaction",  file.path(drt_dir, "sig_drt_genes.txt"),
  "DrtInt_H",   "drought",    "interaction",  file.path(drt_dir, "sig_herb_genes.txt"),
  "PinInt_D",   "pinoxaden",  "interaction",  file.path(pin_dir, "sig_drt_genes.txt"),
  "PinInt_H",   "pinoxaden",  "interaction",  file.path(pin_dir, "sig_herb_genes.txt"),
  "PinBase_D",  "pinoxaden",  "baseline",     file.path(pin_dir, "sig_baseline_drt.txt"),
  "PinBase_H",  "pinoxaden",  "baseline",     file.path(pin_dir, "sig_baseline_hrb.txt")
)

# Which sets go into which UpSet figure (each is drawn for any / up / down)
panels <- list(
  response        = c("Drt_all", "Drt_D", "Drt_H", "Pin_all", "Pin_D", "Pin_H"),
  lineage_effects = c("DrtInt_D", "DrtInt_H", "PinInt_D", "PinInt_H", "PinBase_D", "PinBase_H")
)
# ==============================================================================


# ---- 1. Read DE gene sets ----------------------------------------------------
read_sig <- function(f) {
  # write.table() wrote a header with one fewer field than the data rows,
  # so read.table() automatically uses column 1 as row names (= gene ids)
  df <- tryCatch(
    read.table(f, header = TRUE, sep = "\t", quote = "", comment.char = "", check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(df) || nrow(df) == 0) {
    return(tibble(gene_id = character(), log2FoldChange = numeric()))
  }
  tibble(gene_id = rownames(df), log2FoldChange = df$log2FoldChange)
}

missing <- !file.exists(sets_cfg$file)
if (any(missing)) message("Skipping (file not found): ", paste(sets_cfg$set[missing], collapse = ", "))
sets_cfg <- sets_cfg[!missing, ]

de_long <- sets_cfg %>%
  mutate(data = map(file, read_sig)) %>%
  unnest(data) %>%
  mutate(direction = if_else(log2FoldChange > 0, "up", "down"))

set_sizes <- de_long %>%
  count(set, direction) %>%
  pivot_wider(names_from = direction, values_from = n, values_fill = 0) %>%
  mutate(total = up + down)
print(set_sizes)
write_csv(set_sizes, file.path(out_dir, "set_sizes.csv"))

sets_cfg <- filter(sets_cfg, set %in% de_long$set)   # drop sets with zero DE genes


# ---- 2. UpSet plots ----------------------------------------------------------
# "up"/"down" panels only count a gene in a set if it moves in that direction,
# so opposite-direction calls are never counted as "shared".
make_upset <- function(panel_name, keys, dir_filter = c("any", "up", "down")) {
  dir_filter <- match.arg(dir_filter)
  d <- de_long %>% filter(set %in% keys)
  if (dir_filter != "any") d <- filter(d, direction == dir_filter)

  lst <- d %>% group_by(set) %>% summarise(g = list(unique(gene_id)), .groups = "drop")
  lst <- setNames(lst$g, lst$set)
  lst <- lst[intersect(keys, names(lst))]            # config order, empty sets dropped
  if (length(lst) < 2) {
    message("Skipping upset ", panel_name, "/", dir_filter, ": <2 non-empty sets")
    return(invisible(NULL))
  }

  pdf(file.path(out_dir, sprintf("upset_%s_%s.pdf", panel_name, dir_filter)), width = 11, height = 6)
  print(upset(fromList(lst), sets = names(lst), keep.order = TRUE,
              nintersects = 40, order.by = "freq",
              mb.ratio = c(0.6, 0.4), text.scale = 1.3,
              mainbar.y.label = paste0("Intersection size (", dir_filter, ")"),
              sets.x.label = "DE genes"))
  dev.off()
}

for (p in names(panels)) for (dr in c("any", "up", "down")) make_upset(p, panels[[p]], dr)


# ---- 3. Wide membership table + classification -------------------------------
set_cols <- sets_cfg$set

wide <- de_long %>%
  select(gene_id, set, log2FoldChange) %>%
  pivot_wider(names_from = set, values_from = log2FoldChange) %>%
  select(gene_id, all_of(set_cols))

mat <- as.matrix(wide[set_cols])

# Was each gene expressed enough to be tested in each experiment?
# (same filter as the DESeq2 scripts: >= 10 counts in >= 3 samples)
tested_genes <- function(cts_file) {
  cts <- read.table(cts_file, header = TRUE, sep = "\t", quote = "", check.names = FALSE)
  rownames(cts)[rowSums(cts >= 10) >= 3]
}
tested_drt <- tested_genes(file.path(drt_dir, "cts.txt"))
tested_pin <- tested_genes(file.path(pin_dir, "cts.txt"))

drt_resp <- sets_cfg$set[sets_cfg$type == "response" & sets_cfg$experiment == "drought"]
pin_resp <- sets_cfg$set[sets_cfg$type == "response" & sets_cfg$experiment == "pinoxaden"]

any_sign <- function(m, cols, test) {
  if (length(cols) == 0) return(rep(FALSE, nrow(m)))
  rowSums(test(m[, cols, drop = FALSE]), na.rm = TRUE) > 0
}
up   <- function(x) x > 0
down <- function(x) x < 0

membership <- wide %>%
  mutate(
    n_sets = rowSums(!is.na(mat)),
    sets   = apply(!is.na(mat), 1, function(r) paste(set_cols[r], collapse = "+")),
    drt_up = any_sign(mat, drt_resp, up),  drt_dn = any_sign(mat, drt_resp, down),
    pin_up = any_sign(mat, pin_resp, up),  pin_dn = any_sign(mat, pin_resp, down),
    in_drt = drt_up | drt_dn,
    in_pin = pin_up | pin_dn,
    tested_drought  = gene_id %in% tested_drt,
    tested_pinoxaden = gene_id %in% tested_pin,
    class = case_when(
      in_drt & in_pin & !drt_dn & !pin_dn ~ "shared_up",
      in_drt & in_pin & !drt_up & !pin_up ~ "shared_down",
      in_drt & in_pin                     ~ "discordant",
      in_drt & tested_pinoxaden           ~ "drought_specific",
      in_drt                              ~ "drought_only_not_tested_in_pinoxaden",
      in_pin & tested_drought             ~ "pinoxaden_specific",
      in_pin                              ~ "pinoxaden_only_not_tested_in_drought",
      TRUE                                ~ "lineage_effects_only"   # only interaction/baseline sets
    )
  ) %>%
  select(-drt_up, -drt_dn, -pin_up, -pin_dn, -in_drt, -in_pin)


# ---- 4. InterProScan annotation ----------------------------------------------
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

# Long GO table, ready for topGO / clusterProfiler (enricher TERM2GENE)
gene2go <- ipr %>%
  filter(!is.na(go)) %>%
  mutate(go = str_extract_all(go, "GO:\\d{7}")) %>%
  unnest(go) %>%
  distinct(go, gene_id)
write_tsv(gene2go, file.path(out_dir, "gene2go.tsv"))

# Sanity check on the id mapping. Proteins with no hits are absent from iprscan
# output, so <100% is normal; ~0% means to_gene() doesn't match your gene ids.
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
  arrange(factor(class, levels = c("shared_up", "shared_down", "discordant",
                                   "drought_specific", "pinoxaden_specific")),
          desc(n_sets))

write_csv(annotated, file.path(out_dir, "DE_membership_annotated.csv"))

annotated %>%
  count(class, name = "n_genes") %>%
  write_csv(file.path(out_dir, "class_summary.csv"))

# One file per class for quick browsing
annotated %>%
  group_by(class) %>%
  group_walk(~ write_csv(.x, file.path(out_dir, paste0("genes_", .y$class, ".csv"))))

print(count(annotated, class))
message("Done. Outputs in ", out_dir)