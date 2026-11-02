################################################################################
# topGO enrichment for the drought / pinoxaden overlap analysis
#
# Inputs (all written by DE_overlap_upset_v2.R):
#   DE_membership_annotated.csv  - DE genes with class + per-experiment status
#   gene2go.tsv                  - InterProScan GO mappings (gene_id, go)
#   cache_Drt_all.rds / cache_Pin_all.rds - used only to define the tested universe
#
# For each gene list and ontology (BP/MF/CC):
#   - weight01 + Fisher (primary; accounts for GO hierarchy)
#   - classic Fisher + BH FDR (reference)
# NOTE: weight01 p-values are NOT multiple-testing corrected (terms are not
# independent). Use them to rank and explore; lean on fdr_classic for strict claims.
################################################################################

suppressPackageStartupMessages({
  # Bioconductor packages first, tidyverse LAST so dplyr verbs are not masked
  library(topGO)      # BiocManager::install("topGO")
  library(DESeq2)     # needed to read the cached DESeqResults objects
  library(tidyverse)
})

# ============================== CONFIG ========================================
setwd("/Users/jmlemas/Library/CloudStorage/OneDrive-MichiganStateUniversity/Dissertation/Chapter3/B.hybridum_BARD")

overlap_dir <- "03_outputs/05_overlap"
out_dir     <- "03_outputs/06_topGO"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ontologies <- c("BP", "MF", "CC")
node_size  <- 10     # ignore GO terms with fewer annotated genes than this
min_list   <- 10     # skip lists with fewer ANNOTATED genes than this
plot_p     <- 0.05   # weight01 p cutoff for dotplots
plot_top   <- 10     # max terms per list in dotplots
# ==============================================================================


# ---- 1. Inputs ----------------------------------------------------------------
m <- read_csv(file.path(overlap_dir, "DE_membership_annotated.csv"), show_col_types = FALSE)

gene2go <- read_tsv(file.path(overlap_dir, "gene2go.tsv"), show_col_types = FALSE)
geneID2GO <- split(gene2go$go, gene2go$gene_id)
message(length(geneID2GO), " genes with at least one GO term")

# Tested universe = genes with a non-NA padj in the overall-treatment results
tested <- function(f) {
  r <- as.data.frame(readRDS(f))
  rownames(r)[!is.na(r$padj)]
}
tested_drt <- tested(file.path(overlap_dir, "cache_Drt_all.rds"))
tested_pin <- tested(file.path(overlap_dir, "cache_Pin_all.rds"))

universes <- list(
  drought   = tested_drt,
  pinoxaden = tested_pin,
  both      = intersect(tested_drt, tested_pin)   # genes that COULD be classed specific/shared
)

# ---- 2. Gene lists -------------------------------------------------------------
# universe: the background each list is tested against.
#  - class-based lists (specific / shared) use "both": a gene can only be called
#    specific or shared if it was tested in both experiments.
#  - experiment-level DE lists use that experiment's tested genes.
mk <- function(universe, cond) list(universe = universe, genes = m$gene_id[which(cond)])

gene_lists <- list(
  drought_specific_up   = mk("both",      m$class == "drought_specific" & m$drought_status == "up"),
  drought_specific_down = mk("both",      m$class == "drought_specific" & m$drought_status == "down"),
  shared_up             = mk("both",      m$class == "shared_up"),
  shared_down           = mk("both",      m$class == "shared_down"),
  drought_DE_up         = mk("drought",   m$drought_status == "up"),
  drought_DE_down       = mk("drought",   m$drought_status == "down"),
  pinoxaden_DE_up       = mk("pinoxaden", m$pinoxaden_status == "up"),
  pinoxaden_DE_down     = mk("pinoxaden", m$pinoxaden_status == "down")
)


# ---- 3. Run topGO --------------------------------------------------------------
run_topgo <- function(lst, ont, genes, bg) {
  geneList <- factor(as.integer(bg %in% genes), levels = c(0, 1))
  names(geneList) <- bg

  go_data <- new("topGOdata",
                 description = lst, ontology = ont,
                 allGenes = geneList,
                 annot = annFUN.gene2GO, gene2GO = geneID2GO,
                 nodeSize = node_size)

  res_w <- runTest(go_data, algorithm = "weight01", statistic = "fisher")
  res_c <- runTest(go_data, algorithm = "classic",  statistic = "fisher")

  p_w <- score(res_w)
  p_c <- score(res_c)

  tab <- GenTable(go_data, weight01 = res_w, topNodes = length(p_w), numChar = 200) %>%
    as_tibble() %>%
    dplyr::select(GO.ID, Term, Annotated, Significant, Expected) %>%
    mutate(
      p_weight01  = p_w[GO.ID],
      p_classic   = p_c[GO.ID],
      fdr_classic = p.adjust(p_c, method = "BH")[GO.ID],
      fold        = Significant / pmax(Expected, 1e-9)
    ) %>%
    filter(Significant > 0)

  # Which of the list's genes drive each term
  sig  <- sigGenes(go_data)
  gt   <- genesInTerm(go_data)
  tab$genes <- map_chr(tab$GO.ID, ~ paste(intersect(gt[[.x]], sig), collapse = ";"))

  arrange(tab, p_weight01)
}

all_res <- list()
summ    <- list()

for (lst in names(gene_lists)) {
  gl    <- gene_lists[[lst]]
  bg    <- intersect(universes[[gl$universe]], names(geneID2GO))
  genes <- intersect(gl$genes, bg)

  summ[[lst]] <- tibble(
    list = lst, universe = gl$universe,
    n_genes = length(unique(gl$genes)),
    n_annotated = length(genes),
    n_background_annotated = length(bg)
  )

  if (length(genes) < min_list) {
    message("Skipping ", lst, ": only ", length(genes), " annotated genes (< ", min_list, ")")
    next
  }

  for (ont in ontologies) {
    message("topGO: ", lst, " / ", ont, " (", length(genes), " genes vs ", length(bg), " background)")
    tab <- run_topgo(lst, ont, genes, bg) %>%
      mutate(list = lst, ontology = ont, .before = 1)
    write_csv(tab, file.path(out_dir, sprintf("topGO_%s_%s.csv", lst, ont)))
    all_res[[paste(lst, ont)]] <- tab
  }
}

summ <- bind_rows(summ)
print(summ)
write_csv(summ, file.path(out_dir, "list_summary.csv"))

all_res <- bind_rows(all_res)
if (nrow(all_res) == 0) stop("No lists passed the minimum size; lower `min_list` or check gene lists.")
write_csv(all_res, file.path(out_dir, "topGO_all_results.csv"))


# ---- 4. Dotplots ---------------------------------------------------------------
for (ont in ontologies) {
  d <- all_res %>%
    filter(ontology == ont, p_weight01 < plot_p, Significant >= 2) %>%
    group_by(list) %>%
    slice_min(p_weight01, n = plot_top, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(list, fold) %>%
    mutate(key = factor(paste(list, Term, sep = "___"), levels = unique(paste(list, Term, sep = "___"))))

  if (nrow(d) == 0) {
    message("No terms below p < ", plot_p, " for ", ont, "; skipping plot")
    next
  }

  p <- ggplot(d, aes(fold, key, size = Significant, colour = p_weight01)) +
    geom_point() +
    facet_grid(list ~ ., scales = "free_y", space = "free_y") +
    scale_y_discrete(labels = function(x) str_wrap(sub("^.*___", "", x), 45)) +
    scale_colour_viridis_c(direction = -1, name = "weight01 p") +
    labs(x = "Fold enrichment (observed / expected)", y = NULL,
         title = paste("topGO", ont),
         subtitle = sprintf("Top %d terms per list with weight01 p < %s", plot_top, plot_p)) +
    theme_bw() +
    theme(strip.text.y = element_text(angle = 0))

  ggsave(file.path(out_dir, sprintf("topGO_dotplot_%s.pdf", ont)), p,
         width = 9, height = max(4, 0.35 * nrow(d) + 2), limitsize = FALSE)
}

message("Done. Outputs in ", out_dir)