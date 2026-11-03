################################################################################
# Heatmap of shared-stress genes (supported_both shared_up / shared_down)
# across the samples of BOTH experiments.
#
# - Expression = variance-stabilized counts (vst, blind = FALSE) from each
#   experiment's dds_full.rds, z-scored PER GENE WITHIN EACH EXPERIMENT. The two
#   experiments are drawn as side-by-side panels with the same gene order; they
#   are not merged into one matrix because they are separate experiments.
# - Columns are ordered by group (not clustered) so treatment patterns are visible.
# - Rows are split into shared_up / shared_down and clustered within each block.
################################################################################

suppressPackageStartupMessages({
  library(ComplexHeatmap)   # BiocManager::install("ComplexHeatmap")
  library(DESeq2)
  library(tidyverse)        # last, so dplyr verbs are not masked
})

message("shared_heatmap.R v3: unified lineage/treatment labels, blank column titles")

# ============================== CONFIG ========================================
setwd("/Users/jmlemas/Library/CloudStorage/OneDrive-MichiganStateUniversity/Dissertation/Chapter3/B.hybridum_BARD")

robust_dir <- "03_outputs/07_robustness"
drt_dir    <- "03_outputs/04_DESeq2/Drought_results"
pin_dir    <- "03_outputs/04_DESeq2/Pinoxaden_results"
out_dir    <- "03_outputs/08_heatmap"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

classes_to_plot <- c("shared_up", "shared_down")
clip        <- 2       # z-scores clipped to +/- this value
label_genes <- FALSE   # TRUE prints gene ids on the right (readable up to ~100 genes)

# dds colData columns to use as annotation / column order (group 1, then group 2)
experiments <- list(
  drought   = list(dds = file.path(drt_dir, "dds_full.rds"),
                   group_cols = c("lineage", "trt"), title = "Drought"),
  pinoxaden = list(dds = file.path(pin_dir, "dds_full.rds"),
                   group_cols = c("line", "trt"),    title = "Pinoxaden")
)

# The same inbred lines appear in both experiments (C/D/H = Control/Drought/
# Herbicide), and "untreated" is well-watered (drought experiment) or Control
# (pinoxaden experiment). Map both experiments onto ONE set of labels/colours.
# Left side = values found in colData; right side = label shown in the plot.
lineage_map   <- c(C = "Control", Control = "Control",
                   D = "Drought", Drought = "Drought",
                   H = "Herbicide", Herbicide = "Herbicide")
treatment_map <- c(WW = "Untreated / well-watered", Control = "Untreated / well-watered",
                   DW = "Drought", Pinoxaden = "Pinoxaden")

# Lineage and treatment use different colour families so the two annotation
# rows are not confused; the three treatment colours are strongly different.
lineage_cols   <- c(Control = "#1B9E77", Drought = "#7570B3", Herbicide = "#E7298A")
treatment_cols <- c("Untreated / well-watered" = "#FFD92F",
                    Drought = "#1F78B4", Pinoxaden = "#E6550D")
# ==============================================================================


# ---- 1. Genes ------------------------------------------------------------------
m <- read_csv(file.path(robust_dir, "DE_membership_annotated_supported.csv"),
              show_col_types = FALSE, guess_max = Inf) %>%
  filter(support == "supported_both", class %in% classes_to_plot)
message(nrow(m), " supported_both genes: ",
        paste(names(table(m$class)), table(m$class), collapse = ", "))

# ---- 2. VST + z-score per experiment ---------------------------------------------
get_z <- function(ex) {
  dds <- readRDS(ex$dds)
  cd  <- as.data.frame(colData(dds))
  if (!all(ex$group_cols %in% colnames(cd))) {
    stop("colData of ", ex$dds, " has columns: ", paste(colnames(cd), collapse = ", "),
         "\nSet `group_cols` for ", ex$title, " in CONFIG.")
  }
  mat <- assay(vst(dds, blind = FALSE))
  mat <- mat[intersect(m$gene_id, rownames(mat)), , drop = FALSE]

  z <- t(scale(t(mat)))
  z <- z[apply(z, 1, function(r) all(is.finite(r))), , drop = FALSE]   # drop zero-variance rows
  z[z >  clip] <-  clip
  z[z < -clip] <- -clip

  grp <- cd[colnames(z), ex$group_cols, drop = FALSE]
  grp[] <- lapply(grp, as.character)
  colnames(grp) <- c("lineage", "treatment")

  bad_l <- setdiff(unique(grp$lineage), names(lineage_map))
  bad_t <- setdiff(unique(grp$treatment), names(treatment_map))
  if (length(bad_l) || length(bad_t)) {
    stop(ex$title, ": unmapped values. lineage: ", paste(bad_l, collapse = ", "),
         " | treatment: ", paste(bad_t, collapse = ", "),
         "\nAdd them to `lineage_map` / `treatment_map` in CONFIG.")
  }
  grp$lineage   <- factor(unname(lineage_map[grp$lineage]),     levels = names(lineage_cols))
  grp$treatment <- factor(unname(treatment_map[grp$treatment]), levels = names(treatment_cols))

  list(z = z, grp = grp)
}

obj <- lapply(experiments, get_z)

# Check: both experiments should show the same lineage labels (Control/Drought/Herbicide)
for (nm in names(obj)) {
  message("Samples per lineage x treatment, ", nm, ":")
  print(table(obj[[nm]]$grp$lineage, obj[[nm]]$grp$treatment))
}

# Same genes, same order, in both panels
common <- Reduce(intersect, lapply(obj, function(o) rownames(o$z)))
m <- m %>%
  filter(gene_id %in% common) %>%
  mutate(class = factor(class, levels = classes_to_plot)) %>%
  arrange(class)
common <- m$gene_id
row_cls <- m$class
message(length(common), " genes plotted")

col_fun <- circlize::colorRamp2(c(-clip, 0, clip), c("#2166AC", "white", "#B2182B"))

# ---- 3. Heatmaps ------------------------------------------------------------------
make_ht <- function(o, title, main = FALSE) {
  z   <- o$z[common, , drop = FALSE]
  grp <- o$grp
  ord <- do.call(order, unname(as.list(grp)))     # group 1, then group 2
  z   <- z[, ord, drop = FALSE]
  grp <- grp[ord, , drop = FALSE]

  lab <- paste(grp[[1]], grp[[2]], sep = " / ")
  col_split <- factor(lab, levels = unique(lab))

  ann_df <- grp
  colnames(ann_df) <- c("Lineage", "Treatment")
  ann <- HeatmapAnnotation(
    df = ann_df,
    col = list(Lineage = lineage_cols, Treatment = treatment_cols),
    annotation_name_side = "left",
    show_annotation_name = main,   # names and legends drawn once (first panel)
    show_legend = main
  )

  Heatmap(
    z, name = paste0("z_", title), col = col_fun,
    top_annotation = ann,
    column_split = col_split, cluster_columns = FALSE, show_column_names = FALSE,
    column_title = rep("", nlevels(col_split)),   # blank slice titles (long labels overlapped)
    row_split = if (main) row_cls else NULL,
    cluster_rows = main, cluster_row_slices = FALSE,
    show_row_names = (!main) && label_genes, row_names_gp = grid::gpar(fontsize = 5),
    show_heatmap_legend = main,
    heatmap_legend_param = list(title = "Row z-score\n(per experiment)"),
    row_title_rot = 0,
    use_raster = FALSE
  )
}

ht <- make_ht(obj$drought, "Drought", main = TRUE) +
      make_ht(obj$pinoxaden, "Pinoxaden")

pdf(file.path(out_dir, "shared_genes_heatmap.pdf"), width = 15, height = 9)
ht_drawn <- draw(ht, merge_legends = TRUE)
dev.off()

# Gene order as drawn, with annotations, for reference
ro <- unlist(row_order(ht_drawn))
gene_order <- common[ro]
m %>%
  arrange(match(gene_id, gene_order)) %>%
  dplyr::select(gene_id, class, any_of(c("sets", "n_sets_strict", "interpro", "pfam", "go_terms"))) %>%
  write_csv(file.path(out_dir, "heatmap_gene_order.csv"))

message("Done. Outputs in ", out_dir)