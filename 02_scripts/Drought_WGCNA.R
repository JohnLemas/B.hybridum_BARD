################################################################################
# Date: 5-14-2026
# Author: John Lemas (Optimized & Debugged using Gemini; further fixed by Claude)
# Description: BFT project WGCNA with Integrated DESeq2 Output
#
# GOAL OF THIS VERSION: identify hub genes whose drought (WW -> DW) response
# is SHARED between the D (drought recurrent selection) and H (herbicide
# recurrent selection) lineages. Network/module detection still uses all
# samples (all lineages x both treatments) so the co-expression structure is
# estimated from the full dataset, as is standard. What changed is which
# DESeq2 contrasts get joined onto the module results and how "shared hub
# gene" is defined.
#
# FIX LOG (see chat for full explanation):
#  - Removed references to `resLFC_line`, `deseq_gen_df`, and `resLFC_trt`,
#    none of which were ever created anywhere in this pipeline. WGCNA.R runs
#    as its own script/session, so any DESeq2 result it needs must be loaded
#    from the .rds files DESeq2.R already saved, not assumed to be in memory.
#  - Reused dds_full.rds only for the network/module construction (as
#    before). For the biological question of interest, this version instead
#    loads the WITHIN-LINEAGE drought-response contrasts DESeq2.R already
#    computed from dds_int: resLFC_trt_drt.rds (D lineage's WW->DW response)
#    and resLFC_trt_hrb.rds (H lineage's WW->DW response). These are the
#    correct objects for "how does each lineage handle drought" -- a plain
#    lineage-vs-control contrast (what the previous fix used) does NOT
#    answer that question, since it ignores treatment entirely.
#  - Added an explicit "shared hub gene" definition: high module membership
#    (MM_OwnModule >= 0.8) AND significant drought response in BOTH D and H
#    lineages, using the same padj < 0.01 & |log2FC| > 2 threshold DESeq2.R
#    already uses for sig_trt_drt / sig_trt_hrb.
################################################################################

# 1. SETUP & PACKAGE LOADING
options(stringsAsFactors = FALSE)

required_packages <- c("DESeq2", "tidyverse", "WGCNA", "ggdendro")
for (package in required_packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    if (package %in% c("DESeq2", "WGCNA")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(package)
    } else {
      install.packages(package)
    }
  }
  library(package, character.only = TRUE)
}

# Set environment
setwd("/Users/jmlemas/Library/CloudStorage/OneDrive-MichiganStateUniversity/Dissertation/Chapter3/BARD_drought/RNA_seq/03_outputs/04_DESeq2/")

wgcna_dir <- "wgcna"

# 2. LOAD PRE-COMPUTED DATA
dds           <- readRDS("dds_full.rds")
coldata       <- as.data.frame(colData(dds))

# Within-lineage drought-response contrasts, already computed & shrunk in
# DESeq2.R from the interaction model (dds_int). These represent each
# lineage's OWN WW -> DW response, which is what "how are lineages dealing
# with treatment" actually requires -- not a lineage-vs-control comparison.
resLFC_trt_drt <- readRDS("resLFC_trt_drt.rds")   # D lineage: drought response
resLFC_trt_hrb <- readRDS("resLFC_trt_hrb.rds")   # H lineage: drought response

# 3. TRANSFORMATION & QUALITY CHECK
dds_vst   <- varianceStabilizingTransformation(dds, blind = FALSE)
geneCount <- t(assay(dds_vst))

gsg <- goodSamplesGenes(geneCount, verbose = 3)
if (!gsg$allOK) {
  geneCount <- geneCount[gsg$goodSamples, gsg$goodGenes]
}

# 4. SAMPLE CLUSTERING DENDROGRAM
sampleTree    <- hclust(dist(geneCount), method = "average")
sampleTree_gg <- dendro_data(sampleTree)

pdf(file = file.path(wgcna_dir, "Tree_DESeq2_means_combined.pdf"),
    width = 18, height = 9, family = "Times", pointsize = 12)
print(
  ggplot() +
    geom_segment(data = segment(sampleTree_gg), aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_text(data = label(sampleTree_gg), aes(x = x, y = y, label = label), hjust = -0.1, size = 3) +
    coord_flip() + scale_y_reverse() + theme_minimal() +
    labs(title = "Hierarchical Clustering Dendrogram")
)
dev.off()

# 5. CHOOSE SOFT-THRESHOLDING POWER
powers <- c(1:10, seq(from = 12, to = 24, by = 2))
disableWGCNAThreads()
sft <- pickSoftThreshold(geneCount, powerVector = powers, verbose = 3)
write.csv(sft$fitIndices, file.path(wgcna_dir, "power_table.csv"), row.names = FALSE)

pdf(file = file.path(wgcna_dir, "soft_power.pdf") ,
    width = 12, height = 5, family = "Times", pointsize = 12)
par(mfrow = c(1,2), cex = 0.9)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit,signed R^2", type="n", main="Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, col="red")
abline(h=0.90, col="red"); abline(h=0.80, col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n", main="Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, col="red")
dev.off()

# 6. BLOCKWISE NETWORK CONSTRUCTION
softPower <- 5
enableWGCNAThreads(nThreads = 12)

# Protect and override namespace safely
tmp_cor <- cor; cor <- WGCNA::cor

net <- blockwiseModules(geneCount, 
                        power = softPower, 
                        corType = "pearson",
                        networkType = "signed",
                        TOMType = "signed", 
                        minModuleSize = 30,
                        maxBlockSize = 20000,
                        reassignThreshold = 0, 
                        mergeCutHeight = 0.25,
                        numericLabels = TRUE, 
                        pamRespectsDendro = FALSE,
                        saveTOMs = TRUE,
                        saveTOMFileBase = "drought_TOM",
                        verbose = 3,
                        nThreads = 12)

cor <- tmp_cor # Restore original cor function
mergedColors <- labels2colors(net$colors)

# Plot Final Network Dendrogram
pdf(file = file.path(wgcna_dir, "dendrogram_final_means.pdf"),
    width = 12, height = 9, family = "Times", pointsize = 12)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

# Save simple tracking module key
tibble(gene_id = names(net$colors), colors = mergedColors) %>% 
  write_delim(file.path(wgcna_dir, "gene_modules.tsv"), delim = "\t")

# 7. MODULE-TRAIT CORRELATION HEATMAP
MEs <- orderMEs(moduleEigengenes(geneCount, colors = mergedColors)$eigengenes)
module_order <- gsub("ME", "", names(MEs))

correlation_metadata <- read.table("../../01_inputs/correlation_metadata.txt", header = TRUE)
# NOTE: correlation_metadata.txt has no sample ID column, so this only works
# if its row order exactly matches BARD_drought_metadata.txt's row order
# (which is what generated the sample order in `dds`/`coldata`/`geneCount`).
# Worth double-checking that assumption still holds for this dataset.

# Add trt (drought treatment) as an explicit numeric trait, since that's the
# variable of primary interest here and correlation_metadata.txt doesn't
# include it. coldata is guaranteed to be in the same sample order as
# geneCount/MEs since it comes directly from `dds`.
correlation_metadata$trt_DW <- as.numeric(coldata$trt == "DW")

moduleTraitcor    <- cor(MEs, correlation_metadata, use = "p")
nSamples          <- nrow(geneCount)
moduleTraitPvalue <- corPvalueStudent(moduleTraitcor, nSamples)

testMatrix <- matrix(paste0(signif(moduleTraitcor, 2), "\n(", signif(moduleTraitPvalue, 1), ")"),
                     nrow = nrow(moduleTraitcor), ncol = ncol(moduleTraitcor), dimnames = dimnames(moduleTraitcor))

pdf(file = file.path(wgcna_dir, "MEheatmap.pdf"),
    width = 12, height = 12, family = "Times", pointsize = 12)
labeledHeatmap(Matrix = moduleTraitcor, xLabels = colnames(correlation_metadata), yLabels = names(MEs), 
               ySymbols = module_order, colorLabels = FALSE, colors = blueWhiteRed(50), 
               textMatrix = testMatrix, setStdMargins = TRUE, yLabelsPosition = "left", cex.text = 1.0, zlim = c(-1, 1))
dev.off()

# 8. EIGENGENE EXPRESSION PLOTS
coldata_augmented <- coldata %>%
  rownames_to_column("Sample") %>%
  mutate(treatment_group = factor(paste(line, trt, sep = "_")))

# Order treatment groups (UT then T)
ordered_groups <- coldata_augmented %>%
  distinct(treatment_group, trt) %>%
  arrange(factor(trt, levels = c("WW", "DW"))) %>%
  pull(treatment_group)

MEs_meta <- MEs %>%
  rownames_to_column("Sample") %>%
  left_join(coldata_augmented, by = "Sample")

# Boxplot of Target Modules
modules_of_interest <- c("MEyellow", "MEgreen", "MEbrown", "MEred") %>% intersect(colnames(MEs_meta))
if(length(modules_of_interest) > 0) {
  pdf(file = file.path(wgcna_dir, "ME_dotplots.pdf"),
      width = 8, height = 6, family = "Times", pointsize = 12)
  print(
    MEs_meta %>%
      dplyr::select(treatment_group, all_of(modules_of_interest)) %>%
      pivot_longer(cols = -treatment_group, names_to = "Modules", values_to = "Eigenvalues") %>%
      mutate(treatment_group = factor(treatment_group, levels = ordered_groups)) %>%
      ggplot(aes(x = treatment_group, y = Eigenvalues, color = treatment_group)) +
      geom_boxplot(outliers = FALSE, show.legend = FALSE) +
      geom_point(size = 1.0, alpha = 0.5, position = position_jitter(width = 0.2), show.legend = FALSE) +
      facet_wrap(~Modules, scales = "free_y") + theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 14)) +
      labs(title = "Eigengene Expression by Treatment Group", x = "Treatment", y = "Eigengene Value")
  )
  dev.off()
}

# 9. ADVANCED METRIC COMPUTATIONS (MM & GS)
gene_ids <- colnames(geneCount)
modNames <- substring(names(MEs), 3)

# Module Membership (MM) across all modules
geneModuleMembership <- as.data.frame(cor(geneCount, MEs, use = "p"))
MMPvalue             <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
colnames(geneModuleMembership) <- paste0("MM.", modNames)
colnames(MMPvalue)             <- paste0("p.MM.", modNames)

# Robust dynamic assignment loop handling "grey" module cases safely
MM_own   <- numeric(length(gene_ids))
p.MM_own <- numeric(length(gene_ids))

for (i in seq_along(gene_ids)) {
  target_mod <- mergedColors[i]
  col_idx    <- match(target_mod, modNames)
  
  if (!is.na(col_idx)) {
    MM_own[i]   <- abs(geneModuleMembership[i, col_idx])
    p.MM_own[i] <- MMPvalue[i, col_idx]
  } else {
    MM_own[i]   <- 0  # Grey module exceptions default gracefully
    p.MM_own[i] <- 1
  }
}

# Gene Significance (GS) vs. EVERY phenotype trait in correlation_metadata.txt
# (final_height, num_of_tillers, dry_shoot_weight, dry_spike_weight,
# dry_seeds_weight, num_of_seeds), so we can see which trait(s) each gene
# actually tracks instead of assuming it's height. trt_DW is excluded here
# since that's a treatment indicator already captured by the DESeq2 columns,
# not a phenotype measurement.
phenotype_traits <- correlation_metadata %>% dplyr::select(-any_of("trt_DW"))

geneTraitSignificance <- as.data.frame(cor(geneCount, phenotype_traits, use = "p"))
GSPvalue              <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
colnames(geneTraitSignificance) <- paste0("GS.", colnames(phenotype_traits))
colnames(GSPvalue)              <- paste0("p.GS.", colnames(phenotype_traits))

# For each gene, identify which single trait it correlates with most strongly
# in absolute value, and that correlation's magnitude -- this is what tells
# you whether "final_height" specifically stands out, or whether some other
# trait (e.g. dry_seeds_weight) is actually the strongest signal for a given
# gene/module.
abs_GS <- abs(geneTraitSignificance)
top_trait_idx   <- max.col(abs_GS, ties.method = "first")
top_trait_name  <- sub("^GS\\.", "", colnames(geneTraitSignificance)[top_trait_idx])
top_trait_GS    <- abs_GS[cbind(seq_len(nrow(abs_GS)), top_trait_idx)]

# 10. GENERATING THE COMBINED TABLE, FOCUSED ON SHARED D/H DROUGHT RESPONSE

wgcna_base <- tibble(
  gene_id = gene_ids,
  WGCNA_module = mergedColors,
  MM_OwnModule = MM_own,
  p.MM_OwnModule = p.MM_own,
  top_phenotype_trait = top_trait_name,
  top_phenotype_GS = top_trait_GS
) %>%
  bind_cols(geneTraitSignificance) %>%
  bind_cols(GSPvalue)

# --- Join the within-lineage drought-response contrasts (loaded in step 2) ---
# resLFC_trt_drt / resLFC_trt_hrb are each lineage's own WW -> DW response,
# already shrunk (ashr) and computed from the interaction model in DESeq2.R.

deseq_trt_drt_df <- as.data.frame(resLFC_trt_drt) %>%
  rownames_to_column("gene_id") %>%
  dplyr::select(gene_id, log2FoldChange, pvalue, padj) %>%
  rename_with(~paste0("DESeq2_D_drought_", .), -gene_id)

deseq_trt_hrb_df <- as.data.frame(resLFC_trt_hrb) %>%
  rownames_to_column("gene_id") %>%
  dplyr::select(gene_id, log2FoldChange, pvalue, padj) %>%
  rename_with(~paste0("DESeq2_H_drought_", .), -gene_id)

# Significance thresholds:
#  - SIG_PADJ / SIG_LFC match sig_trt_drt / sig_trt_hrb in DESeq2.R
#  - MM_HUB_CUT: standard WGCNA "hub gene" cutoff for module membership
#  - GS_PHENOTYPE_CUT: cutoff for "this gene's expression meaningfully
#    tracks a phenotype trait." 0.5 is a reasonable starting point (WGCNA
#    tutorials often use 0.2 as a much looser bar) -- tighten or loosen this
#    depending on how many candidates you want to review.
SIG_PADJ         <- 0.01
SIG_LFC          <- 2
MM_HUB_CUT       <- 0.8
GS_PHENOTYPE_CUT <- 0.5

combined_master_table <- wgcna_base %>%
  left_join(deseq_trt_drt_df, by = "gene_id") %>%
  left_join(deseq_trt_hrb_df, by = "gene_id") %>%
  mutate(
    sig_in_D = !is.na(DESeq2_D_drought_padj) &
               DESeq2_D_drought_padj < SIG_PADJ &
               abs(DESeq2_D_drought_log2FoldChange) > SIG_LFC,
    sig_in_H = !is.na(DESeq2_H_drought_padj) &
               DESeq2_H_drought_padj < SIG_PADJ &
               abs(DESeq2_H_drought_log2FoldChange) > SIG_LFC,
    is_hub_gene            = MM_OwnModule >= MM_HUB_CUT,
    is_shared_drought_hub  = is_hub_gene & sig_in_D & sig_in_H,
    phenotype_relevant     = !is.na(top_phenotype_GS) & top_phenotype_GS >= GS_PHENOTYPE_CUT,
    # The three-way filter: network hub AND drought-responsive in both
    # lineages AND meaningfully tied to a measured phenotype trait.
    is_priority_candidate  = is_shared_drought_hub & phenotype_relevant
  ) %>%
  arrange(desc(is_priority_candidate), desc(is_shared_drought_hub), desc(MM_OwnModule))

# Save full master table
write.csv(combined_master_table, file.path(wgcna_dir, "DESeq2_WGCNA_master_metrics.csv"), row.names = FALSE)

# Save just the genes of primary interest: hub genes whose drought response
# is shared between the D and H lineages
shared_hub_genes <- combined_master_table %>%
  filter(is_shared_drought_hub) %>%
  arrange(WGCNA_module, desc(MM_OwnModule))

write.csv(shared_hub_genes, file.path(wgcna_dir, "shared_D_H_drought_hubgenes.csv"), row.names = FALSE)

# Save the top-priority shortlist: hub + shared D/H drought DE + tied to a
# measured phenotype trait. This is the three-way filter requested.
priority_candidates <- combined_master_table %>%
  filter(is_priority_candidate) %>%
  arrange(WGCNA_module, desc(top_phenotype_GS))

write.csv(priority_candidates, file.path(wgcna_dir, "priority_candidate_hubgenes.csv"), row.names = FALSE)

# Quick summary: which modules are enriched for shared D/H drought hub genes,
# plus which phenotype trait each module's eigengene correlates with most
# strongly (pulled from the module-trait heatmap computed in step 7).
top_module_trait <- apply(abs(moduleTraitcor), 1, function(x) colnames(moduleTraitcor)[which.max(x)])
top_module_trait_r <- sapply(seq_along(top_module_trait), function(i) {
  moduleTraitcor[i, top_module_trait[i]]
})
module_trait_lookup <- tibble(
  WGCNA_module        = module_order,
  top_module_trait    = top_module_trait,
  top_module_trait_r  = round(top_module_trait_r, 3)
)

module_summary <- combined_master_table %>%
  group_by(WGCNA_module) %>%
  summarise(
    n_genes               = n(),
    n_shared_drought_hub  = sum(is_shared_drought_hub),
    n_priority_candidates = sum(is_priority_candidate),
    pct_shared            = round(100 * n_shared_drought_hub / n_genes, 2),
    .groups = "drop"
  ) %>%
  left_join(module_trait_lookup, by = "WGCNA_module") %>%
  arrange(desc(n_priority_candidates), desc(n_shared_drought_hub))

write.csv(module_summary, file.path(wgcna_dir, "module_summary_shared_drought_hubs.csv"), row.names = FALSE)

cat("\nPipeline successful!\n")
cat(" - Full table:            03_outputs/04_DESeq2/wgcna/DESeq2_WGCNA_master_metrics.csv\n")
cat(" - Shared D/H hub genes:  03_outputs/04_DESeq2/wgcna/shared_D_H_drought_hubgenes.csv\n")
cat(" - Priority candidates:   03_outputs/04_DESeq2/wgcna/priority_candidate_hubgenes.csv\n")
cat(" - Module summary:        03_outputs/04_DESeq2/wgcna/module_summary_shared_drought_hubs.csv\n")
cat("Found", nrow(shared_hub_genes), "hub genes with a shared drought response in the D and H lineages.\n")
cat("Of those,", nrow(priority_candidates), "are also strongly tied to a measured phenotype trait (GS >=", GS_PHENOTYPE_CUT, ").\n")
