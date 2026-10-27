#######################################

#__Date:__ 5-14-2026
#__Author:__ John Lemas
# List of packages required
packages <- c("DESeq2",
              "tidyverse")

# if not installed, install and load it
for (package in packages) {
  if (!(package %in% rownames(installed.packages()))) {
    BiocManager::install(package)
  }
  library(package, character.only = TRUE)
}

#################################################

# Set the working directory:
message("setting working directory")
setwd("/Users/jmlemas/Library/CloudStorage/OneDrive-MichiganStateUniversity/Dissertation/Chapter3/B.hybridum_BARD")
pin_out <- "03_outputs/04_DESeq2/Pinoxaden_results"
if (!dir.exists(pin_out)) {
  dir.create(pin_out, recursive = TRUE)
}

# ===========>------- Full dataset: design ~ gen + trt --------<=========== #

# Import the counts data file
message("importing counts from featureCounts")
countsData <- read.table(file = "03_outputs/03_featureCounts/counts_multi.txt", header = FALSE, row.names = 1, skip = 2)

# Read in the metadata
message("importing metadata")
metadata_full <- read.table(file = "01_inputs/B.hyb_metadata.txt", header = FALSE)

# Add headers
colnames(metadata_full) <- c("fasta1","name1")

# Now add headers for the counts file. the final columns should be the samples in order from the metadata
as.vector(metadata_full$name1)
colnames(countsData) <- c("chr", "start", "stop", "strand", "length", as.vector(metadata_full$name1))

################### COUNT MATRIX INPUT ###################
# Subset the countsData for just the samples
# The pinoxaden study has 33 samples, and they should be
# the first 33 to feed into featureCounts.

message("generateing cts and coldata")
message("subset full counts table to cts")

head(countsData[,6:38])
dim(countsData[,6:38])

# Ensure that this will capture all of the data!
# We only want the sample names and the counts for each gene
# Save just the subset as an object called cts:
cts <- as.matrix(countsData[,6:38])

# Used -M --fraction in featurecounts to deal with multimapping, 
# so wee need to round the counts data to whole numbers
cts <- round(cts)

message("done. generate coldata from pinoxaden experiment metadata..")

# Row names and col names need to be identical between cts and metadata
metadata <- read.table(file = "01_inputs/pinoxaden_metadata.txt", header = FALSE)
colnames(metadata) <- c("fasta1", "line", "trt", "name1", "name2", "rep")

rownames(metadata) <- metadata$name1

# create DESeq2 input and designate factors
# these variables can be used in plots later
coldata <- metadata[,c("line","trt","rep")]
coldata$line <- as.factor(coldata$line)
coldata$trt <- as.factor(coldata$trt)
coldata$rep <- as.factor(coldata$rep)
str(coldata)

message("coldata built. checking..")

# Double check that row names and column names are identical:
rownames(coldata)
colnames(cts)

# double check naming before building dds
all(rownames(coldata) == colnames(cts))

message("saving cts and coldata")

# save these tables as resources in the GitHub repo
write.table(coldata, file = file.path(pin_out, "coldata.txt"), 
            quote = FALSE, sep = "\t", col.names = FALSE)

write.table(cts, file = file.path(pin_out, "cts.txt"), 
            quote = FALSE, sep = "\t")

# Building dds objects for full dataset and for interaction design
# Notes on Design: additive effects look at first item while 
# controlling for the effects of the second item, and vise versa
# If you want to pull specific contrasts from different treatment 
# groups (ie 'results()'), make sure they're in the design

message("building dds object for full treatment effect")
message("design = ~ line + trt")

dds_full <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ line + trt)

message("building dds object with interaction design")
message("design = ~ line + trt + line:trt")

dds_int <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ line + trt + line:trt)

# Alter the below factor levels based on experiment design
# These should be the values in coldata

message("defining factor levels for both dds objects")

dds_full$line <- factor(dds_full$line, levels = c("Control","Drought","Herbicide"))
dds_full$trt  <- factor(dds_full$trt,  levels = c("Control","Pinoxaden"))

dds_int$line <- factor(dds_int$line, levels = c("Control","Drought","Herbicide"))
dds_int$trt  <- factor(dds_int$trt,  levels = c("Control","Pinoxaden"))

# this code is selecting the genes that have normalized counts
# of more than 10 in at least three samples (replicate groups)

message("subsetting for counts greater than 10 for at least three reps")

keep <- rowSums(counts(dds_full) >= 10) >= 3
dds_full <- dds_full[keep, ]

keep <- rowSums(counts(dds_int) >= 10) >= 3
dds_int <- dds_int[keep, ]

# PERFORM DESEQ2 analysis:

message("preforming full analysis")
dds_full <- DESeq(dds_full)

# Save the main dds object
saveRDS(dds_full, file = file.path(pin_out, "dds_full.rds"))

message("preforming interaction analysis")
dds_int <- DESeq(dds_int)

# Save the interaction dds object
saveRDS(dds_int, file = file.path(pin_out, "dds_int.rds"))

# ===========>------- Generating Results --------<=========== #

# in contrast, the last item is the baseline
# this won't work if the contrast you're running wasn't in the initial design
# Main questions:
# 1. What genes are differentially expressed between the two treatments in all samples?
# 2. What genes are differentially expressed between the two treatments in each line?
# 3. Are there differences in expression between the untreated lines?

message("building results for treatment contrast using full dds object")
message('contrast=c("trt", "Pinoxaden", "Control")')
message("lfcThreshold: 2, alpha: 0.01")

res_trt <- results(dds_full,
                   lfcThreshold = 2,
                   alpha = 0.01,
                   contrast=c("trt", "Pinoxaden", "Control"))

summary(res_trt)

### Significant genes for line contrast
sig_trt_genes <- as.data.frame(res_trt)
sig_trt_genes <- sig_trt_genes[which(sig_trt_genes$padj < 0.01 & abs(sig_trt_genes$log2FoldChange) > 2), ]
sig_trt_genes <- sig_trt_genes[order(sig_trt_genes$padj), ]

# View them
print(sig_trt_genes)

# save them
write.table(sig_trt_genes, file = file.path(pin_out, "sig_trt_genes.txt"), 
            quote = FALSE, sep = "\t")

# Save the full results object
saveRDS(res_trt, file = file.path(pin_out, "res_trt.rds"))

# --------<======= Drought Stressed Lines ==========>---------- #

message("generating results for the interaction dds object")
message('name = "lineDrought.trtPinoxaden"')
message("alpha: 0.05")

res_interaction_drt <- results(
  dds_int, 
  name = "lineDrought.trtPinoxaden", 
  alpha = 0.05
)

summary(res_interaction_drt)

# Log fold shrinkage (more robust results)
# The above results are not following the strict H0 of lfc > |±2|
# Using LFC shrinkage will give more robust results, especially for low count genes
# Following up with the post-hoc thresholds below will give a more accurate list of significant genes
# that may have been missed in the strict H0 test above. This is especially important for TE expression

message("preforming log fold shrinkage")

resLFC_interaction_drt <- lfcShrink(dds_int,
                              res = res_interaction_drt,
                              type = "ashr")

summary(resLFC_interaction_drt)

### Significant genes for drought lineage contrast
message("generating significant genes for the treated drought lineage")

sig_drt_genes <- as.data.frame(resLFC_interaction_drt)
sig_drt_genes <- sig_drt_genes[which(sig_drt_genes$padj < 0.01 &
                                      abs(sig_drt_genes$log2FoldChange) > 2), ]
sig_drt_genes <- sig_drt_genes[order(sig_drt_genes$padj), ]

# View them
print(sig_drt_genes)

# save them
write.table(sig_drt_genes, file = file.path(pin_out, "sig_drt_genes.txt"), 
            quote = FALSE, sep = "\t")

# -------<======= Herbicide Stressed Lines ========>-------- #

message("generating results for the treated herbicide lineage")
message('name = "lineHerbicide.trtPinoxaden"')
message("alpha: 0.05")

res_interaction_hrb <- results(
  dds_int,
  name = "lineHerbicide.trtPinoxaden",
  alpha = 0.05
)

summary(res_interaction_hrb)

resLFC_interaction_hrb <- lfcShrink(dds_int,
                              res = res_interaction_hrb,
                              type = "ashr")

summary(resLFC_interaction_hrb)

### Significant genes for herbicide contrast
sig_herb_genes <- as.data.frame(resLFC_interaction_hrb)
sig_herb_genes <- sig_herb_genes[which(sig_herb_genes$padj < 0.01 &
                                              abs(sig_herb_genes$log2FoldChange) > 2), ]
sig_herb_genes <- sig_herb_genes[order(sig_herb_genes$padj), ]

# View them
print(sig_herb_genes)

# save them
write.table(sig_herb_genes, file = file.path(pin_out, "sig_herb_genes.txt"), 
            quote = FALSE, sep = "\t")

# -------<======= Differences between lines in untreated samples ========>-------- #

message("generating baseline expression results for drought lineage")
message('contrast = c("line","Drought","Control")')
message("alpha: 0.05")

res_drt_baseline <- results(
  dds_int,
  contrast = c("line","Drought","Control"),
  alpha = 0.05
)

summary(res_drt_baseline)

resLFC_baseline_drt <- lfcShrink(dds_int,
                              res = res_drt_baseline,
                              type = "ashr")

summary(resLFC_baseline_drt)

### Significant genes for drought lineage contrast
sig_baseline_drt <- as.data.frame(resLFC_baseline_drt)
sig_baseline_drt <- sig_baseline_drt[which(sig_baseline_drt$padj < 0.01 &
                                      abs(sig_baseline_drt$log2FoldChange) > 2), ]
sig_baseline_drt <- sig_baseline_drt[order(sig_baseline_drt$padj), ]

# View them
print(sig_baseline_drt)

# save them
write.table(sig_baseline_drt, file = file.path(pin_out, "sig_baseline_drt.txt"), 
            quote = FALSE, sep = "\t")

message("generating baseline expression results for herbicide lineage")
message('contrast = c("line","Herbicide","Control")')
message("alpha: 0.05")

res_hrb_baseline <- results(
  dds_int,
  contrast = c("line","Herbicide","Control"),
  alpha = 0.05
)

summary(res_hrb_baseline)

resLFC_baseline_hrb <- lfcShrink(dds_int,
                              res = res_hrb_baseline,
                              type = "ashr")

summary(resLFC_baseline_hrb)

### Significant genes for herbicide lineage contrast
sig_baseline_hrb <- as.data.frame(resLFC_baseline_hrb)
sig_baseline_hrb <- sig_baseline_hrb[which(sig_baseline_hrb$padj < 0.01 &
                                      abs(sig_baseline_hrb$log2FoldChange) > 2), ]
sig_baseline_hrb <- sig_baseline_hrb[order(sig_baseline_hrb$padj), ]

# View them
print(sig_baseline_hrb)

# save them
write.table(sig_baseline_hrb, file = file.path(pin_out, "sig_baseline_hrb.txt"), 
            quote = FALSE, sep = "\t")


# Save the subset results objects
saveRDS(resLFC_interaction_hrb, file = file.path(pin_out, "resLFC_interaction_hrb.rds"))
saveRDS(resLFC_interaction_drt, file = file.path(pin_out, "resLFC_interaction_drt.rds"))
saveRDS(resLFC_baseline_drt, file = file.path(pin_out, "resLFC_baseline_drt.rds"))
saveRDS(resLFC_baseline_hrb, file = file.path(pin_out, "resLFC_baseline_hrb.rds"))

library(tidyverse)

# 1. Define a helper function to extract key columns and add custom prefixes
format_contrast <- function(res_obj, prefix, padj_thresh = 0.01, lfc_thresh = 2) {
  df <- as.data.frame(res_obj) %>%
    rownames_to_column(var = "gene_id") %>%
    select(gene_id, log2FoldChange, pvalue, padj) %>%
    mutate(
      is_sig = !is.na(padj) & padj < padj_thresh & abs(log2FoldChange) > lfc_thresh
    )
  
  # Rename columns with the contrast prefix
  colnames(df)[-1] <- paste(prefix, colnames(df)[-1], sep = "_")
  return(df)
}

# 2. Process each of your contrasts
c1 <- format_contrast(res_trt, "trt_main")
c2 <- format_contrast(resLFC_interaction_drt, "int_drt")
c3 <- format_contrast(resLFC_interaction_hrb, "int_hrb")
c4 <- format_contrast(resLFC_baseline_drt, "base_drt")
c5 <- format_contrast(resLFC_baseline_hrb, "base_hrb")

# 3. Merge all contrasts into a single master table
master_table <- list(c1, c2, c3, c4, c5) %>%
  reduce(full_join, by = "gene_id")

# 4. Add a summary column for total contrasts passed per gene
sig_cols <- grep("_is_sig$", colnames(master_table), value = TRUE)
master_table$total_sig_contrasts <- rowSums(master_table[, sig_cols], na.rm = TRUE)

# 5. Export the master matrix
write.table(master_table, 
            file = file.path(pin_out, "master_results_all_contrasts.csv"), 
            quote = FALSE, 
            sep = ",", 
            row.names = FALSE)
