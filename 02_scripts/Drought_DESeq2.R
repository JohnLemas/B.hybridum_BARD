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
  suppressPackageStartupMessages(library(package, character.only = TRUE))
}

#################################################

# Set the working directory:
message(">>> setting working directory")
setwd("/Users/jmlemas/Library/CloudStorage/OneDrive-MichiganStateUniversity/Dissertation/Chapter3/B.hybridum_BARD")
drt_out <- "03_outputs/04_DESeq2/Drought_results"
if (!dir.exists(drt_out)) {
  dir.create(drt_out, recursive = TRUE)
}

# ===========>------- Full dataset: design ~ gen + trt --------<=========== #

# Import the counts data file
message(">>> importing counts from featureCounts")
countsData <- read.table(file = "03_outputs/03_featureCounts/counts_multi.txt", header = FALSE, row.names = 1, skip = 2)

# Read in the metadata
message(">>> importing metadata")
metadata_full <- read.table(file = "01_inputs/B.hyb_metadata.txt", header = FALSE)

# Add headers
colnames(metadata_full) <- c("fasta1","name1")

# Now add headers for the counts file. the final columns should be the samples in order from the metadata
as.vector(metadata_full$name1)
colnames(countsData) <- c("chr", "start", "stop", "strand", "length", as.vector(metadata_full$name1))

################### COUNT MATRIX INPUT ###################
# Subset the countsData for just the samples
# The drought study has 144 samples, and they should be
# the last 144 to feed into featureCounts.

message(">>> generating cts and coldata")
message(">>> subsetting full counts table to cts")

head(countsData[,39:182])
dim(countsData[,39:182])

# Ensure that this will capture all of the data!
# We only want the sample names and the counts for each gene
# Save just the subset as an object called cts:
cts <- as.matrix(countsData[,39:182])

# Used -M --fraction in featurecounts to deal with multimapping, 
# so wee need to round the counts data to whole numbers
cts <- round(cts)

message(">>> done. generate coldata from pinoxaden experiment metadata..")

# Row names and col names need to be identical between cts and metadata
metadata <- read.table(file = "01_inputs/drought_metadata.txt", header = TRUE)

rownames(metadata) <- metadata$name2

# create DESeq2 input and designate factors
# these variables can be used in plots later
coldata <- metadata[,c("lineage","trt","rep", "date","batch","line","plate")]
coldata$line <- as.factor(coldata$line)
coldata$trt <- as.factor(coldata$trt)
coldata$rep <- as.factor(coldata$rep)
coldata$lineage <- as.factor(coldata$lineage)
coldata$batch <- as.factor(coldata$batch)
coldata$date <- as.factor(coldata$date)
coldata$plate <- as.factor(coldata$plate)
str(coldata)

message(">>> done. double check row names and col names..")

# Double check that row names and column names are identical:
rownames(coldata)
colnames(cts)

# double check naming before building dds
all(rownames(coldata) == colnames(cts))

message(">>> saving cts and coldata")

# save these tables as resources in the GitHub repo
write.table(coldata, file.path(drt_out, "coldata.txt"),
            quote = FALSE, sep = "\t", col.names = FALSE)

write.table(cts, file = file.path(drt_out, "cts.txt"),
            quote = FALSE, sep = "\t")

# Building dds objects for full dataset and for interaction design
# Notes on Design: additive effects look at first item while controlling for the effects of the second item, and vise versa
# If you want to pull specific contrasts from different treatment groups (ie 'results()'), make sure they're in the design

message(">>> building dds object for full treatment effect")
message("design = ~ lineage + trt")

dds_full <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ lineage + trt)

message(">>>  building dds object with interaction design")
message("design = ~ lineage + trt + lineage:trt")

dds_int <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ lineage + trt + lineage:trt)

# Alter the below factor levels based on experiment design
# These should be the values in coldata

message(">>> defining factor levels for both dds objects")

dds_full$lineage <- factor(dds_full$lineage, levels = c("C","D","H"))
dds_full$trt  <- factor(dds_full$trt,  levels = c("WW","DW"))

dds_int$lineage <- factor(dds_int$lineage, levels = c("C","D","H"))
dds_int$trt  <- factor(dds_int$trt,  levels = c("WW","DW"))

# this code is selecting the genes that have normalized counts of more than 10 in at least three samples (replicate groups)

message(">>> subsetting for counts greater than 10 for at least three reps")

keep <- rowSums(counts(dds_full) >= 10) >= 3
dds_full <- dds_full[keep, ]

keep <- rowSums(counts(dds_int) >= 10) >= 3
dds_int <- dds_int[keep, ]

# PERFORM DESEQ2 analysis:

message(">>> performing full analysis")
dds_full <- DESeq(dds_full)

# Save the main dds object
saveRDS(dds_full, file = file.path(drt_out, "dds_full.rds"))

message(">>> performing interaction analysis")
dds_int <- DESeq(dds_int)

# Save the interaction dds object
saveRDS(dds_int, file = file.path(drt_out, "dds_int.rds"))

# ===========>------- Generating Results --------<=========== #

# in contrast, the last item is the baseline
# this won't work if the contrast you're running wasn't in the initial design
# Main questions:
# 1. What genes are differentially expressed under drought stress in all samples?
# 2. What genes are differentially expressed between control and treated lineages?
# 3. Are there differences in expression between the untreated lineages?

message(">>> building results for treatment contrast using full dds object")
message('contrast=c("trt", "DW", "WW")')
message("lfcThreshold: 2, alpha: 0.01")

res_trt <- results(dds_full,
                   lfcThreshold = 2,
                   alpha = 0.01,
                   contrast=c("trt", "DW", "WW"))

summary(res_trt)

### Significant genes for line contrast
sig_trt_all <- as.data.frame(res_trt)
sig_trt_all <- sig_trt_all[which(sig_trt_all$padj < 0.01 & abs(sig_trt_all$log2FoldChange) > 2), ]
sig_trt_all <- sig_trt_all[order(sig_trt_all$padj), ]

# View them
print(sig_trt_all)

# save them
write.table(sig_trt_all, file = file.path(drt_out, "sig_trt_all.txt"), 
            quote = FALSE, sep = "\t")

# Save the full results object
saveRDS(res_trt, file = file.path(drt_out, "res_trt.rds"))

# --------<======= Drought Stressed Lines ==========>---------- #

message(">>> generating results for the interaction dds object")
message('name = "lineageD.trtDW"')
message("alpha: 0.05")

res_interaction_drt <- results(
  dds_int,
  name = "lineageD.trtDW",
  alpha = 0.05
)

summary(res_interaction_drt)

# Log fold shrinkage (more robust results)
# The above results are not following the strict H0 of lfc > |±2|
# Using LFC shrinkage will give more robust results, especially for low count genes
# Following up with the post-hoc thresholds below will give a more accurate list of significant genes
# that may have been missed in the strict H0 test above. This is especially important for TE expression

message(">>> performing log fold shrinkage")

resLFC_interaction_drt <- lfcShrink(dds_int,
                              res = res_interaction_drt,
                              type = "ashr")

summary(resLFC_interaction_drt)

### Significant genes for drought lineage contrast
message(">>> generating significant genes for the treated drought lineage")

sig_drt_genes <- as.data.frame(resLFC_interaction_drt)
sig_drt_genes <- sig_drt_genes[which(sig_drt_genes$padj < 0.01 &
                                      abs(sig_drt_genes$log2FoldChange) > 2), ]
sig_drt_genes <- sig_drt_genes[order(sig_drt_genes$padj), ]

# View them
print(sig_drt_genes)

# save them
write.table(sig_drt_genes, file = file.path(drt_out, "sig_drt_genes.txt"), 
            quote = FALSE, sep = "\t")

# -------<======= Herbicide Stressed Lines ========>-------- #

message(">>> generating results for the treated herbicide lineage")
message('name = "lineageH.trtDW"')
message("alpha: 0.05")

res_interaction_hrb <- results(
  dds_int,
  name = "lineageH.trtDW",
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
write.table(sig_herb_genes, file = file.path(drt_out, "sig_herb_genes.txt"), 
            quote = FALSE, sep = "\t")

# -------<======= Within-lineage treatment response ========>-------- #

message(">>> generating within-lineage treatment response for drought lineage")
message('contrast = list(c("trt_DW_vs_WW","lineageD.trtDW"))')
message("alpha: 0.05")

res_drt_trt <- results(
  dds_int,
  contrast = list(c("trt_DW_vs_WW", "lineageD.trtDW")),
  alpha = 0.05
)

summary(res_drt_trt)

resLFC_trt_drt <- lfcShrink(dds_int,
                              res = res_drt_trt,
                              type = "ashr")

summary(resLFC_trt_drt)

### Significant genes for drought lineage contrast
sig_trt_drt <- as.data.frame(resLFC_trt_drt)
sig_trt_drt <- sig_trt_drt[which(sig_trt_drt$padj < 0.01 &
                                      abs(sig_trt_drt$log2FoldChange) > 2), ]
sig_trt_drt <- sig_trt_drt[order(sig_trt_drt$padj), ]

# View them
print(sig_trt_drt)

# save them
write.table(sig_trt_drt, file = file.path(drt_out, "sig_trt_drt.txt"),
            quote = FALSE, sep = "\t")

message(">>> generating within-lineage treatment response for herbicide lineage")
message('contrast = list(c("trt_DW_vs_WW", "lineageH.trtDW"))')
message("alpha: 0.05")

res_hrb_trt <- results(
  dds_int,
  contrast = list(c("trt_DW_vs_WW", "lineageH.trtDW")),
  alpha = 0.05
)

summary(res_hrb_trt)

resLFC_trt_hrb <- lfcShrink(dds_int,
                              res = res_hrb_trt,
                              type = "ashr")

summary(resLFC_trt_hrb)

### Significant genes for herbicide lineage contrast
sig_trt_hrb <- as.data.frame(resLFC_trt_hrb)
sig_trt_hrb <- sig_trt_hrb[which(sig_trt_hrb$padj < 0.01 &
                                      abs(sig_trt_hrb$log2FoldChange) > 2), ]
sig_trt_hrb <- sig_trt_hrb[order(sig_trt_hrb$padj), ]

# View them
print(sig_trt_hrb)

# save them
write.table(sig_trt_hrb, file = file.path(drt_out, "sig_trt_hrb.txt"), 
            quote = FALSE, sep = "\t")


# Save the subset results objects
message(">>> saving results objects")

saveRDS(resLFC_interaction_hrb, file = file.path(drt_out, "resLFC_interaction_hrb.rds"))
saveRDS(resLFC_interaction_drt, file = file.path(drt_out, "resLFC_interaction_drt.rds"))
saveRDS(resLFC_trt_drt, file = file.path(drt_out, "resLFC_trt_drt.rds"))
saveRDS(resLFC_trt_hrb, file = file.path(drt_out, "resLFC_trt_hrb.rds"))
message(">>> DESeq2 analysis complete")