## Methods

### Reference preparation
The *Brachypodium hybridum* assembly (Bhybridum_463_v1.0) and gene annotation (v1.1, Phytozome) were cleaned with AGAT, which added the missing exon features, and converted to a GTF for STAR indexing. All chromosomes (`BhD*`/`BhS*`) and unplaced scaffolds were retained. For counting, each gene model was extended 1,000 bp at its 3' end (strand-aware, clipped at sequence ends; `bedtools slop`) to capture 3'-biased QuantSeq reads, and written as a SAF file. The STAR index was built from the unextended GTF (`--sjdbOverhang 74`, `--genomeSAindexNbases 13`).

### Read processing and alignment
Single-end QuantSeq reads (177 samples: 33 pinoxaden, 144 drought) were trimmed with fastp (Illumina adapter removal, 12 bp trimmed from the 5' end, quality trimming of both ends, minimum length 20 bp). Reads were aligned with STAR v2.7.11b, allowing up to 20 multimapping loci (`--outFilterMultimapNmax 20`, `--winAnchorMultimapNmax 100`), at most 14 mismatches and a mismatch/length ratio of 0.05, with non-canonical unannotated junctions removed. Sorted BAMs were indexed with SAMtools v1.19.2.

### Quantification
Reads were counted per extended gene model with featureCounts (Subread; `-s 1`, `-O --largestOverlap`). Both unique-only and multimapper-inclusive (`-M --fraction`) count tables were generated; the fractional multimapper-inclusive counts were rounded to integers and used for all downstream analyses.

### Differential expression (DESeq2)
Genes with at least 10 counts in at least 3 samples were retained. Each experiment was analyzed with an additive model (main treatment effect) and an interaction model. Fold changes for interaction, baseline and within-lineage contrasts were shrunk with `ashr`. Per-experiment DE lists were called at adjusted p < 0.01 and |log2FC| > 2.

- **Pinoxaden experiment** (`~ line + trt` and `~ line + trt + line:trt`; lines Control, Drought, Herbicide; treatment Control vs. Pinoxaden): overall treatment effect (tested against an LFC threshold of 2 at alpha 0.01), line-by-treatment interactions for the Drought and Herbicide lines, baseline line differences in untreated plants relative to Control, and each line's own Control-to-Pinoxaden response (treatment main effect plus the line interaction term). Results are merged into `master_results_all_contrasts.csv`.
- **Drought experiment** (`~ lineage + trt` and `~ lineage + trt + lineage:trt`; lineages C, D, H; treatment WW vs. DW): overall drought effect, interaction terms for D and H, and each lineage's own WW-to-DW response (treatment main effect plus the lineage interaction term).

### Comparison of drought and pinoxaden responses
To find genes specific to drought and genes shared between the drought and pinoxaden responses, the two experiments were compared using the full shrunken results from both. For consistency, the overall treatment effects (DW vs. WW; Pinoxaden vs. Control) were recomputed from the additive models using a Wald test against zero followed by `ashr` shrinkage, so that every contrast followed the same procedure. Adjusted p-values are from the Wald test of the unshrunken estimate; log2FC values are the shrunken estimates.

Two significance tiers were applied: a lenient tier (adjusted p < 0.05 and |log2FC| > 1), used for classification, and a strict tier (adjusted p < 0.01 and |log2FC| > 2), used to flag high-confidence genes. Overlaps among DE gene sets were visualized as UpSet plots (UpSetR) for each tier and for all genes, upregulated genes only and downregulated genes only; direction-specific plots count a gene only if it moves in that direction, so opposite-direction changes are never counted as shared. Plots were drawn separately for treatment-response sets (overall and within-lineage responses in both experiments) and for lineage effects (interaction and baseline contrasts).

Genes were classified using the treatment-response contrasts only. Within each experiment, a gene was labeled *up* or *down* if it was DE at the lenient tier in that direction, *mixed* if it was up in one response contrast and down in another, *null* if it was tested and clearly not DE (adjusted p > 0.1 and |log2FC| < 0.5 in every contrast), *ambiguous* if tested but neither DE nor clearly null, and *untested* if absent from the results or with an NA adjusted p-value. Genes up (or down) in both experiments were classed as shared up (or down); genes DE in both experiments in opposite directions or with mixed direction as discordant; genes DE in one experiment and clearly null in the other as drought-specific or pinoxaden-specific; and genes DE in one experiment but ambiguous or untested in the other were kept as separate classes. This avoids calling a gene specific merely because it lacked power in the smaller pinoxaden experiment.

As a threshold-free check, Spearman correlations of shrunken log2FC between matched drought and pinoxaden contrasts (D, H and all lineages) were calculated across all genes tested in both experiments, and across genes DE in the drought contrast, together with the fraction of those genes changing in the same direction.

### Functional annotation (InterProScan)
Predicted proteins were annotated with InterProScan [version] [applications used, e.g., Pfam, PANTHER, ProSite, Coils, MobiDB-Lite]. Transcript identifiers were mapped to gene identifiers by removing the isoform number, and annotations were combined across isoforms to give one record per gene: InterPro entries, Pfam domains, GO terms and pathway annotations. Gene-to-GO mappings were exported for downstream enrichment analysis.

### Co-expression network analysis (WGCNA)
Using the drought dataset, counts were variance-stabilized (`vst`, `blind = FALSE`) and a single signed network was built from all samples (all lineages, both treatments) with `blockwiseModules`: Pearson correlation, signed TOM, soft-threshold power 5, minimum module size 30, merge cut height 0.25. Module eigengenes were correlated with treatment and measured phenotypes (final height, tiller number, dry shoot/spike/seed weight, seed number). For each gene we calculated module membership (MM) and gene significance (GS) for each phenotype.

**Candidate definitions**
- *Hub gene*: MM to its own module >= 0.8.
- *Shared drought hub*: hub gene that is significantly drought-responsive (WW to DW) in both the D and H lineages.
- *Priority candidate*: shared drought hub with |GS| >= 0.5 for at least one phenotype.

### Software
fastp 1.3.3, STAR 2.7.11b, SAMtools 1.19.2, Subread/featureCounts [version], AGAT 1.4.3, SeqKit 2.3.0, BEDTools 2.31.0, InterProScan [version], R [version] with DESeq2 [version], ashr [version], WGCNA [version], UpSetR [version], and tidyverse [version]. Pipeline steps were run as SLURM jobs (array jobs for trimming and alignment).