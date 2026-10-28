## Methods

### Reference preparation
The *Brachypodium hybridum* assembly (Bhybridum_463_v1.0) and gene annotation (v1.1, Phytozome) were cleaned with AGAT, which added the missing exon features, and converted to a GTF for STAR indexing. All chromosomes (`BhD*`/`BhS*`) and unplaced scaffolds were retained. For counting, each gene model was extended 1,000 bp at its 3' end (strand-aware, clipped at sequence ends; `bedtools slop`) to capture 3'-biased QuantSeq reads, and written as a SAF file. The STAR index was built from the unextended GTF (`--sjdbOverhang 74`, `--genomeSAindexNbases 13`).

### Read processing and alignment
Single-end QuantSeq reads (177 samples: 33 pinoxaden, 144 drought) were trimmed with fastp (Illumina adapter removal, 12 bp trimmed from the 5' end, quality trimming of both ends, minimum length 20 bp). Reads were aligned with STAR v2.7.11b, allowing up to 20 multimapping loci (`--outFilterMultimapNmax 20`, `--winAnchorMultimapNmax 100`), at most 14 mismatches and a mismatch/length ratio of 0.05, with non-canonical unannotated junctions removed. Sorted BAMs were indexed with SAMtools v1.19.2.

### Quantification
Reads were counted per extended gene model with featureCounts (Subread; `-s 1`, `-O --largestOverlap`). Both unique-only and multimapper-inclusive (`-M --fraction`) count tables were generated; the fractional multimapper-inclusive counts were rounded to integers and used for all downstream analyses.

### Differential expression (DESeq2)
Genes with at least 10 counts in at least 3 samples were retained. Each experiment was analyzed with an additive model (main treatment effect) and an interaction model. Genes were called significant at adjusted p < 0.01 and |log2FC| > 2. Fold changes for interaction and within-lineage contrasts were shrunk with `ashr` before thresholding.

- **Pinoxaden experiment** (`~ line + trt` and `~ line + trt + line:trt`; lines Control, Drought, Herbicide; treatment Control vs. Pinoxaden): overall treatment effect (tested against an LFC threshold of 2 at alpha 0.01), line-by-treatment interactions for the Drought and Herbicide lines, and baseline line differences in untreated plants relative to Control. Results are merged into `master_results_all_contrasts.csv`.
- **Drought experiment** (`~ lineage + trt` and `~ lineage + trt + lineage:trt`; lineages C, D, H; treatment WW vs. DW): overall drought effect, interaction terms for D and H, and each lineage's own WW-to-DW response (treatment main effect plus the lineage interaction term).

### Co-expression network analysis (WGCNA)
Using the drought dataset, counts were variance-stabilized (`vst`, `blind = FALSE`) and a single signed network was built from all samples (all lineages, both treatments) with `blockwiseModules`: Pearson correlation, signed TOM, soft-threshold power 5, minimum module size 30, merge cut height 0.25. Module eigengenes were correlated with treatment and measured phenotypes (final height, tiller number, dry shoot/spike/seed weight, seed number). For each gene we calculated module membership (MM) and gene significance (GS) for each phenotype.

**Candidate definitions**
- *Hub gene*: MM to its own module >= 0.8.
- *Shared drought hub*: hub gene that is significantly drought-responsive (WW to DW) in both the D and H lineages.
- *Priority candidate*: shared drought hub with |GS| >= 0.5 for at least one phenotype.

### Software
fastp 1.3.3, STAR 2.7.11b, SAMtools 1.19.2, Subread/featureCounts , AGAT 1.4.3, SeqKit 2.3.0, BEDTools 2.31.0, R with DESeq2, WGCNA, and tidyverse. Pipeline steps were run as SLURM jobs (array jobs for trimming and alignment).