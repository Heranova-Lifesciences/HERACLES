# HERACLES User Manual

> Hera Engine for small RNA Analysis with CLustering and Expression Signatures

HERACLES is an integrated tsRNA (tRNA-derived small RNA) analysis pipeline covering the full workflow from raw sequencing data to differential expression, target prediction, and enrichment analysis.

---

## Table of Contents

1. [Pipeline Overview](#pipeline-overview)
2. [Dependencies](#dependencies)
3. [Quick Start](#quick-start)
4. [Input File Formats](#input-file-formats)
5. [Command-Line Arguments](#command-line-arguments)
6. [Output Files](#output-files)
7. [Running Modules Individually](#running-modules-individually)
8. [FAQ](#faq)

---

## Pipeline Overview

```
FASTQ
 │
 ├─[Step 1] QC ────────────────── Trim Galore + MultiQC
 │  └── trimmed FASTQ + QC reports
 │
 ├─[Step 2] Collapse ──────────── sequence deduplication
 │  └── collapsed FASTA (per sample)
 │
 ├─[Step 3] tsRNA Annotation ──── Bowtie alignment → tsRNA classification (tRF-5/3/i, tiRNA)
 │  └── annotation.tsv + counts.tsv (per sample)
 │
 ├─[Step 4] Merge ─────────────── merge all samples into counts_matrix.tsv
 │
 ├─[Step 5] Clustering ────────── k-mer + sequence alignment clustering (redundancy reduction)
 │  └── clustered_counts_for_DEG.csv + cluster_summary_*.csv
 │
 ├─[Step 6] DESeq2 ────────────── differential expression analysis
 │  └── deseq2_results.tsv + volcano plot + heatmap
 │
 ├─[Step 7] Extend (optional) ─── extend DEG tsRNA ±N nt + search similar sequences
 │  └── output_result.csv + clustered_tsRNA_similar_sequences.csv
 │
 └─[Step 8] Predict (optional) ── RIsearch2 target prediction + GO/KEGG enrichment
    └── gene_stats.txt + enrichment results
```

---

## Dependencies

### Python Dependencies

```bash
pip install -r requirements.txt
```

### External Tools

| Tool | Purpose | Installation |
|------|---------|-------------|
| **Trim Galore** | Quality trimming (Step 1) | `conda install -c bioconda trim-galore` |
| **Bowtie** | Sequence alignment (Step 3) | `conda install -c bioconda bowtie` |
| **MultiQC** | QC reporting (Step 1) | Included in Python dependencies |
| **RIsearch2** | Target prediction (Step 8) | https://github.com/RTH-tools/RIsearch2 |

### Bowtie Index

The `--index-dir` directory must contain Bowtie index files and the tRNA reference FASTA:

```
/path/to/index/
├── mature.1.ebwt          # Bowtie index files
├── mature.2.ebwt
├── mature.3.ebwt
├── mature.4.ebwt
├── mature.rev.1.ebwt
├── mature.rev.2.ebwt
└── hg38-tRNA.fa           # tRNA reference sequences (recommended)
```

> With other naming conventions, the pipeline auto-detects indices prefixed with `hg38-tRNA*` or `mature*`.

---

## Quick Start

### 0. Installation

```bash
git clone https://github.com/Heranova-Lifesciences/HERACLES
cd HERACLES
conda create -n HERACLES
conda activate HERACLES
pip install -r requirements.txt
```

> Also install external tools: Trim Galore, Bowtie, and optionally RIsearch2 (see [Dependencies](#dependencies)).

### 1. Default Run (QC → Collapse → Annotation → Cluster → DESeq2)

```bash
python main.py \
    --fastq-list samples.txt \
    --index-dir /path/to/bowtie_index/ \
    --metadata metadata.tsv \
    --contrast Treat Control \
    --output-dir ./HERACLES_results
```

### 2. Full Run with Target Prediction

```bash
python main.py \
    --fastq-list samples.txt \
    --index-dir /path/to/bowtie_index/ \
    --metadata metadata.tsv \
    --contrast Treat Control \
    --output-dir ./HERACLES_results \
    --stages qc,collapse,annotation,cluster,deseq2,extend,predict \
    --risearch-path /path/to/RIsearch2 \
    --predict-index /path/to/target.suf
```

### 3. Resume from a Specific Stage

```bash
# QC and Collapse already done — start from Annotation
python main.py \
    --stages annotation,cluster,deseq2 \
    --collapsed-dir ./HERACLES_results/collapse_results \
    --index-dir /path/to/bowtie_index/ \
    --metadata metadata.tsv \
    --contrast Treat Control

# Only re-run DESeq2 and extend
python main.py \
    --stages deseq2,extend \
    --index-dir /path/to/bowtie_index/ \
    --metadata metadata.tsv \
    --contrast Treat Control
```

---

## Input File Formats

### `samples.txt` — FASTQ File List

One FASTQ path per line. Lines starting with `#` are ignored:

```
# One FASTQ path per line (supports .fastq / .fq / .fastq.gz / .fq.gz)
/path/to/sample_Treat_rep1.fastq.gz
/path/to/sample_Treat_rep2.fastq.gz
/path/to/sample_Control_rep1.fastq.gz
/path/to/sample_Control_rep2.fastq.gz
```

### `metadata.tsv` — Experimental Metadata

Two TAB-separated columns: **sample name** and **condition**. No header row. Sample names must match FASTQ filenames after stripping extensions:

```
sample_Treat_rep1	Treat
sample_Treat_rep2	Treat
sample_Control_rep1	Control
sample_Control_rep2	Control
```

> **Sample name matching**: if the FASTQ file is `sample_Treat_rep1_trimmed.fq.gz`, the sample name becomes `sample_Treat_rep1` (the pipeline strips `.fq.gz`, `_trimmed`, `_val_1`, and `_collapsed` suffixes automatically).

### `--contrast` — Comparison Groups

Format: `--contrast <treatment> <control>`. Order matters:
```bash
--contrast Treat Control   # Treat vs Control; positive log2FC = upregulated in Treat
```

---

## Command-Line Arguments

### Core Arguments

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `--fastq-list` | Cond.* | — | Path to FASTQ file list |
| `--fastq-dir` | Cond.* | — | Directory of FASTQ files (alternative to `--fastq-list`) |
| `--index-dir` | Yes | — | Bowtie index directory |
| `--metadata` | Cond.† | — | Metadata TSV with sample conditions |
| `--contrast` | Cond.† | — | Comparison: `TREATMENT CONTROL` |
| `--output-dir` | No | `HERACLES_output` | Root output directory |

\* Required for QC/Collapse stages  
† Required for DESeq2 stage

### `--stages` — Pipeline Stages

A comma-separated list of stages to execute. Only the specified stages will run.

**Available stages:**

| Stage | Description |
|-------|-------------|
| `qc` | Trim Galore quality trimming + MultiQC report |
| `collapse` | Sequence deduplication → collapsed FASTA |
| `annotation` | Bowtie alignment → tsRNA classification (tRF-5/3/i, tiRNA) |
| `cluster` | K-mer + alignment-based sequence clustering |
| `deseq2` | Differential expression analysis with volcano/heatmap |
| `extend` | Extend significant DEG tsRNA sequences by ±N nt |
| `primer` | Search count matrix for similar sequences |
| `predict` | RIsearch2 target prediction + GO/KEGG enrichment |

Default: `qc,collapse,annotation,cluster,deseq2`

**Examples:**

```bash
# Core analysis only (default)
--stages qc,collapse,annotation,cluster,deseq2

# Skip QC, use pre-trimmed FASTQs
--stages collapse,annotation,cluster,deseq2

# Only clustering + DESeq2 (count matrix already exists)
--stages cluster,deseq2

# Full analysis with prediction
--stages qc,collapse,annotation,cluster,deseq2,extend,predict
```

### Tuning Parameters

| Argument | Default | Description |
|----------|---------|-------------|
| `--min-len` | 18 | Minimum tsRNA fragment length |
| `--max-len` | 50 | Maximum tsRNA fragment length |
| `--mismatch` | 0 | Bowtie mismatches allowed |
| `--threads` | 4 | Threads for Bowtie and QC |
| `--cluster-method` | `directional` | Clustering method: `directional` (stricter) or `cluster` |
| `--min-count-deseq2` | 10 | Minimum total count for DESeq2 filtering |
| `--extend-by` | 5 | Nucleotides to extend on each side of the sequence |
| `--energy` | -27 | RIsearch2 free energy threshold (kcal/mol) |
| `--threshold` | 0.5 | Gene frequency threshold for enrichment analysis |

### Misc

| Argument | Description |
|----------|-------------|
| `--collapsed-dir` | Directory of pre-existing collapsed FASTA files (use with `--skip-collapse`) |
| `--keep-temp` | Keep temporary intermediate files |
| `--tRNA-fasta` | Manually specify the tRNA reference FASTA (for the extend step) |

---

## Output Files

After a successful run, the following directory structure is produced under `--output-dir`:

```
HERACLES_output/
│
├── qc_results/                         # Step 1: QC
│   ├── trimmed/                        #  trimmed FASTQ files
│   ├── reports/                        #  MultiQC report + QC summary
│   │   ├── multiqc_report.html
│   │   └── qc_summary.csv
│   └── ...
│
├── collapse_results/                   # Step 2: Collapse
│   ├── <sample>_collapsed.fasta        #  deduplicated FASTA per sample
│   └── collapse_statistics.csv         #  aggregate statistics
│
├── tsRNA_results/                      # Step 3: Annotation
│   ├── <sample>_annotation.tsv         #  annotation per sample
│   ├── <sample>_counts.tsv             #  counts per sample
│   └── ...
│
├── counts_matrix.tsv                   # Step 4: Merged matrix
│
├── cluster_results/                    # Step 5: Clustering
│   ├── clustered_counts_for_DEG.csv    #  clustered matrix (input for DESeq2)
│   ├── cluster_summary_detailed.csv    #  detailed cluster report
│   └── cluster_summary_simplified.csv  #  simplified cluster summary
│
├── deseq2_results/                     # Step 6: DESeq2
│   ├── deseq2_results.tsv             #  full differential expression results
│   ├── volcano_plot.png               #  volcano plot
│   ├── heatmap.png                    #  heatmap (top 50 significant genes)
│   └── sample_counts/                 #  per-sample clustered counts
│
├── extend_primer_results/ (--stages extend)  # Step 7: Extend + Primer
│   ├── significant_tsRNAs.txt         #  list of significant DEG tsRNAs
│   ├── output_result.csv             #  extended tsRNA sequences
│   └── clustered_tsRNA_similar_sequences.csv
│
└── prediction_results/ (--stages predict) # Step 8: Prediction
    ├── all_results.txt                #  raw RIsearch2 output
    ├── gene_stats.txt                 #  gene target statistics
    ├── selected_genes.txt             #  high-frequency target genes
    └── enrichment/                    #  GO/KEGG enrichment results
        ├── GO_Biological_Process_2025_enrichment.csv
        ├── GO_Biological_Process_2025_dotplot.png
        ├── GO_Cellular_Component_2025_enrichment.csv
        ├── GO_Molecular_Function_2025_enrichment.csv
        └── KEGG_enrichment.csv
```

---

## Running Modules Individually

Each module can also be run standalone.

### QC Module

```bash
python modules/qc.py -i /path/to/fastq_dir/ -o qc_results -t 8
```

### Collapse Module

```bash
python modules/collapse.py trimmed_sample.fq.gz
```

### tsRNA Annotation Module

```bash
# Single sample
python modules/tsRNA_annotation.py -i sample.fasta -d /path/to/index/ -n Sample1

# Multiple samples (file format: sample_name<tab>fasta_path)
python modules/tsRNA_annotation.py -s sample_list.txt -d /path/to/index/ -t 8
```

### Clustering Module

```bash
python modules/cluster.py -i counts_matrix.tsv -o cluster_results -m directional
```

### DESeq2 Module

```bash
python modules/run_deseq2.py \
    -m metadata.tsv \
    -o deseq2_results \
    --contrast Treat Control \
    --min-count 10
```

### Extend & Primer Module

```bash
# Extend tsRNA sequences only
python modules/extend.py extend \
    --fa hg38-tRNA.fa \
    --list significant_tsRNAs.txt \
    --csv cluster_summary_simplified.csv \
    --out output_result.csv

# Primer search only
python modules/extend.py primer \
    --extend output_result.csv \
    --count-matrix counts_matrix.tsv \
    --out similar_sequences.csv

# Both together
python modules/extend.py all \
    --fa hg38-tRNA.fa \
    --list significant_tsRNAs.txt \
    --csv cluster_summary_simplified.csv \
    --count-matrix counts_matrix.tsv \
    --out-dir extend_primer_results
```

### Prediction Module

```bash
python modules/prediction_enrichment.py \
    --list significant_tsRNAs.txt \
    --index target.suf \
    --energy -27 \
    --threads 8 \
    --risearch_path /path/to/RIsearch2 \
    --output_dir prediction_results
```

---

## FAQ

### Q: "No FASTQ files found"
**A:** Check that the `--fastq-list` file path is correct and that each line points to an existing FASTQ file.

### Q: "Bowtie index not found"
**A:** Ensure the `--index-dir` contains all `.ebwt` index files (typically 6 files).

### Q: "DESeq2 requires --metadata and --contrast"
**A:** Both `--metadata` and `--contrast` are required when running the DESeq2 stage.

### Q: Sample names in metadata don't match FASTQ filenames
**A:** The pipeline extracts sample names by:
1. Stripping file extensions (`.fastq.gz` → `.fastq` → `.fq.gz` → `.fq` → `.gz`)
2. Stripping `_trimmed` and `_val_1` suffixes
3. Stripping `_collapsed` suffix

For example, `Liver_Treat_S1.fastq.gz` yields sample name `Liver_Treat_S1`.

### Q: How do I skip clustering and run DESeq2 on raw counts?
**A:** Omit `cluster` from `--stages`, e.g. `--stages qc,collapse,annotation,deseq2`. The pipeline falls back to the raw `counts_matrix.tsv`.

### Q: Extend step fails with "tRNA FASTA not found"
**A:** Use `--tRNA-fasta` to explicitly specify the tRNA reference FASTA file path. You can also point `--index-dir` to a directory containing `hg38-tRNA.fa`.

### Q: DESeq2 reports "conditions have 0 samples"
**A:** The sample names in `--metadata` do not match the column names in the count matrix. Verify that metadata sample names match the FASTQ filenames after suffix stripping.

---

## Citing

If you use HERACLES, please cite the underlying tools:
- **Trim Galore**: https://github.com/FelixKrueger/TrimGalore
- **Bowtie**: Langmead et al., 2009
- **pydeseq2**: Muzellec et al., 2023
- **RIsearch2**: Alkan et al., 2017
- **gseapy**: https://github.com/zqfang/gseapy
- **MultiQC**: Ewels et al., 2016
