# HERACLES User Manual

> Hera Engine for small RNA Analysis with Clustering and Expression Signatures

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

### External Tools

| Tool | Purpose | Installation |
|------|---------|-------------|
| **Trim Galore** | Quality trimming (Step 1) | `conda install -c bioconda trim-galore=0.6.10` |
| **Bowtie** | Sequence alignment (Step 3) | `conda install -c bioconda bowtie` |
| **MultiQC** | QC reporting (Step 1) | Included in Python dependencies |
| **RIsearch2** | Target prediction (Step 8) | `conda install -c rthtools risearch` |

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

### RIsearch2 Index

For target prediction (`--stages predict`), RIsearch2 requires a target sequence index. Two reference databases are included under `RIsearch2_index/`:

```
RIsearch2_index/
├── GRCh38.CDS.fa           # Human coding sequences (CDS)
└── GRCh38.3UTR.fa          # Human 3' UTR sequences
```

**How it works:**

- Use `--predict-index CDS` or `--predict-index 3UTR` to select the built-in database. The pipeline automatically builds a `.suf` index on first use (cached for subsequent runs).
- You can also provide a custom `.suf` index: `--predict-index /path/to/custom.suf`.

**Building a custom index:**

```bash
RIsearch2 -c your_targets.fa -o your_targets.suf
```

---

## Quick Start

### 0. Installation

```bash
git clone https://github.com/Heranova-Lifesciences/HERACLES
cd HERACLES
conda create -n HERACLES -c conda-forge -c bioconda -c rthtools python=3.10 pip trim-galore=0.6.10 bowtie risearch
conda activate HERACLES
pip install -r requirements.txt
```

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
    --stages full \
    --predict-index CDS
```

### 3. Full Run Without Clustering

```bash
# DESeq2 directly on raw count matrix (counts_matrix.tsv), no clustering
python main.py \
    --fastq-list samples.txt \
    --index-dir /path/to/bowtie_index/ \
    --metadata metadata.tsv \
    --contrast Treat Control \
    --output-dir ./HERACLES_results \
    --stages full_without_cluster
```

### 4. Resume from a Specific Stage

If a pipeline run breaks midway, simply re-run with the same `--output-dir` and specify only the remaining stages. The pipeline auto-detects existing results from completed stages within the output directory — no need to manually specify intermediate file paths.

```bash
# Pipeline broke during Annotation — re-run from Annotation onwards
python main.py \
    --fastq-list samples.txt \
    --index-dir /path/to/bowtie_index/ \
    --metadata metadata.tsv \
    --contrast Treat Control \
    --output-dir ./HERACLES_results \
    --stages annotation,cluster,deseq2

# QC and Collapse already done, re-run from Annotation
python main.py \
    --output-dir ./HERACLES_results \
    --stages annotation,cluster,deseq2 \
    --index-dir /path/to/bowtie_index/ \
    --metadata metadata.tsv \
    --contrast Treat Control

# Only re-run DESeq2 and extend
python main.py \
    --output-dir ./HERACLES_results \
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
| `extend` | Extend significant DEG tsRNA sequences by ±N nt + primer search |
| `predict` | RIsearch2 target prediction + GO/KEGG enrichment |

Default: `qc,collapse,annotation,cluster,deseq2`

**Stage aliases:**

| Alias | Equivalent to |
|-------|---------------|
| `full` | `qc,collapse,annotation,cluster,deseq2,extend,predict` |
| `full_without_cluster` | `qc,collapse,annotation,deseq2,extend,predict` (DESeq2 uses raw `counts_matrix.tsv` on individual tsRNAs, skipping clustering) |

> `full_without_cluster` skips the clustering step. In this mode, DESeq2 runs directly on the raw count matrix (`counts_matrix.tsv`) where each row is an individual tsRNA, rather than on clustered groups.

**Examples:**

```bash
# Core analysis only (default)
--stages qc,collapse,annotation,cluster,deseq2

# Full analysis with prediction (alias)
--stages full

# Full analysis without clustering (alias)
--stages full_without_cluster

# Skip QC, use pre-trimmed FASTQs
--stages collapse,annotation,cluster,deseq2

# Only clustering + DESeq2 (count matrix already exists)
--stages cluster,deseq2

# Custom subset
--stages qc,collapse,annotation,cluster,deseq2,extend,predict
```

### Tuning Parameters

| Argument | Default | Description |
|----------|---------|-------------|
| `--qc-threads` | 4 | Threads for Trim Galore |
| `--qc-quality` | 20 | Phred quality trimming threshold |
| `--qc-length` | 18 | Minimum read length after trimming |
| `--trim-galore-path` | `trim_galore` | Path to trim_galore executable |
| `--adapter` | `""` | Adapter sequence (auto-detected if empty) |
| `--min-count` | 1 | Minimum read count for collapse retention |
| `--bowtie-path` | `bowtie` | Path to bowtie executable |
| `--min-len` | 18 | Minimum tsRNA fragment length |
| `--max-len` | 50 | Maximum tsRNA fragment length |
| `--mismatch` | 0 | Bowtie mismatches allowed |
| `--threads` | 4 | Threads for Bowtie and QC |
| `--cluster-method` | `directional` | Clustering method: `directional` (stricter) or `cluster` |
| `--design` | `condition` | Design formula factor name |
| `--min-count-deseq2` | 10 | Minimum total count for DESeq2 filtering |
| `--top-n-heatmap` | 50 | Top N genes shown in heatmap |
| `--pvalue-thresh` | 0.05 | P-value threshold for significance |
| `--output-normalized` | `counts_matrix_normalized.tsv` | Path to save DESeq2 normalized count matrix (TSV) |
| `--extend-by` | 5 | Nucleotides to extend on each side of the sequence |
| `--energy` | -27 | RIsearch2 free energy threshold (kcal/mol) |
| `--predict-index` | `CDS` | RIsearch2 index: CDS, 3UTR, or path to .suf |
| `--threshold` | — | Gene frequency threshold for enrichment |
| `--top-percent` | 10 | Select top N% of genes by target count (default) |

### Misc

| Argument | Description |
|----------|-------------|
| `--collapsed-dir` | Directory of pre-existing collapsed FASTA files (used when `collapse` is not in `--stages`) |
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

Each module can also be run standalone, in order or independently. Below are detailed usage instructions including input requirements, all CLI options, and output files.

---

### QC Module (`modules/qc.py`)

**Purpose:** Quality trimming via Trim Galore + MultiQC report aggregation.

**Prerequisites:** Trim Galore must be installed and available on PATH (or specified via `--trim_galore_path`).

**Input:** A FASTQ directory, a single FASTQ file, or a `.txt` list file (one FASTQ path per line).

**Input file formats:**

```
# Option A: .txt list file
# samples.txt
/path/to/sample1.fastq.gz
/path/to/sample2.fastq.gz

# Option B: directory containing .fastq / .fastq.gz / .fq / .fq.gz files
/path/to/fastq_dir/

# Option C: single FASTQ file
/path/to/sample.fastq.gz
```

**Arguments:**

| Flag | Default | Description |
|------|---------|-------------|
| `-i`, `--input` | (required) | Input path: directory, single FASTQ, or `.txt` list file |
| `-o`, `--output_dir` | `qc_results` | Output directory |
| `-t`, `--threads` | 4 | Threads for Trim Galore |
| `-q`, `--quality` | 20 | Phred quality trimming threshold |
| `-l`, `--length` | 18 | Minimum read length after trimming |
| `--trim_galore_path` | `trim_galore` | Path to trim_galore executable |
| `--multiqc_path` | `multiqc` | Path to multiqc executable |
| `--adapter` | `""` | Adapter sequence (auto-detected if empty) |

**Output:**

```
qc_results/
├── trimmed/                 # Trimmed FASTQ files (*_trimmed.fq.gz)
├── reports/
│   ├── multiqc_report.html  # MultiQC HTML report
│   └── qc_summary.csv       # QC statistics summary
└── ...
```

**Examples:**

```bash
# Process a directory of FASTQ files with 8 threads
python modules/qc.py -i /path/to/fastq_dir/ -o qc_results -t 8

# Process a list file with custom adapter
python modules/qc.py -i samples.txt -o qc_results --adapter TGGAATTCTCGGGTGCCAAGG

# Single FASTQ, stricter trimming
python modules/qc.py -i sample.fastq.gz -o qc_results -q 25 -l 20
```

---

### Collapse Module (`modules/collapse.py`)

**Purpose:** Deduplicate reads from a trimmed FASTQ file, output a collapsed FASTA with unique sequences and their counts.

**Prerequisites:** Trimmed FASTQ file (from QC module or pre-existing).

**Input:** A single trimmed FASTQ file (`.fastq` / `.fq` / `.fastq.gz` / `.fq.gz`).

**Arguments:**

| Flag | Default | Description |
|------|---------|-------------|
| *(positional)* | (required) | Path to trimmed FASTQ file |
| *(hardcoded)* | `collapse_results` | Output directory |
| *(hardcoded)* | 1 | Minimum read count |

**Output:**

```
collapse_results/
└── <sample>_collapsed.fasta    # FASTA with unique sequences and counts in headers
```

> FASTA header format: `>sequence_X-<count>` where X is a running index and count is the number of reads.

**Example:**

```bash
python modules/collapse.py qc_results/trimmed/sample1_trimmed.fq.gz
```

---

### tsRNA Annotation Module (`modules/tsRNA_annotation.py`)

**Purpose:** Align collapsed FASTA sequences to the tRNA reference via Bowtie, classify tsRNA types (tRF-5, tRF-3, tRF-i, tiRNA-5, tiRNA-3, tiRNA-5L, tiRNA-3L), and generate annotated count tables.

**Prerequisites:**
- Bowtie must be installed and on PATH (or specified via `--bowtie-path`).
- Bowtie index directory (see [Bowtie Index](#bowtie-index)).
- Collapsed FASTA file(s) (from collapse module).

**Input:** Single collapsed FASTA file or a sample list file.

```
# Sample list file format (TAB-separated, no header):
Sample1    /path/to/sample1_collapsed.fasta
Sample2    /path/to/sample2_collapsed.fasta
```

**Arguments:**

| Flag | Default | Description |
|------|---------|-------------|
| `-i`, `--input` | — | Single input FASTA file |
| `-s`, `--samples` | — | Sample list file (mutually exclusive with `-i`) |
| `-d`, `--index-dir` | (required) | tRNA Bowtie index directory |
| `-n`, `--sample-name` | (auto) | Sample name (used with `-i`; auto-detected from filename) |
| `-o`, `--output-dir` | `tsRNA_results` | Output directory |
| `--min-len` | 18 | Minimum fragment length |
| `--max-len` | 45 | Maximum fragment length |
| `-v`, `--mismatch` | 0 | Bowtie mismatches allowed |
| `-t`, `--threads` | 4 | Number of threads |
| `--bowtie-path` | `bowtie` | Path to bowtie executable |
| `--keep-temp` | off | Keep temporary intermediate files |
| `--debug` | off | Enable debug logging |

**Output (per sample):**

```
tsRNA_results/
├── <sample>_annotation.tsv    # Per-locus annotation
├── <sample>_counts.tsv        # Final tsRNA count table
└── <sample>_aligned.sam       # Raw Bowtie alignment (SAM)
```

> `<sample>_counts.tsv` header format: each tsRNA ID follows `tRNA:start-end:type:sequence`.

**Examples:**

```bash
# Single sample
python modules/tsRNA_annotation.py \
    -i collapse_results/Sample1_collapsed.fasta \
    -d /path/to/bowtie_index/ \
    -n Sample1

# Multiple samples via list file
python modules/tsRNA_annotation.py \
    -s sample_list.txt \
    -d /path/to/bowtie_index/ \
    -t 8

# Allow 1 mismatch
python modules/tsRNA_annotation.py \
    -i sample.fasta -d /path/to/index/ -v 1
```

---

### Clustering Module (`modules/cluster.py`)

**Purpose:** Cluster tsRNA sequences by k-mer overlap and terminal alignment similarity to reduce redundancy. Outputs a clustered count matrix and cluster summary tables.

**Prerequisites:** Count matrix `counts_matrix.tsv` (produced by merging per-sample `*_counts.tsv` from annotation module).

**Input:** `counts_matrix.tsv` — a TAB-separated file with `tsRNA_id` as the first column (index) and sample names as remaining columns.

```
# counts_matrix.tsv format:
tsRNA_id	Sample1	Sample2	Sample3
tRNA-Gly-GCC:1-32:tRF-5:GCATTGGTGGT...	120	85	210
tRNA-Val-CAC:5-28:tRF-5:GTTTCCGTAGT...	45	32	67
...
```

**Arguments:**

| Flag | Default | Description |
|------|---------|-------------|
| `-i`, `--input` | (required) | Count matrix TSV (`counts_matrix.tsv`) |
| `-o`, `--output-dir` | `.` | Output directory |
| `-m`, `--method` | `directional` | Clustering method: `directional` (stricter) or `cluster` |
| `--min-shared-kmers` | 5 | Minimum shared k-mers for candidate pairs |
| `--kmer-size` | 10 | K-mer size |
| `--min-length-ratio` | 0.94 | Minimum length ratio for alignment validation |

**Output:**

```
cluster_results/
├── clustered_counts_for_DEG.csv    # Clustered count matrix → input for DESeq2
├── cluster_summary_detailed.csv    # Per-cluster details (cluster_id, type, seq, count, members)
└── cluster_summary_simplified.csv  # Simplified summary (cluster_id, tsRNA_type, seq, total_count)
```

**Examples:**

```bash
# Directional clustering (stricter, recommended)
python modules/cluster.py -i counts_matrix.tsv -o cluster_results -m directional

# Standard clustering with custom k-mer parameters
python modules/cluster.py -i counts_matrix.tsv -o cluster_results -m cluster --kmer-size 8
```

---

### DESeq2 Module (`modules/run_deseq2.py`)

**Purpose:** Run differential expression analysis with pydeseq2, generate volcano plot and heatmap.

**Prerequisites:**
- Per-sample count files (individual `*_counts.tsv`, typically from `split_clustered_matrix`).
- A metadata file listing each count file path and its condition.

**Metadata file format (TAB-separated, no header):**

```
# metadata_deseq2.tsv
/path/to/deseq2_results/sample_counts/Sample1_counts.tsv	Treat
/path/to/deseq2_results/sample_counts/Sample2_counts.tsv	Treat
/path/to/deseq2_results/sample_counts/Sample3_counts.tsv	Control
/path/to/deseq2_results/sample_counts/Sample4_counts.tsv	Control
```

**Arguments:**

| Flag | Default | Description |
|------|---------|-------------|
| `-m`, `--metadata` | (required) | Metadata file (see format above) |
| `-o`, `--output` | `deseq2_results` | Output directory |
| `--design` | `condition` | Design formula factor name |
| `--contrast` | (required) | Comparison: `TREATMENT CONTROL` (order matters) |
| `--min-count` | 10 | Minimum total count to retain a gene |
| `--top-n` | 50 | Top N genes to show in heatmap |
| `--pvalue-thresh` | 0.05 | P-value threshold for significance |
| `--output-normalized` | `counts_matrix_normalized.tsv` | Path to save DESeq2 normalized count matrix (TSV) |

**Output:**

```
deseq2_results/
├── deseq2_results.tsv      # Full DESeq2 results (baseMean, log2FoldChange, pvalue, padj)
├── volcano_plot.png        # Volcano plot
├── heatmap.png             # Expression heatmap of top significant genes
└── sample_counts/          # (input) Per-sample count files
```

**Examples:**

```bash
python modules/run_deseq2.py \
    -m metadata_deseq2.tsv \
    -o deseq2_results \
    --contrast Treat Control \
    --min-count 10

# Custom thresholds
python modules/run_deseq2.py \
    -m metadata_deseq2.tsv \
    -o deseq2_results \
    --contrast Treat Control \
    --min-count 5 \
    --top-n 100 \
    --pvalue-thresh 0.01
```

---

### Extend & Primer Module (`modules/extend.py`)

**Purpose:** Extend significant DEG tsRNA sequences by ±N nucleotides using the tRNA reference genome, then search the count matrix for sequences similar to the extended tsRNAs (primer search).

**Prerequisites:**
- DESeq2 results (`deseq2_results.tsv`) or a list of significant tsRNA IDs.
- tRNA reference FASTA (e.g. `hg38-tRNA.fa`).
- Count matrix (`counts_matrix.tsv`) for primer search.
- Cluster summary CSV (optional, for count lookup).

**tsRNA ID list format (`significant_tsRNAs.txt`):**

```
tRNA-Gly-GCC:1-32:tRF-5:GCATTGGTGGT...
tRNA-Val-CAC:5-28:tRF-5:GTTTCCGTAGT...
```

**Sub-commands:**

```
extend    — Extend tsRNA sequences only
primer    — Primer/similar sequence search only
all       — Extend + primer search together
```

**Arguments (sub-command: `extend`):**

| Flag | Default | Description |
|------|---------|-------------|
| `--fa` | (required) | tRNA reference FASTA |
| `--list` | (required) | tsRNA ID list (one per line) |
| `--csv` | — | Cluster summary CSV for count lookup (optional) |
| `--out` | `output_result.csv` | Output CSV path |
| `--extend-by` | 5 | Nucleotides to extend on each side |

**Arguments (sub-command: `primer`):**

| Flag | Default | Description |
|------|---------|-------------|
| `--extend` | (required) | Extended tsRNA CSV (from `extend` step) |
| `--count-matrix` | (required) | Count matrix TSV (`counts_matrix.tsv`) |
| `--out` | `clustered_tsRNA...csv` | Output CSV path |

**Arguments (sub-command: `all`):**

| Flag | Default | Description |
|------|---------|-------------|
| `--fa` | (required) | tRNA reference FASTA |
| `--list` | (required) | tsRNA ID list |
| `--csv` | — | Cluster summary CSV (optional) |
| `--count-matrix` | (required) | Count matrix TSV |
| `--out-dir` | `extend_primer_results` | Output directory |
| `--extend-by` | 5 | Nucleotides to extend per side |

**Output:**

```
extend_primer_results/
├── output_result.csv                       # Extended tsRNA sequences with counts
└── clustered_tsRNA_similar_sequences.csv   # Similar sequences from count matrix
```

**Examples:**

```bash
# Step 1: Extend tsRNA sequences
python modules/extend.py extend \
    --fa hg38-tRNA.fa \
    --list significant_tsRNAs.txt \
    --csv cluster_summary_simplified.csv \
    --out output_result.csv \
    --extend-by 5

# Step 2: Search similar sequences
python modules/extend.py primer \
    --extend output_result.csv \
    --count-matrix counts_matrix.tsv \
    --out similar_sequences.csv

# All in one
python modules/extend.py all \
    --fa hg38-tRNA.fa \
    --list significant_tsRNAs.txt \
    --csv cluster_summary_simplified.csv \
    --count-matrix counts_matrix.tsv \
    --out-dir extend_primer_results
```

---

### Prediction Module (`modules/prediction_enrichment.py`)

**Purpose:** Run RIsearch2 target prediction on significant tsRNAs, then perform GO and KEGG enrichment analysis on high-frequency target genes.

**Prerequisites:**
- RIsearch2 must be installed (https://github.com/RTH-tools/RIsearch2).
- Target index: built-in `CDS`/`3UTR` (from `RIsearch2_index/`) or custom `.suf` file.
- gseapy Python package (`pip install gseapy`).
- List of significant tsRNA IDs (from DESeq2 results).

**tsRNA ID list format (`significant_tsRNAs.txt`):**

```
tRNA-Gly-GCC:1-32:tRF-5:GCATTGGTGGT...
tRNA-Val-CAC:5-28:tRF-5:GTTTCCGTAGT...
```

**Arguments:**

| Flag | Default | Description |
|------|---------|-------------|
| `--list` | (required) | tsRNA ID list file |
| `--index` | `CDS` | Target index: `CDS`, `3UTR`, or path to custom `.suf` |
| `--energy` | -27 | Free energy threshold (kcal/mol, more negative = stricter) |
| `--threads` | 8 | Threads for RIsearch2 |
| `--output_dir` | `./pipeline_output` | Output directory |
| `--threshold` | — | Fraction of tsRNAs targeting a gene for selection |
| `--top-percent` | 10 | Select top N% of genes by target count (default) |
| `--risearch_path` | `RIsearch2` | Path to RIsearch2 executable |

**Output:**

```
pipeline_output/
├── query.fa                          # FASTA of input tsRNA sequences
├── all_results.txt                   # Raw RIsearch2 results (merged)
├── gene_stats.txt                    # Per-gene target frequency statistics
├── selected_genes.txt                # High-frequency target genes (above threshold)
└── enrichment/                       # GO/KEGG enrichment results
    ├── GO_Biological_Process_2025_enrichment.csv
    ├── GO_Biological_Process_2025_dotplot.png
    ├── GO_Cellular_Component_2025_enrichment.csv
    ├── GO_Molecular_Function_2025_enrichment.csv
    └── KEGG_enrichment.csv
```

**Examples:**

```bash
# Predict targets in CDS regions
python modules/prediction_enrichment.py \
    --list significant_tsRNAs.txt \
    --index CDS \
    --energy -27 \
    --threads 8 \
    --output_dir prediction_results

# Custom index with stricter energy threshold
python modules/prediction_enrichment.py \
    --list significant_tsRNAs.txt \
    --index /path/to/custom_target.suf \
    --energy -30 \
    --threshold 0.3 \
    --risearch_path /path/to/RIsearch2 \
    --output_dir prediction_results
```

---

### Module Dependency Chain

```
FASTQ (.fastq.gz)
  │
  ├── modules/qc.py ──────────────────> trimmed FASTQ
  │     │
  │     └── modules/collapse.py ──────> collapsed FASTA
  │           │
  │           └── modules/tsRNA_annotation.py ──> *_counts.tsv (per sample)
  │                 │
  │                 └── (merge) ──────────────> counts_matrix.tsv
  │                       │
  │                       ├── modules/cluster.py ──> clustered_counts_for_DEG.csv
  │                       │     │
  │                       │     └── modules/run_deseq2.py ──> deseq2_results.tsv
  │                       │           │
  │                       │           ├── modules/extend.py ──> extended + primer results
  │                       │           │
  │                       │           └── modules/prediction_enrichment.py ──> enrichment results
  │                       │
  │                       └── modules/run_deseq2.py (skip cluster) ──> deseq2_results.tsv
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
**A:** Use `--stages full_without_cluster` to run the full pipeline without clustering, or omit `cluster` from a custom stage list, e.g. `--stages qc,collapse,annotation,deseq2`. In both cases, DESeq2 falls back to the raw `counts_matrix.tsv`.

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
