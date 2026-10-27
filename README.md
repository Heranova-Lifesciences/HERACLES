# HERACLES

Hera Engine for small RNA Analysis with Clustering and Expression Signatures

tsRNA analysis pipeline: QC → collapse → annotation → cluster → DESeq2 → extend → predict.

## Quick Install

**Option A — conda, no Docker:**

```bash
git clone https://github.com/Heranova-Lifesciences/HERACLES
cd HERACLES
conda env create -f docker/environment.yml
conda activate heracles
```

**Option B — Docker, pinned stack in an image:**

```bash
git clone https://github.com/Heranova-Lifesciences/HERACLES
cd HERACLES
docker build -f docker/Dockerfile -t heracles:1.0 .
```

Both paths install the same versions: the image is built from the same `environment.yml` that Option A hands to conda, so a local install and the container cannot drift apart. See [docker/README.md](docker/README.md) for build, run and Singularity instructions.

## Quick Start

```bash
# Default (QC → DESeq2)
python main.py \
    --fastq-list samples.txt \
    --index-dir /path/to/tRNA_index/ \
    --metadata metadata.tsv \
    --contrast Treat Control

# Full pipeline with prediction
python main.py \
    --fastq-list samples.txt \
    --index-dir /path/to/tRNA_index/ \
    --metadata metadata.tsv \
    --contrast Treat Control \
    --stages full

# Skip clustering, DESeq2 on raw counts
python main.py \
    --fastq-list samples.txt \
    --index-dir /path/to/tRNA_index/ \
    --metadata metadata.tsv \
    --contrast Treat Control \
    --stages full_without_cluster
```

## Pipeline Stages

| Stage | Description |
|-------|-------------|
| `qc` | Trim Galore quality trimming + MultiQC |
| `collapse` | Sequence deduplication |
| `annotation` | Bowtie alignment → tsRNA classification |
| `cluster` | k-mer + alignment clustering |
| `deseq2` | Differential expression + volcano/heatmap |
| `extend` | Extend DEG tsRNAs ±N nt + primer search |
| `predict` | RIsearch2 target prediction + GO/KEGG enrichment |

Stage aliases: `full` (all 7 stages), `full_without_cluster` (skip cluster).

## Resume from Break

Pass the same `--output-dir` and only the remaining stages:

```bash
python main.py \
    --output-dir ./HERACLES_results \
    --stages deseq2,extend,predict \
    --metadata metadata.tsv \
    --contrast Treat Control \
    --index-dir /path/to/tRNA_index/
```

## Input Files

- `samples.txt` — one FASTQ path per line
- `metadata.tsv` — `sample_name<TAB>condition`, no header
- `--index-dir` — Bowtie index directory with tRNA reference FASTA

## More

See `HERACLES_User_Manual.md` for full documentation including per-module usage, output structure, and FAQ.
