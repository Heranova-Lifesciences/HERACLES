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

## Optional QC filter: adapter-unsupported reads

By default, QC keeps every read that Trim Galore retains. In single-end libraries
whose inserts can fill the full sequencing cycle, a read that is still exactly
`L` bases long after trimming has had no bases removed at all — no 3' adapter was
detected, so its final nucleotide may be a read-length truncation rather than the
authentic RNA 3' terminus. The optional filter removes such reads before
collapsing:

```bash
python main.py \
    --fastq-list samples.txt \
    --index-dir /path/to/tRNA_index/ \
    --metadata metadata.tsv \
    --contrast Treat Control \
    --drop-unsupported
```

- **Off by default.** With the flag off, output is identical to previous releases.
- `L` is inferred **per library** from its raw FASTQ (longest read among the first
  100,000 records). Use `--cycle-length 50` to force a fixed cycle count.
- Original trimmed files stay in `qc_results/trimmed/`; filtered copies are
  written to `qc_results/trimmed_unsupported_filtered/` under the **same file
  names**, so downstream sample names never change.
- Per-library statistics (cycle length, reads kept/removed, percentage) are
  written to `qc_results/reports/unsupported_reads_summary.csv`.

This is a read-level QC filter, applied before collapse and annotation. It is not
the same as deleting 50-nt features from an existing count matrix: reads shorter
than `L` (including 49-nt reads whose terminal A was removed by a 1-bp adapter
overlap) are deliberately kept.

## Pipeline Stages

| Stage | Description |
|-------|-------------|
| `qc` | Trim Galore quality trimming + MultiQC; optional removal of adapter-unsupported full-cycle reads |
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

## Recent changes (November 2026)

- **Optional QC filter** `--drop-unsupported` / `--cycle-length` (see above); off by default, so existing runs are unchanged.
- The clustering abundance rule is exposed as `abund_ratio` in `modules/cluster.py`; the default of `2.0` reproduces the previous hard-coded behaviour.
- Local sequence alignment moved from the deprecated `Bio.pairwise2` to `Bio.Align.PairwiseAligner` with identical scoring (match 2, mismatch −1, gap open −5, gap extend −0.1); requires `biopython>=1.80`.
- Enrichment accepts local GMT files via the `HERACLES_GMT_DIR` environment variable, so GO/KEGG enrichment can run without network access. Without it, the online Enrichr libraries are used as before.

## More

See `HERACLES_User_Manual.md` for full documentation including per-module usage, output structure, and FAQ.

## License

MIT — see [LICENSE](LICENSE). Copyright (c) 2026 Heranova Lifesciences.
