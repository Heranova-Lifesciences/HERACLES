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

## Bundled resources

| Path | Contents |
|------|----------|
| `resources/gmt/` | The four Enrichr libraries used in this study as plain GMT files (KEGG_2021_Human, GO_Biological_Process_2025, GO_Cellular_Component_2025, GO_Molecular_Function_2025), downloaded from Enrichr and shipped so the `predict` stage runs offline |
| `RIsearch2_index/*_black_list.txt` | The curated CDS and 3′ UTR blacklists applied after target prediction |
| `RIsearch2_index/*.fa` | The CDS and 3′ UTR reference sequences from which the RIsearch2 suffix arrays are built; the arrays themselves are not shipped (see `docker/build_index.sh`) |

## Recent changes (November 2026)

- **Optional QC filter** `--drop-unsupported` / `--cycle-length`; off by default, so existing runs are unchanged. See the QC Module section of `HERACLES_User_Manual.md` for usage and output details.
- **`--block-cross-parent`** stops two fragments from being joined when their parental tRNAs differ; off by default. On the in-house matrix it raises the cluster count from 16,931 to 17,559. See the Clustering Module section of the manual.
- The clustering abundance rule is exposed as `abund_ratio` in `modules/cluster.py`; the default of `2.0` reproduces the previous hard-coded behaviour.
- Local sequence alignment moved from the deprecated `Bio.pairwise2` to `Bio.Align.PairwiseAligner` with identical scoring (match 2, mismatch −1, gap open −5, gap extend −0.1); requires `biopython>=1.80`.
- **GO and KEGG libraries ship with the repository** as GMT files under `resources/gmt/` (KEGG_2021_Human, GO_Biological_Process_2025, GO_Cellular_Component_2025, GO_Molecular_Function_2025), so the `predict` stage runs offline by default. `HERACLES_GMT_DIR` points at an alternative directory; if a library is found in neither place the online Enrichr lookup is used as a fallback. Over-representation is evaluated against a fixed background of 20,000 genes (`HERACLES_ENRICHR_BACKGROUND` overrides), which reproduces the online Enrichr p-values.

## More

See `HERACLES_User_Manual.md` for full documentation including per-module usage, output structure, and FAQ.

## License

PolyForm Noncommercial License 1.0.0 — see [LICENSE](LICENSE).
Copyright (c) 2026 Heranova Lifesciences.

This licence permits any noncommercial purpose, including personal research,
study and use by educational, public-research, public-health and government
institutions. Commercial use requires a separate licence from the copyright
holder. The software is provided as is, without warranty of any kind.

The required notice to be carried with copies and modified versions is:

```
Required Notice: Copyright (c) 2026 Heranova Lifesciences
```

The full licence text is at
https://polyformproject.org/licenses/noncommercial/1.0.0.
