# HERACLES container

A pre-built image and a single environment file, so that none of the external tools — Trim Galore, cutadapt, Bowtie, RIsearch2 — has to be compiled or configured by hand.

## What you need first

Docker (or Singularity/Apptainer on an HPC cluster), and your FASTQ files.

Nothing else. The reference FASTA, the Bowtie tRNA index and the blacklists are already in the repository; the only thing built on your side is a suffix-array index, and that takes about twenty seconds.

## 1. Get the code

```bash
git clone https://github.com/Heranova-Lifesciences/HERACLES.git
cd HERACLES
```

## 2. Build the image

```bash
docker build -f docker/Dockerfile -t heracles:1.0 .
```

Roughly three minutes cold, under a minute once the conda packages are in the build cache. The build ends with a self-check: the four external tools and the eight Python libraries are imported inside the image, and the build fails if any of them is missing.

## 3. Lay out your data

```bash
mkdir -p data/{fastq,index,out}
# copy your FASTQ files into data/fastq/
# write data/metadata.tsv: sample<TAB>condition, no header
```

```
data/
├── fastq/            # raw or already-trimmed FASTQ
├── index/            # created in step 4
├── metadata.tsv      # sample<TAB>condition, one line per sample
└── out/              # results
```

## 4. Build the target index (once, only needed for `predict`)

The repository ships `RIsearch2_index/GRCh38.CDS.fa` and `GRCh38.3UTR.fa` but not the `.suf` suffix arrays: the CDS array alone is 541 MB, too large for git. Build them into your data directory; the step is skipped if the files already exist.

```bash
docker run --rm -v "$PWD/data:/data" --entrypoint bash heracles:1.0 \
  /opt/heracles/docker/build_index.sh /data/index
```

About twenty seconds for CDS. The results live in `data/index/` and are reused on every later run.

## 5. Run

```bash
docker run --rm -v "$PWD/data:/data" heracles:1.0 \
  --fastq-dir /data/fastq \
  --metadata /data/metadata.tsv \
  --contrast Tumor Control \
  --index-dir /opt/heracles/tRNA_index \
  --predict-index /data/index/GRCh38.CDS.suf \
  --output-dir /data/out
```

`tRNA_index` comes from the image itself, so it needs no preparation. Everything written goes to `data/out/` on your host.

The default stages end at DESeq2, so `--predict-index` is not used by the command above. Add the prediction and primer-extension stages when you want them, and only then is the index from step 4 required:

```bash
docker run --rm -v "$PWD/data:/data" heracles:1.0 \
  --stages full \
  --fastq-dir /data/fastq \
  --metadata /data/metadata.tsv \
  --contrast Tumor Control \
  --index-dir /opt/heracles/tRNA_index \
  --predict-index /data/index/GRCh38.CDS.suf \
  --output-dir /data/out
```

To run a subset of stages:

```bash
docker run --rm -v "$PWD/data:/data" heracles:1.0 \
  --stages QC,Collapse --fastq-dir /data/fastq --output-dir /data/out
```

List all options with:

```bash
docker run --rm heracles:1.0 --help
```

## Without Docker

```bash
conda env create -f docker/environment.yml
conda activate heracles
bash docker/build_index.sh ./RIsearch2_index
python main.py --help
```

The same environment the image uses, from the same file.

## Singularity / Apptainer on an HPC cluster

```bash
singularity build heracles_1.0.sif docker-daemon://heracles:1.0
```

Singularity does not execute `ENTRYPOINT`, so name the program explicitly:

```bash
singularity exec --bind "$PWD/data:/data" heracles_1.0.sif \
  python /opt/heracles/main.py \
    --fastq-dir /data/fastq \
    --metadata /data/metadata.tsv \
    --contrast Tumor Control \
    --index-dir /opt/heracles/tRNA_index \
    --predict-index /data/index/GRCh38.CDS.suf \
    --output-dir /data/out
```

## What is in the image

| Component | Version | Channel |
|---|---|---|
| Python | 3.10 | conda-forge |
| Bowtie | 1.3.1 | bioconda |
| cutadapt | 5.2 | bioconda |
| Trim Galore | 0.6.10 | bioconda |
| RIsearch | 2.1 | **rthtools** (public channel on anaconda.org) |
| pydeseq2 | 0.5.2 | PyPI |
| gseapy | 1.2.1 | PyPI |

RIsearch 2.1 is not on bioconda; it comes from `rthtools`, which conda resolves without any compilation step.

Also inside: the code, `tRNA_index/` (17 MB, so annotation works immediately) and `docker/smoke_test.sh`.

Not inside, by design: any `.suf` index. Those are half a gigabyte each and belong in a mounted volume.

## Checking an installation

```bash
docker run --rm --entrypoint bash heracles:1.0 -c 'bash /opt/heracles/docker/smoke_test.sh'
```

Eighteen checks: four external tools, eight Python libraries, the entry point, and the reference files. Pass a data directory as the first argument to include the reference checks.

## Verification status

| Item | Status |
| --- | --- |
| Every conda package and version exists in its channel | **checked** against the anaconda.org API |
| `environment.yml` resolves | **tested** with `conda create --dry-run` |
| Environment created from the file | **tested** on a clean Linux host, about 8 minutes |
| `smoke_test.sh` in that environment | **18 of 18 passed** |
| `docker build` | **built in CI**, 158 s ([run 40266387452](https://github.com/Heranova-Lifesciences/HERACLES/actions/runs/40266387452)) |
| Smoke test inside the built image | **passed** in the same run |
| Suffix-array index built by `build_index.sh` | **tested**, 18 s for CDS, 541 MB output; 25 s inside CI |
| Full pipeline inside the image | **exercised in CI** end to end: QC, collapse, annotation, cluster and DESeq2 |
| Full pipeline on real libraries | **exercised in CI** on four real libraries truncated to 10,000 reads each: 12 s, 1.9 MB of output ([run 40279774012](https://github.com/Heranova-Lifesciences/HERACLES/actions/runs/40279774012)) |

The four libraries in `tests/data/` are real reads truncated to 10,000 each so they fit in the repository. The synthetic fallback set is slower (686 s against 12 s) because synthetic reads almost never collapse, which inflates the fragment count.

## Five pitfalls, already handled here

0. **`docker run image bash script.sh` does not do what it looks like.** The image has `ENTRYPOINT ["python", "main.py"]`, so the arguments are appended to that command and you get `python main.py bash script.sh`. Use `--entrypoint bash` for anything that is not the pipeline itself. The end-to-end workflow caught this in the documentation.

1. **`RIsearch2 --help` exits with status 1, not 0.** A presence check written as `RIsearch2 --help || exit 1` fails the build outright. The check matches the banner text `RIsearch v2` instead.
2. **micromamba 1.5.x runs the `pip:` section before creating the environment**, failing with `Cannot activate, prefix does not exist`. The image therefore uses conda, which orders the two correctly. This was observed in CI, not guessed.
3. **`bash -l` re-sources `/etc/profile` and discards the `PATH` set by the image**, so `bowtie` and friends appear to be missing. Run the smoke test with `-c`, not `-lc`.
4. **`smoke_test.sh` must keep LF line endings.** Edited on Windows and uploaded with CRLF it fails with `set: -: invalid option`. `.gitattributes` pins `*.sh` to LF.

## If the build fails

| Symptom | Cause | Fix |
| --- | --- | --- |
| Pulling `continuumio/miniconda3` times out | network or registry mirror | configure a mirror, or pre-pull the image on the target host |
| `risearch` cannot be solved | `rthtools` not active | keep `rthtools` in `channels`, alongside `conda-forge` |
| `COPY` cannot find `docker/environment.yml` | wrong build context | run `docker build -f docker/Dockerfile .` from the repository root |
