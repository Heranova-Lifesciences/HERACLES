#!/usr/bin/env python3
"""
Generate a tiny, deterministic FASTQ set for the end-to-end workflow.

Reads are drawn from the tRNA reference that ships with the repository, with the
Illumina small-RNA 3' adapter appended, so the QC, collapse, annotation and
cluster stages all have something real to work on.

Group A and group B draw from disjoint halves of the tRNA set, which gives
DESeq2 a difference to find. The data are synthetic: this checks that the
pipeline runs end to end, not that any biological result is meaningful.

Usage:
    python3 tests/make_test_data.py --fastq-dir run/fastq --metadata run/metadata.tsv
"""
import argparse
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))
TRNA_FASTA = os.path.join(HERE, "..", "tRNA_index", "hg38-tRNA.fa")
ADAPTER = "TGGAATTCTCGGGTGCCAAGG"          # Illumina TruSeq small RNA 3' adapter
SAMPLES = [("test_A_1", "A", 1), ("test_A_2", "A", 2),
           ("test_B_1", "B", 3), ("test_B_2", "B", 4)]


def load_trna():
    seqs, seen_header, buf = [], False, []
    with open(TRNA_FASTA) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if buf:
                    seqs.append("".join(buf))
                seen_header, buf = True, []
            elif seen_header:
                buf.append(line.upper())
    if buf:
        seqs.append("".join(buf))
    return [s for s in seqs if len(s) >= 40]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fastq-dir", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--n", type=int, default=20000)
    args = ap.parse_args()

    os.makedirs(args.fastq_dir, exist_ok=True)
    trna = load_trna()
    print(f"tRNA reference: {len(trna)} sequences")

    half = len(trna) // 2
    meta = []
    for name, group, seed in SAMPLES:
        rng = random.Random(seed)
        pool = trna[:half] if group == "A" else trna[half:]
        path = os.path.join(args.fastq_dir, f"{name}.fq")
        with open(path, "w") as fh:
            for i in range(args.n):
                src = rng.choice(pool)
                length = rng.randint(18, 36)
                start = rng.randint(0, max(0, len(src) - length))
                read = src[start:start + length]
                fh.write(f"@{name}_{i}\n{read}{ADAPTER}\n+\n{'I' * (len(read) + len(ADAPTER))}\n")
        meta.append(f"{name}\t{group}")
        print(f"  {name}: {args.n} reads from {len(pool)} tRNAs")

    os.makedirs(os.path.dirname(os.path.abspath(args.metadata)), exist_ok=True)
    with open(args.metadata, "w") as fh:
        fh.write("\n".join(meta) + "\n")
    print(f"metadata written to {args.metadata}")


if __name__ == "__main__":
    main()
