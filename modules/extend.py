"""
tsRNA Extension & Primer Search Module
=======================================
1. extend_tsrna_sequences()   — extend tsRNA sequences by ±N nt using tRNA reference FASTA
2. run_primer_search()         — search count matrix for sequences similar to extended tsRNAs
3. run_extend_and_primer()     — combined: extend + search in one call
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import csv
import logging
from pathlib import Path
from typing import Optional

import pandas as pd

from modules.utilities import adaptive_subset, align_to_ref

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# ============================================================
#  Extension
# ============================================================

def extend_tsrna_sequences(
    fa_file: str,
    list_file: str,
    count_csv: str,
    out_file: str,
    extend_by: int = 5
):
    fasta_dict = {}
    with open(fa_file, 'r') as f:
        header = None
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                header = line[1:]
                fasta_dict[header] = ""
            elif header:
                fasta_dict[header] += line.lower()

    count_dict = {}
    if count_csv and Path(count_csv).exists():
        with open(count_csv, 'r') as f:
            reader = csv.reader(f)
            header = next(reader)
            try:
                cluster_id_idx = header.index('cluster_id')
            except ValueError:
                cluster_id_idx = -1
            try:
                total_count_idx = header.index('total_count')
            except ValueError:
                total_count_idx = -1
            if cluster_id_idx < 0 or total_count_idx < 0:
                logger.warning("CSV missing expected columns; skipping count lookup")

            for row in reader:
                if cluster_id_idx >= 0 and len(row) > max(cluster_id_idx, total_count_idx):
                    count_dict[row[cluster_id_idx]] = row[total_count_idx]

    extended_count = 0
    with open(list_file, 'r') as f_list, open(out_file, 'w', newline='') as f_out:
        writer = csv.writer(f_out)
        writer.writerow(['tsRNA_ID', 'tsRNA_extended', 'count'])

        for line in f_list:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split(':')
            if len(parts) < 3:
                logger.warning(f"Skipping malformed line: {line}")
                continue

            trna_id = parts[0]
            pos_str = parts[1]
            tsrna_type = parts[2]

            start_str, end_str = pos_str.split('-')
            start = int(start_str)
            end = int(end_str)

            full_seq = fasta_dict.get(trna_id, "")
            seq_len = len(full_seq)

            if seq_len == 0:
                logger.warning(f"tRNA {trna_id} not found in FASTA, skipping")
                continue

            new_start = max(1, start - extend_by)
            new_end = min(seq_len, end + extend_by)

            new_seq = full_seq[new_start - 1: new_end]

            tsrna_extended = f"{trna_id}:{new_start}-{new_end}:{tsrna_type}:{new_seq}"

            writer.writerow([line, tsrna_extended, count_dict.get(line, "0")])
            extended_count += 1

    logger.info(f"Extended {extended_count} tsRNAs → {out_file}")
    return out_file


# ============================================================
#  Primer / Similar Sequence Search
# ============================================================

def run_primer_search(extend_csv, count_matrix_path, output_file):
    logger.info(f"Loading extended DEGs from: {extend_csv}")
    degs_matrix = pd.read_csv(extend_csv)
    degs_matrix["tsRNA_orig_ID"] = degs_matrix["tsRNA_ID"].str.split(":").str[0]
    degs_matrix["tsRNA_orig_start"] = degs_matrix["tsRNA_ID"].str.split(":").str[1].str.split("-").str[0]
    degs_matrix["tsRNA_orig_end"] = degs_matrix["tsRNA_ID"].str.split(":").str[1].str.split("-").str[1]
    degs_matrix["tsRNA_orig_seq"] = degs_matrix["tsRNA_ID"].str.split(":").str[3]
    degs_matrix["tsRNA_ext_start"] = degs_matrix["tsRNA_extended"].str.split(":").str[1].str.split("-").str[0]
    degs_matrix["tsRNA_ext_end"] = degs_matrix["tsRNA_extended"].str.split(":").str[1].str.split("-").str[1]
    degs_matrix["tsRNA_ext_seq"] = degs_matrix["tsRNA_extended"].str.split(":").str[3]

    n_before = len(degs_matrix)
    degs_matrix = degs_matrix.dropna(subset=["tsRNA_orig_seq", "tsRNA_ext_seq"])
    n_dropped = n_before - len(degs_matrix)
    if n_dropped > 0:
        logger.warning(f"Dropped {n_dropped} rows with malformed tsRNA IDs (missing sequence field)")

    logger.info(f"Loading count matrix from: {count_matrix_path}")
    count_matrix = pd.read_table(count_matrix_path)
    count_matrix["tsRNA_seq"] = count_matrix["tsRNA_id"].str.split(":").str[3]
    count_matrix["tsRNA_type"] = count_matrix["tsRNA_id"].str.split(":").str[2]
    count_matrix = count_matrix.dropna(subset=["tsRNA_seq"])

    sample_cols = [c for c in count_matrix.columns if not c.startswith("tsRNA")]
    count_matrix[sample_cols] = count_matrix[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    count_matrix["rowsum"] = count_matrix[sample_cols].sum(axis=1)

    seq_to_ids = count_matrix.groupby("tsRNA_seq")["tsRNA_id"].apply(list).to_dict()
    seq_dict_all = dict(zip(count_matrix["tsRNA_seq"], count_matrix["rowsum"]))

    all_results = []

    for idx, row in degs_matrix.iterrows():
        tsRNA_orig = row["tsRNA_orig_seq"]
        if pd.isna(tsRNA_orig) or not tsRNA_orig:
            continue
        seq_kmer = adaptive_subset(tsRNA_orig)
        query = row["tsRNA_ext_seq"]

        matches = {seq: count for seq, count in seq_dict_all.items() if seq_kmer in seq}
        results = align_to_ref(query, matches, tsRNA_orig)
        all_results.extend(results)

    df = pd.DataFrame(all_results)

    if df.empty:
        logger.warning("No similar sequences found; writing empty output.")
        pd.DataFrame(columns=["original_tsRNA", "query_seq", "aligned_tsRNA",
                              "aligned_tsRNA_similar_cleaned", "matched_tsRNA_ids",
                              "target_count", "similarity", "alignment_score"]).to_csv(
            output_file, index=False)
        return output_file

    target_col = "aligned_tsRNA_similar_cleaned"

    def get_matched_ids(seq):
        ids = seq_to_ids.get(seq, [])
        return ",".join(ids) if ids else ""

    df["matched_tsRNA_ids"] = df[target_col].apply(get_matched_ids)

    cols = df.columns.tolist()
    if target_col in cols:
        idx = cols.index(target_col)
        cols.insert(idx + 1, cols.pop(cols.index("matched_tsRNA_ids")))
        df = df[cols]

    df_display = df.copy()
    df_display.loc[df_display.duplicated("original_tsRNA"), "original_tsRNA"] = ""
    df_display.loc[df_display.duplicated(target_col), target_col] = ""

    df_display.to_csv(output_file, index=False)
    logger.info(f"Result saved to {output_file}")
    return output_file


# ============================================================
#  Combined: Extend → Primer
# ============================================================

def run_extend_and_primer(
    fa_file: str,
    list_file: str,
    count_csv: str,
    count_matrix_path: str,
    extend_dir: str,
    extend_by: int = 5
):
    import os
    os.makedirs(extend_dir, exist_ok=True)

    extend_csv = os.path.join(extend_dir, "output_result.csv")
    primer_csv = os.path.join(extend_dir, "clustered_tsRNA_similar_sequences.csv")

    extend_tsrna_sequences(
        fa_file=fa_file,
        list_file=list_file,
        count_csv=count_csv,
        out_file=extend_csv,
        extend_by=extend_by
    )

    run_primer_search(
        extend_csv=extend_csv,
        count_matrix_path=count_matrix_path,
        output_file=primer_csv
    )

    return extend_csv, primer_csv


# ============================================================
#  CLI
# ============================================================

if __name__ == "__main__":
    import argparse
    import os

    parser = argparse.ArgumentParser(
        description="tsRNA Extension & Primer Search Module"
    )
    sub = parser.add_subparsers(dest="command")

    # Sub-command: extend
    p_ext = sub.add_parser("extend", help="Extend tsRNA sequences by ±N nt")
    p_ext.add_argument("--fa", required=True, help="tRNA reference FASTA")
    p_ext.add_argument("--list", required=True, help="tsRNA ID list (one per line)")
    p_ext.add_argument("--csv", help="Cluster summary CSV for count lookup")
    p_ext.add_argument("--out", default="output_result.csv", help="Output CSV")
    p_ext.add_argument("--extend-by", type=int, default=5, help="nt to extend per side")

    # Sub-command: primer
    p_pri = sub.add_parser("primer", help="Search similar sequences for extended tsRNAs")
    p_pri.add_argument("--extend", required=True, help="Extended tsRNA CSV (output_result.csv)")
    p_pri.add_argument("--count-matrix", required=True, help="Count matrix TSV")
    p_pri.add_argument("--out", default="clustered_tsRNA_similar_sequences.csv", help="Output CSV")

    # Sub-command: all (extend + primer)
    p_all = sub.add_parser("all", help="Run extend + primer together")
    p_all.add_argument("--fa", required=True, help="tRNA reference FASTA")
    p_all.add_argument("--list", required=True, help="tsRNA ID list")
    p_all.add_argument("--csv", help="Cluster summary CSV for count lookup")
    p_all.add_argument("--count-matrix", required=True, help="Count matrix TSV")
    p_all.add_argument("--out-dir", default="extend_primer_results", help="Output directory")
    p_all.add_argument("--extend-by", type=int, default=5, help="nt to extend per side")

    args = parser.parse_args()

    if args.command == "extend":
        extend_tsrna_sequences(
            fa_file=args.fa, list_file=args.list, count_csv=args.csv,
            out_file=args.out, extend_by=args.extend_by
        )
    elif args.command == "primer":
        run_primer_search(args.extend, args.count_matrix, args.out)
    elif args.command == "all":
        run_extend_and_primer(
            fa_file=args.fa, list_file=args.list, count_csv=args.csv,
            count_matrix_path=args.count_matrix, extend_dir=args.out_dir,
            extend_by=args.extend_by
        )
    else:
        parser.print_help()
