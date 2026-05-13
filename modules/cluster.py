import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pandas as pd
import logging
from collections import defaultdict
import networkx as nx
from tqdm import tqdm
from modules.utilities import kmers, valid_terminal_alignment_tRFi, collapse_cluster_tRFi

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def run_clustering(
    input_matrix_path,
    output_dir=".",
    method="directional",
    min_shared_kmers=5,
    max_seq_per_kmer=1000,
    kmer_size=10,
    min_length_ratio=0.94
):
    logger.info(f"Loading count matrix from: {input_matrix_path}")
    count_matrix = pd.read_table(input_matrix_path)

    count_matrix["tsRNA_seq"] = count_matrix["tsRNA_id"].str.split(":").str[3]
    count_matrix["tsRNA_type"] = count_matrix["tsRNA_id"].str.split(":").str[2]

    sample_cols = [c for c in count_matrix.columns if not c.startswith("tsRNA")]
    count_matrix[sample_cols] = count_matrix[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    count_matrix["rowsum"] = count_matrix[sample_cols].sum(axis=1)

    seq_to_tsRNA_id = {}
    for _, row in count_matrix.iterrows():
        seq = row["tsRNA_seq"]
        tsRNA_id = row["tsRNA_id"]
        seq_to_tsRNA_id[seq] = tsRNA_id

    tR_seq_to_id = dict(zip(count_matrix["tsRNA_seq"], count_matrix["tsRNA_id"]))
    seq_dict_all = dict(zip(count_matrix["tsRNA_seq"], count_matrix["rowsum"]))
    seq_to_type = dict(zip(count_matrix["tsRNA_seq"], count_matrix["tsRNA_type"]))

    tsRNA_dict = count_matrix.set_index("tsRNA_id")[sample_cols + ["rowsum", "tsRNA_seq", "tsRNA_type"]].to_dict(orient="index")

    logger.info("Extracting k-mers...")
    seq_kmers = {s: kmers(s, kmer_size) for s in seq_dict_all}

    kmer_index = defaultdict(set)
    for seq, kmset in seq_kmers.items():
        for kmer in kmset:
            kmer_index[kmer].add(seq)

    logger.info("Finding candidate pairs...")
    pair_counts = defaultdict(int)

    for kmer, seqs in kmer_index.items():
        if len(seqs) > max_seq_per_kmer:
            continue
        seqs_list = sorted(seqs)
        for i in range(len(seqs_list)):
            seqA = seqs_list[i]
            for j in range(i + 1, len(seqs_list)):
                seqB = seqs_list[j]
                pair_counts[(seqA, seqB)] += 1

    logger.info("Validating alignments and building graph...")
    tR_G = nx.Graph()
    tR_G.add_nodes_from((seq, {'count': count}) for seq, count in seq_dict_all.items())

    edges_to_add = []
    for (seqA, seqB), n_shared in tqdm(pair_counts.items(), desc="Aligning candidates"):
        if n_shared < min_shared_kmers:
            continue

        min_len = min(len(seqA), len(seqB))
        max_len = max(len(seqA), len(seqB))

        if max_len == 0 or (min_len / max_len) < min_length_ratio:
            continue

        if valid_terminal_alignment_tRFi(seqA, seqB):
            edges_to_add.append((seqA, seqB))

    tR_G.add_edges_from(edges_to_add)

    logger.info("Clustering tsRNA sequences...")
    tR_clusters = {}
    tR_visited = set()

    for seq in tqdm(sorted(seq_dict_all, key=lambda x: seq_dict_all[x], reverse=True), desc="Clustering tsRNA sequences"):
        if seq in tR_visited:
            continue
        if method == "cluster":
            neighbors = nx.node_connected_component(tR_G, seq)
        elif method == "directional":
            neighbors = {seq}
            stack = [seq]
            counts = seq_dict_all

            while stack:
                current = stack.pop()
                current_count = counts[current]

                for n in tR_G.neighbors(current):
                    if n in tR_visited:
                        continue
                    if current_count >= 2 * counts[n] - 1:
                        neighbors.add(n)
                        stack.append(n)
                        tR_visited.add(n)
        else:
            raise ValueError("Method must be 'cluster' or 'directional'")

        for n in neighbors:
            tR_visited.add(n)

        rep, total_count = collapse_cluster_tRFi(neighbors, tR_G)
        tR_clusters[rep] = {
            "count": total_count,
            "members": neighbors,
            "member_counts": {m: seq_dict_all[m] for m in neighbors},
            "rep_id": tR_seq_to_id.get(rep, rep)
        }

    import os
    os.makedirs(output_dir, exist_ok=True)  # noqa: F811 (os already imported at top)

    seq_to_cluster = {}

    def update_map(clusters, prefix):
        if clusters is None:
            return
        for rep_seq, info in clusters.items():
            cluster_id = info.get("rep_id", f"{prefix}_{rep_seq}")
            for member_seq in info['members']:
                member_id = seq_to_tsRNA_id.get(member_seq, member_seq)
                seq_to_cluster[member_id] = cluster_id

    update_map(tR_clusters, "tR")

    clustered_ids = set(seq_to_cluster.keys())

    for tsRNA_id in count_matrix["tsRNA_id"]:
        if tsRNA_id not in clustered_ids:
            seq_to_cluster[tsRNA_id] = tsRNA_id

    count_matrix['cluster_id'] = count_matrix['tsRNA_id'].map(seq_to_cluster)
    count_matrix['cluster_id'] = count_matrix['cluster_id'].fillna(count_matrix['tsRNA_id'])

    real_sample_cols = [c for c in sample_cols if c != 'rowsum']
    count_matrix['cluster_id'] = count_matrix['cluster_id'].astype(str)

    final_counts = count_matrix.groupby('cluster_id')[real_sample_cols].sum()

    deg_matrix_path = os.path.join(output_dir, "clustered_counts_for_DEG.csv")
    final_counts.to_csv(deg_matrix_path)
    logger.info(f"Clustered count matrix saved to: {deg_matrix_path}")

    summary_data = []

    def add_cluster_summary(cluster_dict, default_type="tsRNA"):
        for rep_seq, info in cluster_dict.items():
            cluster_id = info.get("rep_id", f"{default_type}_{rep_seq}")
            rep_tsRNA_id = cluster_id
            real_type = seq_to_type.get(rep_seq, default_type)
            total_rowsum = info["count"]

            cluster_row = {
                "cluster_id": cluster_id,
                "tsRNA_type": real_type,
                "tsRNA_id": "CLUSTER_TOTAL",
                "rep_tsRNA_id": rep_tsRNA_id,
                "seq": rep_seq,
                "num_members": len(info["members"]),
                "total_count": total_rowsum,
                "is_cluster_total": True
            }
            summary_data.append(cluster_row)

            for member_seq, member_count in info["member_counts"].items():
                member_tsRNA_id = seq_to_tsRNA_id.get(member_seq, member_seq)
                member_real_type = seq_to_type.get(member_seq, real_type)

                member_row = {
                    "cluster_id": cluster_id,
                    "tsRNA_type": member_real_type,
                    "tsRNA_id": member_tsRNA_id,
                    "rep_tsRNA_id": rep_tsRNA_id,
                    "seq": member_seq,
                    "num_members": 1,
                    "total_count": member_count,
                    "is_cluster_total": False
                }
                summary_data.append(member_row)

    add_cluster_summary(tR_clusters, "tsRNA")

    processed_tsRNA_ids = set()
    for row in summary_data:
        if row["tsRNA_id"] != "CLUSTER_TOTAL":
            processed_tsRNA_ids.add(row["tsRNA_id"])

    for _, row in count_matrix.iterrows():
        tsRNA_id = row["tsRNA_id"]
        if tsRNA_id not in processed_tsRNA_ids:
            tsRNA_type = row["tsRNA_type"]
            seq = row["tsRNA_seq"]
            total_count = row["rowsum"]
            cluster_id = tsRNA_id

            cluster_row = {
                "cluster_id": cluster_id, "tsRNA_type": tsRNA_type, "tsRNA_id": "CLUSTER_TOTAL",
                "rep_tsRNA_id": tsRNA_id, "seq": seq, "num_members": 1,
                "total_count": total_count, "is_cluster_total": True
            }
            summary_data.append(cluster_row)

            member_row = {
                "cluster_id": cluster_id, "tsRNA_type": tsRNA_type, "tsRNA_id": tsRNA_id,
                "rep_tsRNA_id": tsRNA_id, "seq": seq, "num_members": 1,
                "total_count": total_count, "is_cluster_total": False
            }
            summary_data.append(member_row)

    summary_simplified = os.path.join(output_dir, "cluster_summary_simplified.csv")

    if summary_data:
        summary_df = pd.DataFrame(summary_data)

        for col in real_sample_cols:
            summary_df[col] = 0

        tsRNA_id_to_sample_counts = {}
        for tsRNA_id, data in tsRNA_dict.items():
            sample_counts = {col: data[col] for col in real_sample_cols if col in data}
            tsRNA_id_to_sample_counts[tsRNA_id] = sample_counts

        for idx, row in summary_df.iterrows():
            if not row["is_cluster_total"]:
                tsRNA_id = row["tsRNA_id"]
                if tsRNA_id in tsRNA_id_to_sample_counts:
                    for sample_col, count in tsRNA_id_to_sample_counts[tsRNA_id].items():
                        summary_df.at[idx, sample_col] = count

        cluster_totals = {}
        unique_cluster_ids = summary_df[summary_df["is_cluster_total"]]["cluster_id"].unique()

        for cluster_id in unique_cluster_ids:
            member_rows = summary_df[(summary_df["cluster_id"] == cluster_id) & (~summary_df["is_cluster_total"])]
            if len(member_rows) > 0:
                sample_totals = {col: member_rows[col].sum() for col in real_sample_cols}
                cluster_totals[cluster_id] = sample_totals

        for idx, row in summary_df[summary_df["is_cluster_total"]].iterrows():
            cluster_id = row["cluster_id"]
            if cluster_id in cluster_totals:
                for col, total in cluster_totals[cluster_id].items():
                    summary_df.at[idx, col] = total

        column_order = ["cluster_id", "tsRNA_type", "tsRNA_id", "rep_tsRNA_id", "seq",
                        "num_members", "total_count", "is_cluster_total"] + real_sample_cols
        summary_df = summary_df[column_order]
        summary_df = summary_df.sort_values(by=["cluster_id", "is_cluster_total"], ascending=[True, False])

        summary_detailed = os.path.join(output_dir, "cluster_summary_detailed.csv")
        summary_df.to_csv(summary_detailed, index=False)
        logger.info(f"Cluster summary saved to: {summary_detailed}")

        cluster_summary_only = summary_df[summary_df["is_cluster_total"]].copy()
        columns_to_drop = ["tsRNA_id", "is_cluster_total"]
        cluster_summary_only = cluster_summary_only.drop(columns=[c for c in columns_to_drop if c in cluster_summary_only.columns])
        cluster_summary_only = cluster_summary_only.rename(columns={"seq": "rep_seq"})

        cluster_summary_only.to_csv(summary_simplified, index=False)
        logger.info(f"Simplified cluster summary saved to: {summary_simplified}")
    else:
        pd.DataFrame().to_csv(summary_simplified)

    return deg_matrix_path, summary_simplified


if __name__ == "__main__":
    import argparse
    import os

    parser = argparse.ArgumentParser(description="tsRNA sequence clustering")
    parser.add_argument("-i", "--input", required=True, help="Count matrix TSV (counts_matrix.tsv)")
    parser.add_argument("-o", "--output-dir", default=".", help="Output directory")
    parser.add_argument("-m", "--method", default="directional", choices=["cluster", "directional"],
                        help="Clustering method (default: directional)")
    parser.add_argument("--min-shared-kmers", type=int, default=5, help="Min shared k-mers for candidates")
    parser.add_argument("--kmer-size", type=int, default=10, help="K-mer size")
    parser.add_argument("--min-length-ratio", type=float, default=0.94, help="Min length ratio for alignment")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    run_clustering(
        input_matrix_path=args.input,
        output_dir=args.output_dir,
        method=args.method,
        min_shared_kmers=args.min_shared_kmers,
        kmer_size=args.kmer_size,
        min_length_ratio=args.min_length_ratio
    )
