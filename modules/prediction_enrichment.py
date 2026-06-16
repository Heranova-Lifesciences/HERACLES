#!/usr/bin/env python3
"""
Integrated pipeline for tsRNA target prediction, filtering, and enrichment analysis.
Usage: python pipeline.py --list list.txt --index target.suf --energy -27 --threads 8 --risearch_path path_to_risearch2 --output_dir ./results
"""

import argparse
import subprocess
import sys
import os
import gzip
import glob
from collections import defaultdict
try:
    import gseapy as gp
    from gseapy.plot import dotplot
except ImportError:
    print("Error: gseapy module not found. Please install: pip install gseapy")
    sys.exit(1)

def count_tsRNAs_in_list(input_list_file):
    """Count the number of tsRNAs in the input list file."""
    count = 0
    with open(input_list_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            count += 1
    return count

# ---------- Step 1: Convert list to FASTA ----------
def list_to_fasta(input_list, fasta_file):
    """Convert a list of tsRNA headers (one per line) to FASTA format."""
    with open(input_list, 'r') as fin, open(fasta_file, 'w') as fout:
        for line in fin:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if ':' in line:
                seq = line.split(':')[-1]
            else:
                seq = line
            fout.write(f'>{line}\n{seq}\n')
    print(f"[Step 1] FASTA file written to {fasta_file}")

# ---------- Step 2: Run RIsearch2 ----------
def run_risearch(risearch_path, query_fasta, index_files, energy, threads, out_dir):
    """Execute RIsearch2 and produce .out.gz files."""
    os.makedirs(out_dir, exist_ok=True)
    abs_query = os.path.abspath(query_fasta)
    for suf in index_files:
        idx_name = os.path.splitext(os.path.basename(suf))[0]
        idx_out = os.path.join(out_dir, idx_name)
        os.makedirs(idx_out, exist_ok=True)
        cmd = [
            risearch_path,
            "-q", abs_query,
            "-i", os.path.abspath(suf),
            "-e", str(energy),
            "-t", str(threads)
        ]
        print(f"[Step 2] RIsearch2 ({idx_name}): {' '.join(cmd)}")
        subprocess.run(cmd, check=True, cwd=idx_out)
    print(f"[Step 2] RIsearch2 finished. Output files are in {out_dir}")


# ---------- Step 4: Decompress and merge ----------
def merge_results(out_dir, merged_file):
    """Decompress all .out.gz files and concatenate them into a single file."""
    gz_files = glob.glob(os.path.join(out_dir, "**", "risearch_*.out.gz"), recursive=True)
    if not gz_files:
        sys.exit("Error: No .out.gz files found. Did RIsearch2 run correctly?")
    print(f"[Step 3] Found {len(gz_files)} .out.gz files. Merging...")
    with open(merged_file, 'w') as fout:
        for gz_path in sorted(gz_files):
            with gzip.open(gz_path, 'rt') as fin:
                for line in fin:
                    fout.write(line)
    print(f"[Step 3] Merged results written to {merged_file}")

# ---------- Step 4: Statistics ----------
def parse_gene(target_col):
    """Extract gene symbol from column like ENST00000411466.7|SRCAP"""
    return target_col.split('|')[-1].strip()

def load_blacklist(index_value):
    """Load blacklist for the given index (CDS or 3UTR). Returns set of genes to exclude."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    blacklist_file = os.path.join(script_dir, "..", "RIsearch2_index", f"{index_value}_black_list.txt")
    blacklist_file = os.path.abspath(blacklist_file)
    if not os.path.exists(blacklist_file):
        print(f"[Info] No blacklist found at {blacklist_file}, skipping filter")
        return set()
    with open(blacklist_file, 'r') as f:
        blacklist = {line.strip() for line in f if line.strip()}
    print(f"[Info] Loaded blacklist: {len(blacklist)} genes excluded")
    return blacklist

def count_genes_per_tsRNA(merged_file):
    """
    Returns:
        all_tsRNAs: set of unique tsRNAs found in the results.
        gene_to_tsRNAs: dict mapping gene -> set of tsRNAs that target it.
        gene_to_total_hits: dict mapping gene -> total binding sites (all RIsearch2 hits).
    """
    gene_to_tsRNAs = defaultdict(set)
    gene_to_total_hits = defaultdict(int)
    all_tsRNAs = set()
    with open(merged_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 4:
                continue
            tsRNA = parts[0]
            gene = parse_gene(parts[3])
            gene_to_tsRNAs[gene].add(tsRNA)
            gene_to_total_hits[gene] += 1
            all_tsRNAs.add(tsRNA)
    return all_tsRNAs, gene_to_tsRNAs, gene_to_total_hits

def parse_fasta_lengths(fa_file):
    """Parse FASTA and return {gene_name: sequence_length}.
    
    Handles headers like '>ENST00000002165.11|FUCA2' -> gene='FUCA2'.
    """
    lengths = {}
    current_gene = None
    seq_len = 0
    with open(fa_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_gene:
                    lengths[current_gene] = seq_len
                header = line[1:]
                if '|' in header:
                    current_gene = header.split('|')[-1].split()[0]
                else:
                    current_gene = header.split()[0]
                seq_len = 0
            else:
                seq_len += len(line)
        if current_gene:
            lengths[current_gene] = seq_len
    print(f"  - Loaded lengths for {len(lengths)} genes")
    return lengths


def output_gene_stats(all_tsRNAs, gene_to_tsRNAs, gene_to_total_hits, stats_file, total_input_count, gene_lengths=None):
    """Write full gene statistics to file, sorted by composite ratio."""

    total = total_input_count
    with open(stats_file, 'w') as fout:
        fout.write("gene\tcount\tfraction\tsites\tlength\tratio\n")
        for gene, tsRNAs in sorted(gene_to_tsRNAs.items(), key=lambda x: len(x[1]), reverse=True):
            cnt = len(tsRNAs)
            frac = cnt / total if total > 0 else 0
            sites = gene_to_total_hits.get(gene, 0)
            length = gene_lengths.get(gene, 0) if gene_lengths else 0
            ratio = (frac * sites / length) if length > 0 else 0
            fout.write(f"{gene}\t{cnt}\t{frac:.4f}\t{sites}\t{length}\t{ratio:.6f}\n")
    print(f"[Step 4] Full gene statistics written to {stats_file}")
    print(f"  - Total input tsRNAs: {total}")
    print(f"  - Unique tsRNAs in results: {len(all_tsRNAs)}")
    return total

# ---------- Step 5: Filter genes ----------
def filter_high_frequency_genes(all_tsRNAs, gene_to_tsRNAs, gene_to_total_hits,
                                threshold=0.5, total_input_count=None, gene_lengths=None):
    """Return list of genes ranked by composite ratio = fraction * sites / length."""

    if total_input_count is None:
        total_input_count = len(all_tsRNAs)

    if total_input_count == 0:
        return []

    cutoff = threshold * total_input_count
    selected = []
    for gene, tsRNAs in gene_to_tsRNAs.items():
        cnt = len(tsRNAs)
        if cnt >= cutoff:
            sites = gene_to_total_hits.get(gene, 0)
            length = gene_lengths.get(gene, 0) if gene_lengths else 0
            frac = cnt / total_input_count
            ratio = (frac * sites / length) if length > 0 else 0
            selected.append((gene, cnt, sites, length, ratio))
    selected.sort(key=lambda x: x[4], reverse=True)
    return selected


def filter_top_percent_genes(all_tsRNAs, gene_to_tsRNAs, gene_to_total_hits,
                             top_percent=5, gene_lengths=None):
    """Return top N% genes sorted by composite ratio = fraction * sites / length."""
    if not gene_to_tsRNAs:
        return []

    total_tsRNAs = len(all_tsRNAs)
    def sort_key(item):
        gene, tsRNAs = item
        cnt = len(tsRNAs)
        sites = gene_to_total_hits.get(gene, 0)
        length = gene_lengths.get(gene, 0) if gene_lengths else 0
        frac = cnt / total_tsRNAs if total_tsRNAs > 0 else 0
        return (frac * sites / length) if length > 0 else 0

    sorted_genes = sorted(gene_to_tsRNAs.items(), key=sort_key, reverse=True)
    n_select = max(1, round(len(sorted_genes) * top_percent / 100))
    selected = []
    for gene, tsRNAs in sorted_genes[:n_select]:
        cnt = len(tsRNAs)
        sites = gene_to_total_hits.get(gene, 0)
        length = gene_lengths.get(gene, 0) if gene_lengths else 0
        frac = cnt / total_tsRNAs if total_tsRNAs > 0 else 0
        ratio = (frac * sites / length) if length > 0 else 0
        selected.append((gene, cnt, sites, length, ratio))
    return selected

# ---------- Step 6: Enrichment Analysis ----------
def run_enrichment_analysis(selected_genes_file, out_dir):
    """
    Perform GO and KEGG enrichment analysis on the selected genes.
    """
    print("[Step 6] Starting enrichment analysis...")
    genes = []
    with open(selected_genes_file, 'r') as f:
        next(f)  
        for line in f:
            line = line.strip()
            if not line:
                continue
            gene = line.split('\t')[0]
            genes.append(gene.upper())
    
    if not genes:
        print("[Step 6] No genes found for enrichment analysis. Skipping...")
        return

    print(f"[Step 6] Performing enrichment analysis on {len(genes)} genes...")
    
    go_libraries = [
        'GO_Biological_Process_2025',
        'GO_Cellular_Component_2025',
        'GO_Molecular_Function_2025'
    ]
    kegg_library = 'KEGG_2021_Human'

    enrich_dir = os.path.join(out_dir, "enrichment")
    os.makedirs(enrich_dir, exist_ok=True)

    def run_and_plot(gene_list, gene_sets, label, out_dir):
        try:
            enr = gp.enrichr(gene_list=gene_list,
                             gene_sets=gene_sets,
                             organism='human',
                             outdir=None,
                             cutoff=1.0)
            results = enr.results
            if results is None or results.empty:
                print(f"[Step 6] No enrichment results found for {label}.")
                return

            old_cols = [c for c in results.columns if c.startswith('Old')]
            if old_cols:
                results = results.drop(columns=old_cols)

            csv_path = os.path.join(out_dir, f"{label}_enrichment.csv")
            results.to_csv(csv_path, index=False)
            print(f"[Step 6] Full {label} enrichment results saved to {csv_path}")

            significant = results[results['P-value'] < 0.05]
            if significant.empty:
                print(f"[Step 6] No significant terms (p < 0.05) for {label}.")
                return

            plot_df = significant.rename(columns={
                'Term': 'Term',
                'P-value': 'P-value',
                'Adjusted P-value': 'Adjusted P-value',
                'Odds Ratio': 'Odds Ratio',
                'Combined Score': 'Combined Score',
                'Genes': 'Genes'
            })

            plot_path = os.path.join(out_dir, f"{label}_dotplot.png")
            
            dynamic_height = max(6, len(plot_df) * 0.4)
            
            dotplot(plot_df,
                    title=f"{label} Enrichment (p < 0.05, n={len(plot_df)})",
                    ofname=plot_path,
                    show_ring=False,
                    top_term=len(plot_df),
                    cutoff=1.0,
                    figsize=(8, dynamic_height))
            print(f"[Step 6] {label} dotplot saved to {plot_path} (showing all {len(plot_df)} terms)")

        except Exception as e:
            print(f"[Step 6] Error during {label} enrichment analysis: {e}")

    for go_lib in go_libraries:
        run_and_plot(genes, go_lib, go_lib, enrich_dir)
    run_and_plot(genes, kegg_library, "KEGG", enrich_dir)
    
    print("[Step 6] Enrichment analysis completed.")

# ---------- Main ----------
def main():
    parser = argparse.ArgumentParser(description="tsRNA target prediction pipeline with enrichment analysis")
    parser.add_argument("--list", required=True, help="Input list file (one tsRNA header per line)")
    parser.add_argument("--index", default="CDS",
                        help="Index to use: CDS, 3UTR (from RIsearch2_index/), or path to custom .suf")
    parser.add_argument("--energy", type=float, default=-29, help="Energy threshold (default: -29)")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads (default: 8)")
    parser.add_argument("--output_dir", default="./pipeline_output", help="Output directory (default: ./pipeline_output)")

    filter_group = parser.add_mutually_exclusive_group()
    filter_group.add_argument("--top-percent", type=float, nargs='?', const=5, default=5,
                              help="Select top N%% of genes by target count (default: 10)")
    filter_group.add_argument("--threshold", type=float, default=argparse.SUPPRESS,
                              help="Fraction of tsRNAs required for gene selection")

    parser.add_argument("--risearch_path",
                        default="RIsearch2",
                        help="Path to RIsearch2 executable (default: RIsearch2, from PATH)")
    args = parser.parse_args()

    out_dir = args.output_dir
    os.makedirs(out_dir, exist_ok=True)

    query_fasta = os.path.join(out_dir, "query.fa")
    merged_results = os.path.join(out_dir, "all_results.txt")
    stats_file = os.path.join(out_dir, "gene_stats.txt")
    selected_file = os.path.join(out_dir, "selected_genes.txt")

    if not os.path.exists(args.list):
        print(f"[Error] Input list file not found: {args.list}")
        sys.exit(1)

    total_input_tsRNAs = count_tsRNAs_in_list(args.list)
    print(f"[Info] Total tsRNAs in input list: {total_input_tsRNAs}")
    if total_input_tsRNAs == 0:
        print("[Error] No tsRNAs found in input list. Exiting.")
        sys.exit(1)

    # Step 1: Generate FASTA
    list_to_fasta(args.list, query_fasta)

    # Step 2: Resolve RIsearch2 index
    gene_lengths = None
    blacklist = set()
    if args.index in ("CDS", "3UTR"):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        fa_file = os.path.join(script_dir, "..", "RIsearch2_index", f"GRCh38.{args.index}.fa")
        fa_file = os.path.abspath(fa_file)
        suf_file = fa_file.replace('.fa', '.suf')
        if not os.path.exists(suf_file):
            cmd = [args.risearch_path, "-c", fa_file, "-o", suf_file]
            print(f"[Init] Building index: {' '.join(cmd)}")
            subprocess.run(cmd, check=True)
        gene_lengths = parse_fasta_lengths(fa_file)
        blacklist = load_blacklist(args.index)
        index_files = [suf_file]
    else:
        index_files = [args.index]

    # Step 3: Run RIsearch2
    run_risearch(args.risearch_path, query_fasta, index_files, args.energy, args.threads, out_dir)

    # Step 4: Merge results
    merge_results(out_dir, merged_results)

    # Step 5: Statistics
    all_tsRNAs, gene_to_tsRNAs, gene_to_total_hits = count_genes_per_tsRNA(merged_results)

    total_tsRNAs = output_gene_stats(all_tsRNAs, gene_to_tsRNAs, gene_to_total_hits, stats_file,
                                     total_input_tsRNAs, gene_lengths)

    # Step 5: Filter genes
    if hasattr(args, 'top_percent'):
        selected = filter_top_percent_genes(all_tsRNAs, gene_to_tsRNAs, gene_to_total_hits,
                                            args.top_percent, gene_lengths)
        filter_desc = f"top {args.top_percent}% of genes by composite score"
    else:
        selected = filter_high_frequency_genes(all_tsRNAs, gene_to_tsRNAs, gene_to_total_hits,
                                               args.threshold, total_input_tsRNAs, gene_lengths)
        filter_desc = f"> {args.threshold*100:.0f}% of input tsRNAs"

    # Remove blacklisted genes from selected set
    if blacklist:
        before = len(selected)
        selected = [x for x in selected if x[0] not in blacklist]
        filtered = before - len(selected)
        print(f"[Filter] Removed {filtered} blacklisted genes from selected set, {len(selected)} remaining")

    with open(selected_file, 'w') as f:
        f.write("gene\tcount\tfraction\tsites\tlength\tratio\n")
        total = total_input_tsRNAs
        for gene, cnt, sites, length, ratio in selected:
            frac = cnt / total if total > 0 else 0
            f.write(f"{gene}\t{cnt}\t{frac:.4f}\t{sites}\t{length}\t{ratio:.6f}\n")

    print("\n=== Initial Pipeline Completed ===")
    print(f"Total input tsRNAs: {total_input_tsRNAs}")
    print(f"Unique tsRNAs in results: {len(all_tsRNAs)}")
    print(f"Genes selected ({filter_desc}): {len(selected)}")
    if selected:
        print("Top selected genes:")
        for gene, cnt, sites, length, ratio in selected[:10]:
            frac = cnt / total_input_tsRNAs if total_input_tsRNAs > 0 else 0
            print(f"  {gene}: {cnt} ({frac:.1%}) sites={sites} len={length} score={ratio:.6f}")
    
    if selected:
        run_enrichment_analysis(selected_file, out_dir)
    else:
        print("\nNo genes selected for enrichment analysis.")

    print(f"\nAll output files are in: {out_dir}")

if __name__ == "__main__":
    main()
