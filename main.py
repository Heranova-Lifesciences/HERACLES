#!/usr/bin/env python3
"""
HERACLES - Integrated tsRNA Analysis Pipeline
=============================================
Stages: qc → collapse → annotation → merge → cluster → deseq2 → extend → predict

Usage:
    python main.py --fastq-list samples.txt --index-dir ./index --metadata metadata.tsv
                   --contrast Treat Control --stages cluster,deseq2

    python main.py --fastq-list samples.txt --index-dir ./index --metadata metadata.tsv
                   --contrast Treat Control --stages qc,collapse,annotation,cluster,deseq2,extend,predict
"""

import argparse
import logging
import subprocess
import sys
import os
from pathlib import Path
from typing import List, Dict

import pandas as pd

from modules.qc import QCProcessor, prepare_input_list
from modules.collapse import FastaCollapser
from modules.tsRNA_annotation import SimpleTsRNAAnalyzer
from modules.cluster import run_clustering
from modules.extend import extend_tsrna_sequences, run_primer_search
from modules.run_deseq2 import plot_volcano, plot_heatmap

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger("HERACLES")

ALL_STAGES = ["qc", "collapse", "annotation", "cluster", "deseq2", "extend", "predict"]
DEFAULT_STAGES = "qc,collapse,annotation,cluster,deseq2"


def merge_counts(count_files: Dict[str, str], output_path: str) -> str:
    logger.info(f"Merging {len(count_files)} sample count files...")
    merged = None
    for sample_name, count_file in count_files.items():
        df = pd.read_csv(count_file, sep='\t', index_col=0)
        if 'count' in df.columns:
            df = df[['count']]
        elif len(df.columns) == 1:
            pass
        else:
            df = df.iloc[:, [0]]
        df.columns = [sample_name]
        if merged is None:
            merged = df
        else:
            merged = merged.join(df, how='outer')

    merged = merged.fillna(0).astype(int)
    merged = merged.sort_index()
    merged.index.name = 'tsRNA_id'
    merged.to_csv(output_path, sep='\t')
    logger.info(f"Merged count matrix ({merged.shape[0]} tsRNAs x {merged.shape[1]} samples) → {output_path}")
    return output_path


def split_clustered_matrix(clustered_csv, output_dir, sample_cols=None):
    os.makedirs(output_dir, exist_ok=True)
    sep = '\t' if clustered_csv.endswith('.tsv') else ','
    df = pd.read_csv(clustered_csv, sep=sep, index_col=0)
    if sample_cols is None:
        sample_cols = [c for c in df.columns]
    sample_files = {}
    for col in sample_cols:
        if col in df.columns:
            out_file = os.path.join(output_dir, f"{col}_counts.tsv")
            sub = df[[col]].copy()
            sub.index.name = 'tsRNA_id'
            sub.to_csv(out_file, sep='\t')
            sample_files[col] = out_file
    return sample_files


def generate_significant_list(deseq2_results_tsv, list_output, pvalue_thresh=0.05):
    df = pd.read_csv(deseq2_results_tsv, sep='\t', index_col=0)
    sig = df[df['pvalue'] < pvalue_thresh]
    with open(list_output, 'w') as f:
        for tsrna_id in sig.index:
            f.write(f"{tsrna_id}\n")
    logger.info(f"Wrote {len(sig)} significant tsRNAs (pvalue < {pvalue_thresh}) → {list_output}")
    return list_output


def main():
    parser = argparse.ArgumentParser(
        description="HERACLES - Integrated tsRNA Analysis Pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # --- Pipeline ---
    parser.add_argument('--stages', default=DEFAULT_STAGES,
                        help=f'Comma-separated stages to run. '
                             f'Available: {",".join(ALL_STAGES)}. '
                             f'Aliases: "full" (all stages), "full_without_cluster" (skip cluster). '
                             f'Default: {DEFAULT_STAGES}')
    parser.add_argument('--output-dir', default='HERACLES_output', help='Root output directory (all modules)')

    # --- QC module ---
    parser.add_argument('--fastq-list', help='File listing FASTQ paths, one per line (QC, Collapse modules)')
    parser.add_argument('--fastq-dir', help='Directory containing FASTQ files, alternative to --fastq-list (QC module)')
    parser.add_argument('--qc-threads', type=int, default=4, help='Number of threads for Trim Galore (QC module)')
    parser.add_argument('--qc-quality', type=int, default=20, help='Phred quality trimming threshold (QC module)')
    parser.add_argument('--qc-length', type=int, default=18, help='Minimum read length after trimming (QC module)')
    parser.add_argument('--trim-galore-path', default='trim_galore', help='Path to trim_galore executable (QC module)')
    parser.add_argument('--adapter', default='', help='Adapter sequence for trimming, auto-detected if empty (QC module)')

    # --- Collapse module ---
    parser.add_argument('--min-count', type=int, default=1, help='Minimum read count for sequence retention (Collapse module)')
    parser.add_argument('--collapsed-dir', help='Directory with pre-existing collapsed FASTAs, for resuming skipped collapse (Collapse module)')

    # --- Annotation module ---
    parser.add_argument('--index-dir', default='tRNA_index', help='tRNA Bowtie index directory (Annotation module)')
    parser.add_argument('--min-len', type=int, default=18, help='Min tsRNA fragment length (Annotation module)')
    parser.add_argument('--max-len', type=int, default=50, help='Max tsRNA fragment length (Annotation module)')
    parser.add_argument('--mismatch', type=int, default=0, help='Bowtie mismatches allowed (Annotation module)')
    parser.add_argument('--threads', type=int, default=4, help='Threads for Bowtie (Annotation, QC modules)')
    parser.add_argument('--bowtie-path', default='bowtie', help='Path to bowtie executable (Annotation module)')
    parser.add_argument('--keep-temp', action='store_true', help='Keep temporary intermediate files (Annotation, QC modules)')

    # --- Cluster module ---
    parser.add_argument('--cluster-method', default='directional', choices=['cluster', 'directional'],
                        help='Clustering method: directional (stricter) or cluster (Cluster module)')

    # --- DESeq2 module ---
    parser.add_argument('--metadata', help='Metadata: <sample_name><TAB><condition>, no header (DESeq2 module)')
    parser.add_argument('--contrast', nargs=2, metavar=('TREATMENT', 'CONTROL'),
                        help='DESeq2 contrast: e.g. --contrast Treat Control (DESeq2 module)')
    parser.add_argument('--design', default='condition', help='Design formula factor name (DESeq2 module)')
    parser.add_argument('--min-count-deseq2', type=int, default=10, help='Min total count for DESeq2 filtering (DESeq2 module)')
    parser.add_argument('--top-n-heatmap', type=int, default=50, help='Top N genes shown in heatmap (DESeq2 module)')
    parser.add_argument('--pvalue-thresh', type=float, default=0.05,
                        help='P-value threshold for significant DEGs (DESeq2 module)')
    parser.add_argument('--output-normalized', nargs='?', const='counts_matrix_normalized.tsv', default='',
                        help='Output DESeq2 normalized count matrix to this TSV file (DESeq2 module)')

    # --- Extend module ---
    parser.add_argument('--tRNA-fasta', help='tRNA reference FASTA. Auto-detected from index-dir. (Extend module)')
    parser.add_argument('--extend-by', type=int, default=5, help='Nucleotides to extend on each side (Extend module)')

    # --- Predict module ---
    parser.add_argument('--predict-index', default='CDS', help='RIsearch2 index: CDS, 3UTR, or path to .suf (Predict module)')
    parser.add_argument('--energy', type=float, default=-27, help='Free energy threshold for RIsearch2, kcal/mol (Predict module)')
    parser.add_argument('--threshold', type=float, default=0.5, help='Gene frequency threshold for enrichment (Predict module)')
    parser.add_argument('--top-percent', type=float,
                        help='Select top N%% of genes by target count, alternative to --threshold (Predict module)')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    # Resolve relative defaults to software directory
    if not os.path.isabs(args.index_dir):
        args.index_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), args.index_dir)

    # Parse stages into a set for fast lookup
    stages_input = set(s.strip().lower() for s in args.stages.split(',') if s.strip())

    if 'full' in stages_input:
        stages = set(ALL_STAGES)
    elif 'full_without_cluster' in stages_input:
        stages = set(ALL_STAGES) - {'cluster'}
    else:
        invalid = stages_input - set(ALL_STAGES)
        if invalid:
            logger.error(f"Invalid stage(s): {', '.join(sorted(invalid))}. "
                         f"Valid stages: {', '.join(ALL_STAGES)}. "
                         f"Aliases: full, full_without_cluster")
            sys.exit(1)
        stages = stages_input

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Log which stages will run
    enabled = [s for s in ALL_STAGES if s in stages]
    logger.info(f"Stages to run: {', '.join(enabled)}")

    # --- Validate FASTQ list early ---
    if args.fastq_list and not Path(args.fastq_list).exists():
        logger.error(f"FASTQ list file not found: {args.fastq_list}")
        sys.exit(1)

    # ===================================================================
    #  STAGE: qc
    # ===================================================================
    trimmed_files = {}

    if 'qc' in stages:
        logger.info("=" * 50)
        logger.info("  STAGE: Quality Control (Trim Galore + MultiQC)")
        logger.info("=" * 50)

        qc_output = output_dir / "qc_results"
        qc = QCProcessor(
            output_dir=str(qc_output),
            trim_galore_path=args.trim_galore_path,
            threads=args.qc_threads,
            quality=args.qc_quality,
            length=args.qc_length,
            adapter=args.adapter
        )

        if args.fastq_list and os.path.isfile(args.fastq_list):
            list_file = args.fastq_list
        elif args.fastq_dir:
            list_file = prepare_input_list(args.fastq_dir)
        else:
            logger.error("Provide --fastq-list or --fastq-dir for QC")
            sys.exit(1)

        qc_result = qc.run_pipeline(list_file)
        if qc_result['status'] != 'success':
            logger.error("QC step failed")
            sys.exit(1)

        trimmed_files = qc_result['trimmed_files']
        logger.info(f"QC complete. {len(trimmed_files)} trimmed samples.")
    else:
        logger.info("Skipping QC (not in --stages).")

    # ===================================================================
    #  STAGE: collapse
    # ===================================================================
    collapsed_fastas = {}

    if 'collapse' in stages:
        logger.info("=" * 50)
        logger.info("  STAGE: Sequence Collapsing")
        logger.info("=" * 50)

        collapse_output = output_dir / "collapse_results"
        collapser = FastaCollapser(output_dir=str(collapse_output), min_count=args.min_count)

        fastq_paths = []
        if trimmed_files:
            fastq_paths = [str(fp) for fp in trimmed_files.values()]
        elif args.fastq_list and os.path.isfile(args.fastq_list):
            with open(args.fastq_list, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        parts = line.split()
                        fastq_paths.append(parts[-1] if len(parts) > 1 else line)
        else:
            logger.error("No FASTQ files found. Provide --fastq-list.")
            sys.exit(1)

        if not fastq_paths:
            logger.error("No FASTQ files found")
            sys.exit(1)

        sample_map = {}
        for fp in fastq_paths:
            fp = Path(fp)
            if fp.exists():
                sample_name = collapser._get_base_name(fp)
                sample_name = sample_name.replace('_trimmed', '').replace('_val_1', '')
                sample_map[sample_name] = str(fp)

        logger.info(f"Collapsing {len(sample_map)} samples...")
        result = collapser.process_multiple_samples(sample_map)
        for sn, res in result['samples'].items():
            if res['status'] == 'success':
                collapsed_fastas[sn] = str(res['fasta_file'])
        logger.info(f"Collapse complete. {len(collapsed_fastas)} samples ready.")
    else:
        logger.info("Skipping collapse (not in --stages).")
        # Load pre-collapsed if available
        if args.collapsed_dir:
            coll_dir = Path(args.collapsed_dir)
            for f in coll_dir.glob("*_collapsed.fasta"):
                sn = f.stem.replace('_collapsed', '')
                collapsed_fastas[sn] = str(f)
            if not collapsed_fastas:
                for f in coll_dir.glob("*.fa"):
                    collapsed_fastas[f.stem] = str(f)
            logger.info(f"Loaded {len(collapsed_fastas)} pre-collapsed FASTAs from: {coll_dir}")

    # ===================================================================
    #  STAGE: annotation
    # ===================================================================
    count_files = {}

    if 'annotation' in stages:
        logger.info("=" * 50)
        logger.info("  STAGE: tsRNA Annotation (Bowtie)")
        logger.info("=" * 50)

        if not collapsed_fastas:
            logger.error("No collapsed FASTAs available. Run 'collapse' first or provide --collapsed-dir.")
            sys.exit(1)

        annot_output = output_dir / "tsRNA_results"
        analyzer = SimpleTsRNAAnalyzer(
            index_dir=args.index_dir,
            output_dir=str(annot_output),
            min_length=args.min_len,
            max_length=args.max_len,
            mismatch=args.mismatch,
            threads=args.threads,
            bowtie_path=args.bowtie_path,
            keep_temp=args.keep_temp
        )

        for sample_name, fasta_path in collapsed_fastas.items():
            logger.info(f"Annotating: {sample_name}")
            result = analyzer.analyze_sample(sample_name, Path(fasta_path))
            if result['status'] == 'success':
                count_files[sample_name] = result['count_file']
            else:
                logger.warning(f"Annotation failed for {sample_name}: {result.get('message', 'unknown')}")
        logger.info(f"Annotation complete. {len(count_files)} samples annotated.")
    else:
        logger.info("Skipping annotation (not in --stages).")
        # Try to load existing count files
        existing = list((output_dir / "tsRNA_results").glob("*_counts.tsv"))
        for f in existing:
            count_files[f.stem.replace('_counts', '')] = str(f)
        if count_files:
            logger.info(f"Found {len(count_files)} existing count files.")

    # ===================================================================
    #  STAGE: merge (always runs if we have count_files or need matrix downstream)
    # ===================================================================
    counts_matrix_path = str(output_dir / "counts_matrix.tsv")

    needs_matrix = bool({'cluster', 'deseq2', 'extend'} & stages)
    if count_files:
        logger.info("=" * 50)
        logger.info("  STAGE: Merging Count Matrices")
        logger.info("=" * 50)
        merge_counts(count_files, counts_matrix_path)
    elif os.path.exists(counts_matrix_path):
        logger.info("Using existing counts_matrix.tsv")
    elif needs_matrix:
        logger.warning("No count files to merge; creating empty fallback.")
        pd.DataFrame().to_csv(counts_matrix_path, sep='\t')

    # ===================================================================
    #  STAGE: cluster
    # ===================================================================
    cluster_dir = str(output_dir / "cluster_results")
    cluster_summary = os.path.join(cluster_dir, "cluster_summary_simplified.csv")

    if 'cluster' in stages:
        logger.info("=" * 50)
        logger.info("  STAGE: Sequence Clustering")
        logger.info("=" * 50)

        deg_matrix, cluster_summary = run_clustering(
            input_matrix_path=counts_matrix_path,
            output_dir=cluster_dir,
            method=args.cluster_method
        )
        cluster_matrix_path = deg_matrix
        logger.info(f"Clustering complete. DEG matrix: {deg_matrix}")
    else:
        logger.info("Skipping clustering (not in --stages).")
        cluster_matrix_path = os.path.join(cluster_dir, "clustered_counts_for_DEG.csv")
        if not os.path.exists(cluster_matrix_path):
            logger.warning("Clustered counts not found; falling back to raw counts_matrix.tsv")
            cluster_matrix_path = counts_matrix_path
        else:
            logger.info("Using existing clustered counts for DESeq2")

    # ===================================================================
    #  STAGE: deseq2
    # ===================================================================
    deseq2_dir = str(output_dir / "deseq2_results")
    deseq2_output_path = os.path.join(deseq2_dir, "deseq2_results.tsv")

    if 'deseq2' in stages:
        logger.info("=" * 50)
        logger.info("  STAGE: Differential Expression (DESeq2)")
        logger.info("=" * 50)

        if not args.metadata or not args.contrast:
            logger.error("DESeq2 requires --metadata and --contrast")
            sys.exit(1)

        os.makedirs(deseq2_dir, exist_ok=True)

        # Split clustered matrix into per-sample files
        split_dir = os.path.join(deseq2_dir, "sample_counts")
        sample_files = split_clustered_matrix(cluster_matrix_path, split_dir)
        logger.info(f"Samples found in count matrix: {list(sample_files.keys())}")

        meta_df = pd.read_csv(args.metadata, sep='\t', header=None,
                              names=['sample_name', 'condition'])
        logger.info(f"Samples from --metadata: {meta_df['sample_name'].tolist()}")

        metadata_for_deseq = os.path.join(deseq2_dir, "metadata_deseq2.tsv")
        matched = 0
        with open(metadata_for_deseq, 'w') as f:
            for _, row in meta_df.iterrows():
                sn = row['sample_name']
                cond = row['condition']
                if sn in sample_files:
                    f.write(f"{sample_files[sn]}\t{cond}\n")
                    matched += 1
                else:
                    logger.warning(f"Sample '{sn}' from metadata not in count matrix")

        if matched == 0:
            logger.error(
                "None of the sample names in --metadata match the count matrix columns.\n"
                "  -> Sample names in the count matrix are derived from FASTQ filenames\n"
                "     (stripped of extensions, _trimmed, _val_1, _collapsed suffixes).\n"
                "  -> Update your metadata file to use these exact sample names listed above."
            )
            sys.exit(1)

        try:
            meta_df_deseq = pd.read_csv(metadata_for_deseq, sep='\t', header=None,
                                         names=['file_path', 'condition'])

            if meta_df_deseq.empty:
                logger.error("No matching samples between metadata and count matrix. "
                             "Check that sample names in --metadata match count matrix column names.")
                sys.exit(1)

            count_dfs = []
            clinical_data = []

            for _, row in meta_df_deseq.iterrows():
                f_path = Path(row['file_path'])
                cond = row['condition']
                if not f_path.exists():
                    logger.warning(f"Skipping missing: {f_path}")
                    continue
                sample_id = f_path.stem
                if sample_id.endswith("_counts"):
                    sample_id = sample_id[:-7]
                df = pd.read_csv(f_path, sep='\t', index_col=0)
                if 'count' in df.columns:
                    df = df[['count']]
                elif len(df.columns) == 1:
                    pass
                else:
                    continue
                df.columns = [sample_id]
                count_dfs.append(df)
                clinical_data.append({'sample': sample_id, args.design: cond})

            if not count_dfs:
                logger.error("DESeq2 failed: no valid per-sample count files could be read. "
                             "Check that the count matrix contains samples listed in --metadata.")
                sys.exit(1)

            counts = pd.concat(count_dfs, axis=1, join='outer').fillna(0).astype(int)
            clinical_df = pd.DataFrame(clinical_data).set_index('sample')
            common = [s for s in clinical_df.index if s in counts.columns]
            counts = counts[common]
            clinical_df = clinical_df.loc[common]

            logger.info(f"DESeq2 input: {counts.shape[0]} features x {counts.shape[1]} samples")

            keep = counts.sum(axis=1) >= args.min_count_deseq2
            counts = counts[keep]
            counts = counts.T

            from pydeseq2.dds import DeseqDataSet
            from pydeseq2.ds import DeseqStats
            import inspect

            dds_params = inspect.signature(DeseqDataSet.__init__).parameters
            cpu_kwargs = {'n_cpus': 4} if 'n_cpus' in dds_params else {'n_jobs': 4}
            dds = DeseqDataSet(
                counts=counts,
                metadata=clinical_df,
                design_factors=[args.design],
                refit_cooks=True,
                **cpu_kwargs
            )
            dds.deseq2()

            stat_res = DeseqStats(dds, contrast=[args.design, args.contrast[0], args.contrast[1]])
            stat_res.summary()
            results_df = stat_res.results_df
            results_df.to_csv(deseq2_output_path, sep='\t')

            sig_up = results_df[(results_df['pvalue'] < args.pvalue_thresh) & (results_df['log2FoldChange'] > 0)].shape[0]
            sig_down = results_df[(results_df['pvalue'] < args.pvalue_thresh) & (results_df['log2FoldChange'] < 0)].shape[0]
            logger.info(f"DESeq2 complete. Up: {sig_up}, Down: {sig_down}")

            volcano_path = os.path.join(deseq2_dir, "volcano_plot.png")
            plot_volcano(results_df, volcano_path, pvalue_thresh=args.pvalue_thresh)

            heatmap_path = os.path.join(deseq2_dir, "heatmap.png")
            plot_heatmap(dds, results_df, clinical_df, heatmap_path, top_n=args.top_n_heatmap, pvalue_thresh=args.pvalue_thresh)
            logger.info(f"Plots saved to {deseq2_dir}")

            if args.output_normalized:
                normed_path = args.output_normalized
                if not os.path.isabs(normed_path):
                    normed_path = os.path.join(deseq2_dir, normed_path)
                normed = dds.layers['normed_counts']
                if isinstance(normed, pd.DataFrame):
                    normed.T.to_csv(normed_path, sep='\t')
                else:
                    pd.DataFrame(normed.T, index=dds.var_names, columns=dds.obs_names) \
                        .to_csv(normed_path, sep='\t')
                logger.info(f"Normalized counts saved to {normed_path}")

        except Exception as e:
            logger.error(f"DESeq2 failed: {e}")
            import traceback
            traceback.print_exc()
            sys.exit(1)
    else:
        logger.info("Skipping DESeq2 (not in --stages).")

    # ===================================================================
    #  STAGE: extend (extension + primer search)
    # ===================================================================
    extend_output = None
    sig_list_path = None
    extend_dir = str(output_dir / "extend_primer_results")

    if 'extend' in stages:
        logger.info("=" * 50)
        logger.info("  STAGE: Extending tsRNA Sequences + Primer Search")
        logger.info("=" * 50)

        if not os.path.exists(deseq2_output_path):
            logger.error("DESeq2 results not found. Run 'deseq2' first or provide existing results.")
            sys.exit(1)

        os.makedirs(extend_dir, exist_ok=True)

        sig_list_path = os.path.join(extend_dir, "significant_tsRNAs.txt")
        generate_significant_list(deseq2_output_path, sig_list_path, pvalue_thresh=args.pvalue_thresh)

        tRNA_fasta = args.tRNA_fasta
        if not tRNA_fasta:
            idx_dir = Path(args.index_dir)
            for candidate in ['hg38-tRNA.fa', 'mature.fa', 'tRNA.fa']:
                if (idx_dir / candidate).exists():
                    tRNA_fasta = str(idx_dir / candidate)
                    break
            if not tRNA_fasta:
                fa_files = [f for f in idx_dir.glob("*.fa") if 'collapsed' not in f.name]
                if fa_files:
                    tRNA_fasta = str(fa_files[0])

        if not tRNA_fasta or not os.path.exists(tRNA_fasta):
            logger.error("tRNA FASTA not found. Provide --tRNA-fasta.")
            sys.exit(1)

        # Step 1: Extend tsRNA sequences
        extend_output = os.path.join(extend_dir, "output_result.csv")
        extend_tsrna_sequences(
            fa_file=tRNA_fasta,
            list_file=sig_list_path,
            count_csv=cluster_summary,
            out_file=extend_output,
            extend_by=args.extend_by
        )

        # Step 2: Primer / similar sequence search
        primer_output = os.path.join(extend_dir, "clustered_tsRNA_similar_sequences.csv")
        run_primer_search(
            extend_csv=extend_output,
            count_matrix_path=counts_matrix_path,
            output_file=primer_output
        )

    # ===================================================================
    #  STAGE: predict
    # ===================================================================
    if 'predict' in stages:
        logger.info("=" * 50)
        logger.info("  STAGE: Target Prediction & Enrichment")
        logger.info("=" * 50)

        # Derive significant list if not already generated
        if not sig_list_path or not os.path.exists(sig_list_path):
            sig_list_path = os.path.join(extend_dir, "significant_tsRNAs.txt")
            if not os.path.exists(sig_list_path) and os.path.exists(deseq2_output_path):
                os.makedirs(os.path.dirname(sig_list_path), exist_ok=True)
                generate_significant_list(deseq2_output_path, sig_list_path, pvalue_thresh=args.pvalue_thresh)

        if not sig_list_path or not os.path.exists(sig_list_path):
            logger.error("No significant tsRNA list available. Run 'deseq2' or 'extend' first.")
            sys.exit(1)

        predict_dir = str(output_dir / "prediction_results")
        os.makedirs(predict_dir, exist_ok=True)

        script_dir = os.path.dirname(os.path.abspath(__file__))
        cmd = [
            sys.executable, os.path.join(script_dir, "modules", "prediction_enrichment.py"),
            "--list", sig_list_path,
            "--index", args.predict_index,
            "--energy", str(args.energy),
            "--threads", str(args.threads),
            "--output_dir", predict_dir,
        ]
        if args.top_percent is not None:
            cmd.extend(["--top-percent", str(args.top_percent)])
        else:
            cmd.extend(["--threshold", str(args.threshold)])
        logger.info(f"Running: {' '.join(cmd)}")
        result = subprocess.run(cmd, check=False)
        if result.returncode != 0:
            logger.error(f"Prediction step failed with code {result.returncode}")
        else:
            logger.info("Prediction complete.")

    # ===================================================================
    logger.info("=" * 50)
    logger.info("  HERACLES Pipeline Complete!")
    logger.info(f"  Results in: {output_dir.resolve()}")
    logger.info("=" * 50)


if __name__ == "__main__":
    main()
