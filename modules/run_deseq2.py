#!/usr/bin/env python3
"""
Differential Expression Analysis using pydeseq2
Includes Volcano Plot and Heatmap generation.
"""

import argparse
import logging
import sys
import pandas as pd
import numpy as np
from pathlib import Path

# Plotting libraries
try:
    import matplotlib.pyplot as plt
    import seaborn as sns
    # Set plotting style
    sns.set(style="whitegrid")
except ImportError:
    print("Error: matplotlib or seaborn not installed.")
    print("Please run: pip install matplotlib seaborn")
    sys.exit(1)

# DESeq2 libraries
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def plot_volcano(results_df, output_path, pvalue_thresh=0.05, lfc_thresh=1.0):
    """
    Draw Volcano Plot
    """
    logger.info("Generating Volcano Plot...")
    
    plt.figure(figsize=(10, 8))
    
    pvalue_clean = results_df['pvalue'].replace(0, 1e-300)
    results_df['-log10_pvalue'] = -np.log10(pvalue_clean)
    
    # Define group colors
    up = results_df[(results_df['pvalue'] < pvalue_thresh) & (results_df['log2FoldChange'] >= lfc_thresh)]
    down = results_df[(results_df['pvalue'] < pvalue_thresh) & (results_df['log2FoldChange'] <= -lfc_thresh)]
    ns = results_df[(results_df['pvalue'] >= pvalue_thresh) | (abs(results_df['log2FoldChange']) < lfc_thresh)]
    
    # Plot scatter points
    plt.scatter(ns['log2FoldChange'], ns['-log10_pvalue'], s=20, color='grey', alpha=0.5, label='Not Significant')
    plt.scatter(up['log2FoldChange'], up['-log10_pvalue'], s=20, color='red', alpha=0.7, label='Up-regulated')
    plt.scatter(down['log2FoldChange'], down['-log10_pvalue'], s=20, color='blue', alpha=0.7, label='Down-regulated')
    
    # Add threshold lines
    plt.axhline(-np.log10(pvalue_thresh), color='black', linestyle='--', linewidth=1, alpha=0.6)
    plt.axvline(lfc_thresh, color='black', linestyle='--', linewidth=1, alpha=0.6)
    plt.axvline(-lfc_thresh, color='black', linestyle='--', linewidth=1, alpha=0.6)
    
    # Set title and labels
    plt.title(f'Volcano Plot (P-value < {pvalue_thresh}, |log2FC| > {lfc_thresh})', fontsize=14)
    plt.xlabel('log2 Fold Change', fontsize=12)
    plt.ylabel('-log10 p-value', fontsize=12)
    plt.legend()
    
    # Save figure
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Volcano plot saved to {output_path}")


def plot_heatmap(dds, results_df, clinical_df, output_path, top_n=50):
    """
    Draw Heatmap (using normalized counts)
    Default shows top_n most significant genes
    """
    logger.info(f"Generating Heatmap (Top {top_n} significant genes)...")
    
    # 1. Get normalized counts matrix
    # dds.layers['normed_counts'] has shape (n_obs, n_vars) -> (n_samples, n_genes)
    norm_counts = dds.layers['normed_counts'].copy()
    
    # Convert to DataFrame
    # CRITICAL FIX: index=obs_names (samples), columns=var_names (genes)
    norm_df = pd.DataFrame(norm_counts, index=dds.obs_names, columns=dds.var_names)
    
    # Transpose to (Genes x Samples) for standard heatmap orientation (Genes on Y-axis)
    norm_df = norm_df.T
    
    # 2. Filter significantly differentially expressed genes
    sig_df = results_df[results_df['pvalue'] < 0.05].copy()
    
    if sig_df.empty:
        logger.warning("No significant genes found for heatmap. Skipping heatmap generation.")
        return

    # 3. Sort by absolute log2FoldChange, take top_n
    sig_df = sig_df.reindex(sig_df['log2FoldChange'].abs().sort_values(ascending=False).index)
    top_genes = sig_df.index[:top_n]
    
    # Check if top_genes exist in norm_df (might have been filtered out during low count filtering if logic differs, though usually consistent)
    valid_top_genes = [g for g in top_genes if g in norm_df.index]
    if not valid_top_genes:
        logger.warning("Top genes not found in normalized matrix. Skipping heatmap.")
        return
        
    # 4. Extract these genes from normalized matrix
    heatmap_data = norm_df.loc[valid_top_genes]
    
    # 5. Log2 transform for visualization (add 1 to avoid log0)
    heatmap_data = np.log2(heatmap_data + 1)
    
    # 6. Plot heatmap
    # Use clustermap for clustering
    # row_cluster=True: cluster genes
    # col_cluster=True: cluster samples
    # cmap='vlag': blue-white-red colormap, suitable for up/down regulation
    g = sns.clustermap(
        heatmap_data, 
        cmap='vlag', 
        center=0, 
        linewidths=.5, 
        figsize=(12, 8 + 0.2 * len(valid_top_genes)), # Adjust height dynamically
        row_cluster=True, 
        col_cluster=True,
        yticklabels=True, 
        xticklabels=True
    )
    
    # Adjust x-axis label rotation
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right')
    
    # Set title
    g.fig.suptitle(f'Expression Heatmap of Top {len(valid_top_genes)} Significant Genes', y=1.02)
    
    # Save figure
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Heatmap saved to {output_path}")


def main():
    parser = argparse.ArgumentParser(description='Run DESeq2 analysis on tsRNA count files')
    
    parser.add_argument('-m', '--metadata', required=True, help='Metadata file: <path_to_count_file>\\t<condition>')
    parser.add_argument('-o', '--output', default='deseq2_results', help='Output directory')
    parser.add_argument('--design', default='condition', help='Design formula factor name (default: condition)')
    parser.add_argument('--contrast', nargs=2, required=True, metavar=('TREATMENT', 'CONTROL'), 
                        help='Comparison groups: e.g., --contrast Treat Control')
    parser.add_argument('--min-count', type=int, default=10, help='Minimum total count to keep a gene')
    parser.add_argument('--top-n', type=int, default=50, help='Number of top genes to show in heatmap')
    
    args = parser.parse_args()
    
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # --- 1. Load Metadata and Data ---
    logger.info("Loading data...")
    
    if not Path(args.metadata).exists():
        logger.error(f"Metadata file not found: {args.metadata}")
        sys.exit(1)

    meta_df = pd.read_csv(args.metadata, sep='\t', header=None, names=['file_path', 'condition'])
    
    count_dfs = []
    clinical_data = [] 
    
    for idx, row in meta_df.iterrows():
        f_path = Path(row['file_path'])
        cond = row['condition']
        
        if not f_path.exists():
            logger.warning(f"Skipping missing file: {f_path}")
            continue
            
        sample_id = f_path.stem
        if sample_id.endswith("_counts"):
            sample_id = sample_id[:-7] 
            
        try:
            df = pd.read_csv(f_path, sep='\t', index_col=0)
            
            if 'count' in df.columns:
                df = df[['count']]
            elif len(df.columns) == 1:
                pass 
            else:
                logger.warning(f"Unexpected format in {f_path}, skipping")
                continue
                
            df.columns = [sample_id]
            count_dfs.append(df)
            
            clinical_data.append({
                'sample': sample_id,
                args.design: cond
            })
            
        except Exception as e:
            logger.error(f"Error reading {f_path}: {e}")
            continue

    if not count_dfs:
        logger.error("No data loaded.")
        sys.exit(1)

    # --- Merge Data ---
    counts_matrix = pd.concat(count_dfs, axis=1, join='outer').fillna(0).astype(int)
    
    # Build Clinical DataFrame
    clinical_df = pd.DataFrame(clinical_data)
    clinical_df = clinical_df.set_index('sample')
    
    samples_in_clinical = clinical_df.index.tolist()
    common_samples = [s for s in samples_in_clinical if s in counts_matrix.columns]
    counts_matrix = counts_matrix[common_samples]
    clinical_df = clinical_df.loc[common_samples]
    
    logger.info(f"Data loaded: {counts_matrix.shape[0]} tsRNAs x {counts_matrix.shape[1]} samples")
    logger.info(f"Conditions: {clinical_df[args.design].unique()}")
    
    # --- 2. Filtering ---
    logger.info(f"Filtering genes with total counts < {args.min_count}...")
    keep_genes = counts_matrix.sum(axis=1) >= args.min_count
    counts_matrix = counts_matrix[keep_genes]
    logger.info(f"Remaining genes after filtering: {len(counts_matrix)}")
    
    count_matrix_output = output_dir / "counts_matrix.tsv"
    counts_matrix.to_csv(count_matrix_output, sep='\t')
    logger.info(f"Filtered count matrix saved to: {count_matrix_output}")

    # --- 3. Transpose Matrix ---
    logger.info("Transposing matrix for pydeseq2...")
    counts_matrix = counts_matrix.T
    logger.info(f"Counts matrix shape after transpose: {counts_matrix.shape}")

    # --- 4. Run DESeq2 ---
    logger.info("Initializing DESeq2...")
    
    try:
        dds = DeseqDataSet(
            counts=counts_matrix,
            metadata=clinical_df,
            design_factors=[args.design],
            refit_cooks=True,
            n_cpus=4
        )
        
        logger.info("Running DESeq2...")
        dds.deseq2()
        
        logger.info("Extracting results...")
        stat_res = DeseqStats(
            dds,
            contrast=[args.design, args.contrast[0], args.contrast[1]]
        )
        stat_res.summary()
        
        results_df = stat_res.results_df
        results_df = results_df.dropna(subset=['pvalue'])
        
        # Save results
        output_file = output_dir / "deseq2_results.tsv"
        results_df.to_csv(output_file, sep='\t')
        logger.info(f"Results saved to {output_file}")
        
        significant_up = results_df[(results_df['pvalue'] < 0.05) & (results_df['log2FoldChange'] > 0)].shape[0]
        significant_down = results_df[(results_df['pvalue'] < 0.05) & (results_df['log2FoldChange'] < 0)].shape[0]
        
        logger.info(f"Analysis Complete.")
        logger.info(f"  Significant Up-regulated (pvalue < 0.05): {significant_up}")
        logger.info(f"  Significant Down-regulated (pvalue < 0.05): {significant_down}")
        
        # --- 5. Plotting ---
        logger.info("Generating plots...")
        
        volcano_path = output_dir / "volcano_plot.png"
        plot_volcano(results_df, volcano_path)
        
        heatmap_path = output_dir / "heatmap.png"
        plot_heatmap(dds, results_df, clinical_df, heatmap_path, top_n=args.top_n)
        
    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
