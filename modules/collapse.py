import gzip
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Union, Sequence
from collections import Counter
import pandas as pd

# Setting up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class FastaCollapser:
    """
    Module for collapsing identical sequences from trimmed fastq files.
    Reads trimmed fastq files, counts identical sequences, and generates fasta files with counts.
    """
    
    def __init__(self, 
                 output_dir: str = "collapse_results",
                 min_count: int = 1):
        """
        Initialize the Collapser
        
        Args:
            output_dir: Output directory for collapsed fasta files
            min_count: Minimum sequence count threshold, sequences below this will be filtered out
        """
        self.output_dir = Path(output_dir)
        self.min_count = min_count
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"FastaCollapser initialized, output directory: {self.output_dir}")
    
    def _get_base_name(self, file_path: Union[str, Path]) -> str:
        """
        Get base name from file path without extensions (handles .fq.gz, .fastq.gz properly)
        """
        name = Path(file_path).name
        extensions = ['.fastq.gz', '.fq.gz', '.fastq', '.fq', '.gz']
        for ext in extensions:
            if name.endswith(ext):
                name = name[:-len(ext)]
                break
        return name

    def read_fastq(self, fastq_path: Union[str, Path]) -> List[str]:
        """
        Read fastq file and return list of sequences
        
        Args:
            fastq_path: Path to fastq file (can be gzipped)
            
        Returns:
            List of sequences
        """
        sequences = []
        fastq_path = Path(fastq_path)
        
        logger.info(f"Reading fastq file: {fastq_path}")
        
        if not fastq_path.exists():
            logger.error(f"Fastq file does not exist: {fastq_path}")
            raise FileNotFoundError(f"Fastq file does not exist: {fastq_path}")
        
        try:
            # Safely open file based on extension
            open_func = gzip.open if str(fastq_path).endswith('.gz') else open
            mode = 'rt' if str(fastq_path).endswith('.gz') else 'r'
            
            with open_func(fastq_path, mode) as f:
                # Safer parsing: read exactly 4 lines per record
                while True:
                    header = f.readline()
                    if not header:
                        break  # End of file
                    
                    seq = f.readline().strip()
                    f.readline()  # Ignore '+' line
                    f.readline()  # Ignore quality line
                    
                    if seq:  # Ensure sequence is not empty
                        sequences.append(seq)
            
            logger.info(f"Successfully read {len(sequences)} sequences")
            return sequences
            
        except Exception as e:
            logger.error(f"Failed to read fastq file {fastq_path}: {str(e)}")
            raise
    
    def collapse_sequences(self, sequences: List[str]) -> Dict[str, int]:
        """
        Count occurrences of identical sequences
        
        Args:
            sequences: List of sequences
            
        Returns:
            Dictionary: {sequence: count}
        """
        logger.info("Starting sequence counting...")
        
        # Use Counter to count sequence occurrences
        sequence_counter = Counter(sequences)
        
        # Filter out sequences below threshold
        collapsed_seqs = {
            seq: count for seq, count in sequence_counter.items() 
            if count >= self.min_count
        }
        
        total_unique_seqs = len(collapsed_seqs)
        total_original_seqs = len(sequences)
        
        logger.info(f"Sequence counting completed:")
        logger.info(f"  Original sequences: {total_original_seqs}")
        logger.info(f"  Unique sequences: {total_unique_seqs}")
        
        # [FIX] Prevent ZeroDivisionError if file was empty
        if total_original_seqs > 0:
            logger.info(f"  Compression rate: {total_unique_seqs/total_original_seqs:.2%}")
        else:
            logger.info("  Compression rate: N/A (no sequences found)")
        
        return collapsed_seqs
    
    def write_collapsed_fasta(self, 
                            collapsed_seqs: Dict[str, int], 
                            sample_name: str,
                            output_prefix: str = "") -> Path:
        """
        Write collapsed sequences to fasta file
        
        Args:
            collapsed_seqs: Dictionary {sequence: count}
            sample_name: Sample name
            output_prefix: Prefix for output filename
            
        Returns:
            Path to output fasta file
        """
        if output_prefix:
            output_file = self.output_dir / f"{output_prefix}_{sample_name}_collapsed.fasta"
        else:
            output_file = self.output_dir / f"{sample_name}_collapsed.fasta"
        
        logger.info(f"Writing collapsed fasta file: {output_file}")
        
        # Sort sequences by count (descending)
        sorted_seqs = sorted(collapsed_seqs.items(), key=lambda x: x[1], reverse=True)
        
        with open(output_file, 'w') as f:
            for i, (seq, count) in enumerate(sorted_seqs, 1):
                # Fasta header format: >sample_sequenceID_count
                header = f">{sample_name}_seq{i}_{count}"
                f.write(f"{header}\n")
                f.write(f"{seq}\n")
        
        logger.info(f"Successfully wrote {len(sorted_seqs)} unique sequences")
        return output_file
    
    def create_statistics_report(self, 
                               collapsed_seqs: Dict[str, int], 
                               sample_name: str) -> pd.DataFrame:
        """
        Create statistics report for collapsed sequences
        
        Args:
            collapsed_seqs: Dictionary {sequence: count}
            sample_name: Sample name
            
        Returns:
            DataFrame containing statistics
        """
        logger.info(f"Creating statistics report for sample {sample_name}")
        
        # Calculate statistics
        total_unique_seqs = len(collapsed_seqs)
        total_reads = sum(collapsed_seqs.values())
        
        if collapsed_seqs:
            max_count = max(collapsed_seqs.values())
            min_count = min(collapsed_seqs.values())
            avg_count = total_reads / total_unique_seqs
            
            # Count sequences in different count ranges
            count_ranges = {
                '1': 0, '2-5': 0, '6-10': 0, '11-50': 0, 
                '51-100': 0, '101-1000': 0, '>1000': 0
            }
            
            for count in collapsed_seqs.values():
                if count == 1:
                    count_ranges['1'] += 1
                elif 2 <= count <= 5:
                    count_ranges['2-5'] += 1
                elif 6 <= count <= 10:
                    count_ranges['6-10'] += 1
                elif 11 <= count <= 50:
                    count_ranges['11-50'] += 1
                elif 51 <= count <= 100:
                    count_ranges['51-100'] += 1
                elif 101 <= count <= 1000:
                    count_ranges['101-1000'] += 1
                else:
                    count_ranges['>1000'] += 1
        else:
            max_count = min_count = avg_count = 0
            count_ranges = {k: 0 for k in [
                '1', '2-5', '6-10', '11-50', '51-100', '101-1000', '>1000'
            ]}
        
        # Create DataFrame
        stats_data = {
            'Sample': sample_name,
            'Total_Unique_Sequences': total_unique_seqs,
            'Total_Reads': total_reads,
            'Max_Count': max_count,
            'Min_Count': min_count,
            'Average_Count': avg_count,
            'Sequences_with_Count_1': count_ranges['1'],
            'Sequences_with_Count_2-5': count_ranges['2-5'],
            'Sequences_with_Count_6-10': count_ranges['6-10'],
            'Sequences_with_Count_11-50': count_ranges['11-50'],
            'Sequences_with_Count_51-100': count_ranges['51-100'],
            'Sequences_with_Count_101-1000': count_ranges['101-1000'],
            'Sequences_with_Count_>1000': count_ranges['>1000']
        }
        
        return pd.DataFrame([stats_data])
    
    def process_sample(self, 
                      fastq_path: Union[str, Path], 
                      sample_name: str,
                      output_prefix: str = "") -> Dict:
        """
        Process single sample
        
        Args:
            fastq_path: Path to fastq file
            sample_name: Sample name
            output_prefix: Prefix for output files
            
        Returns:
            Dictionary with processing results
        """
        logger.info(f"Processing sample: {sample_name}")
        
        try:
            # 1. Read fastq file
            sequences = self.read_fastq(fastq_path)
            
            # 2. Collapse sequences
            collapsed_seqs = self.collapse_sequences(sequences)
            
            # 3. Write fasta file
            fasta_file = self.write_collapsed_fasta(
                collapsed_seqs, sample_name, output_prefix
            )
            
            # 4. Generate statistics report
            stats_df = self.create_statistics_report(collapsed_seqs, sample_name)
            
            result = {
                'status': 'success',
                'sample_name': sample_name,
                'fasta_file': fasta_file,
                'statistics': stats_df,
                'total_sequences': len(sequences),
                'unique_sequences': len(collapsed_seqs),
                'collapsed_sequences': collapsed_seqs
            }
            
            logger.info(f"Sample {sample_name} processing completed")
            return result
            
        except Exception as e:
            logger.error(f"Failed to process sample {sample_name}: {str(e)}")
            return {
                'status': 'failed',
                'sample_name': sample_name,
                'error': str(e)
            }
    
    def process_multiple_samples(self, 
                               fastq_files: Dict[str, Union[Path, Sequence[Path]]],
                               output_prefix: str = "") -> Dict:
        """
        Process multiple samples
        
        Args:
            fastq_files: Dictionary of samples, format: {sample_name: fastq_path}
                         or {sample_name: (R1_path, R2_path)}
            output_prefix: Prefix for output files
            
        Returns:
            Dictionary with processing results
        """
        logger.info(f"Starting to process {len(fastq_files)} samples")
        
        results = {}
        all_stats = []
        
        for sample_name, fastq_path in fastq_files.items():
            logger.info(f"Processing sample: {sample_name}")
            
            # [FIX] Handle paired-end data or lists: only process R1
            if isinstance(fastq_path, (tuple, list)):
                logger.info(f"Sample {sample_name} appears to be paired-end or a list, processing index 0 only")
                fastq_path = fastq_path[0]  # Take only R1
            
            # Process single sample
            result = self.process_sample(fastq_path, sample_name, output_prefix)
            results[sample_name] = result
            
            # Collect statistics
            if result['status'] == 'success':
                all_stats.append(result['statistics'])
        
        # Combine all statistics
        if all_stats:
            combined_stats = pd.concat(all_stats, ignore_index=True)
            stats_file = self.output_dir / "collapse_statistics.csv"
            combined_stats.to_csv(stats_file, index=False)
            logger.info(f"Statistics saved to: {stats_file}")
        else:
            combined_stats = pd.DataFrame()
            logger.warning("No successful sample statistics available")
        
        overall_result = {
            'status': 'success' if any(r['status'] == 'success' for r in results.values()) else 'failed',
            'samples': results,
            'combined_statistics': combined_stats,
            'output_dir': str(self.output_dir)
        }
        
        logger.info("All samples processed")
        return overall_result


def main():
    """
    Example usage
    """
    import sys
    
    collapser = FastaCollapser(
        output_dir="collapse_results",
        min_count=1
    )
    
    # Process a trimmed fastq file
    if len(sys.argv) > 1:
        fastq_file = sys.argv[1]
        
        # [FIX] Use the robust base name extractor instead of Path.stem
        sample_name = collapser._get_base_name(fastq_file)
        sample_name = sample_name.replace('_trimmed', '').replace('_val_1', '')
        
        result = collapser.process_sample(fastq_file, sample_name)
        
        if result['status'] == 'success':
            print(f"Processing completed!")
            print(f"FASTA file: {result['fasta_file']}")
            print(f"Statistics:\n{result['statistics'].to_string()}")
    else:
        print("Usage: python collapse.py <fastq_file>")


if __name__ == "__main__":
    main()
