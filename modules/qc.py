import os
import gzip
import subprocess
import tempfile
import logging
import argparse
from pathlib import Path
from typing import List, Dict, Optional
import pandas as pd
import multiqc

# Setting up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Number of raw reads inspected when inferring the sequencing cycle count of a library.
CYCLE_SAMPLE_READS = 100_000


class QCProcessor:
    """
    Quality control processing module for single-end sequencing, 
    using Trim Galore for quality trimming and MultiQC for report generation.
    """
    
    def __init__(self, 
                 output_dir: str = "qc_results",
                 trim_galore_path: str = "trim_galore",
                 multiqc_path: str = "multiqc",
                 threads: int = 4,
                 quality: int = 20,
                 length: int = 18,
                 adapter: str = "",
                 drop_unsupported: bool = False,
                 cycle_length: Optional[int] = None):
        """
        Initialize QC module

        Args:
            output_dir: Output directory
            trim_galore_path: trim_galore path
            multiqc_path: multiqc path
            threads: Number of threads
            quality: Quality threshold
            length: Minimum read length
            adapter: Adapter sequence for trimming (e.g. 'TGGAATTCTCGGGTGCCAAGG').
                     If empty, Trim Galore auto-detects.
            drop_unsupported: Optional. Drop trimmed reads whose length still equals
                     the library's sequencing cycle count L. Such reads had no bases
                     removed by Trim Galore, i.e. no detectable 3' adapter, so their
                     final nucleotide may be a read-length truncation rather than the
                     authentic RNA 3' terminus (reviewer 1, comment 3). Default OFF:
                     with the flag off the pipeline output is byte-for-byte unchanged.
            cycle_length: Fixed cycle count L for every library. When None (default),
                     L is inferred per library from the raw FASTQ (longest read among
                     the first CYCLE_SAMPLE_READS reads), which accommodates libraries
                     sequenced with different cycle counts.
        """
        self.output_dir = Path(output_dir)
        self.trim_galore_path = trim_galore_path
        self.multiqc_path = multiqc_path
        self.threads = threads
        self.quality = quality
        self.length = length
        self.adapter = adapter
        self.drop_unsupported = drop_unsupported
        self.cycle_length = cycle_length
        self.unsupported_stats: List[Dict] = []

        # Setting output directories
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.trimmed_dir = self.output_dir / "trimmed"
        self.reports_dir = self.output_dir / "reports"
        self.filtered_dir = self.output_dir / "trimmed_unsupported_filtered"
        self.trimmed_dir.mkdir(exist_ok=True)
        self.reports_dir.mkdir(exist_ok=True)

        logger.info(f"QC module initialization complete, output directory: {self.output_dir}")
        if self.drop_unsupported:
            logger.info("Unsupported-endpoint filter ENABLED: reads of length L with no "
                        "detectable 3' adapter will be removed after trimming.")
    
    def read_fastq_list(self, input_file: str) -> List[Path]:
        """
        Read fastq file list from file
        
        Args:
            input_file: File containing fastq paths
            
        Returns:
            List of fastq file paths
        """
        fastq_files = []
        
        try:
            with open(input_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        path = Path(line)
                        if path.exists():
                            fastq_files.append(path)
                            logger.info(f"Located fastq file: {path}")
                        else:
                            logger.warning(f"File does not exist, skipping: {line}")
            
            if not fastq_files:
                raise ValueError("No valid fastq files found")
                
            logger.info(f"Total of {len(fastq_files)} single-end fastq files found.")
            return fastq_files
            
        except FileNotFoundError:
            logger.error(f"Input file does not exist: {input_file}")
            raise
        except Exception as e:
            logger.error(f"Error reading fastq list: {str(e)}")
            raise

    def _organize_single_end(self, fastq_files: List[Path]) -> Dict[str, Path]:
        """
        Organize single-end fastq files
        Automatically remove residual suffixes like _R1, _1 from filenames to obtain clean sample names
        
        Returns:
            Dictionary {sample_name: fastq_file_path}
        """
        se_files = {}
        
        for file in fastq_files:
            base_name = self._get_base_name(file)
            
            # Remove common single-end/paired-end identifiers to prevent verbose names
            for suffix in ['_R1', '_R2', '_1', '_2', '_SE', '_se']:
                if base_name.endswith(suffix):
                    base_name = base_name[:-len(suffix)]
                    break
            
            # Add suffix if duplicate names occur
            orig_name = base_name
            counter = 1
            while base_name in se_files:
                base_name = f"{orig_name}_{counter}"
                counter += 1
                
            se_files[base_name] = file
            logger.info(f"Single-end sample: {base_name} -> {file.name}")
        
        return se_files
    
    def _get_base_name(self, file_path: Path) -> str:
        """
        Get base name from file path without extensions
        
        Args:
            file_path: Path to fastq file
            
        Returns:
            Base name without extensions
        """
        name = file_path.name
        
        # Remove common fastq extensions
        extensions = ['.fastq.gz', '.fq.gz', '.fastq', '.fq', '.gz']
        for ext in extensions:
            if name.endswith(ext):
                name = name[:-len(ext)]
                break
        
        return name

    @staticmethod
    def _open_text(file_path: Path):
        """Open a FASTQ in text mode, transparently handling gzip compression."""
        if str(file_path).endswith('.gz'):
            return gzip.open(file_path, 'rt')
        return open(file_path, 'r')

    def detect_cycle_length(self, fastq_path: Path,
                            max_reads: int = CYCLE_SAMPLE_READS) -> int:
        """Infer the sequencing cycle count L of one library from its raw FASTQ.

        Raw single-end reads are all L bases long, so the longest read among the
        first `max_reads` records is L. Sampling rather than scanning every read
        keeps this cheap on large libraries while remaining robust to a few
        already-truncated records.
        """
        longest = 0
        with self._open_text(fastq_path) as handle:
            for i, line in enumerate(handle):
                if i % 4 == 1:
                    longest = max(longest, len(line.strip()))
                    if i // 4 >= max_reads:
                        break
        if longest == 0:
            raise ValueError(f"Cannot infer cycle length from empty FASTQ: {fastq_path}")
        return longest

    def filter_unsupported_reads(self, trimmed_file: Path, cycle_len: int,
                                 out_file: Path) -> Dict:
        """Remove adapter-unsupported reads from one trimmed FASTQ.

        A read that still has length L after Trim Galore has had no bases removed
        at all: no 3' adapter was detected (and no quality trimming occurred).
        Its last base may therefore be a read-length truncation instead of the
        authentic RNA 3' end. Reads shorter than L have adapter support and are
        kept. The original trimmed file is left untouched; the filtered copy is
        written to `out_file`.

        Note: this is a per-read QC filter. It is not identical to the
        annotation-level exclusion of 50-nt fragments reported in Table S32,
        and the 49-nt reads whose terminal A was removed by a 1-bp adapter
        overlap are out of scope here (they are shorter than L).
        """
        out_file.parent.mkdir(parents=True, exist_ok=True)
        total = kept = removed = 0
        with self._open_text(trimmed_file) as src, gzip.open(out_file, 'wt') as dst:
            while True:
                header = src.readline()
                if not header:
                    break
                seq = src.readline()
                plus = src.readline()
                qual = src.readline()
                if not qual:
                    raise ValueError(f"Truncated FASTQ record in {trimmed_file}")
                header, seq, plus, qual = (x.rstrip('\n\r') for x in
                                           (header, seq, plus, qual))
                if not header.startswith('@') or not plus.startswith('+'):
                    raise ValueError(f"Malformed FASTQ record in {trimmed_file}: {header}")
                if len(seq) != len(qual):
                    raise ValueError(f"Sequence/quality length mismatch in {trimmed_file}")
                total += 1
                if len(seq) == cycle_len:
                    removed += 1
                    continue
                dst.write(f"{header}\n{seq}\n{plus}\n{qual}\n")
                kept += 1
        return {
            'total_trimmed_reads': total,
            'removed_reads': removed,
            'kept_reads': kept,
            'removed_pct': round(100.0 * removed / total, 4) if total else 0.0,
            'cycle_length': cycle_len,
            'filtered_file': str(out_file),
        }

    def run_trim_galore(self, fastq_files: List[Path]) -> Dict[str, Path]:
        """
        Run Trim Galore for single-end quality trimming
        
        Args:
            fastq_files: List of fastq file paths
            
        Returns:
            Dictionary: {sample_name: trimmed_fastq_file_path}
        """
        logger.info("Running Trim Galore (single-end mode)...")
        
        # Organize single-end files
        se_files = self._organize_single_end(fastq_files)
        trimmed_files = {}
        
        for sample_name, file_path in se_files.items():
            logger.info(f"Processing sample: {sample_name}")
            
            # Build trim_galore command (strictly single-end)
            cmd = [
                self.trim_galore_path,
                '--quality', str(self.quality),
                '--length', str(self.length),
                '--cores', str(self.threads),
                '--output_dir', str(self.trimmed_dir),
                '--gzip',
            ]
            if self.adapter:
                cmd.extend(['--adapter', self.adapter])
            cmd.append(str(file_path))  # Pass single file directly
            
            # Get base name and define expected output filename
            base_name = self._get_base_name(file_path)
            trimmed_file = self.trimmed_dir / f"{base_name}_trimmed.fq.gz"
            trimmed_files[sample_name] = trimmed_file
            
            logger.debug(f"Expected output file: {trimmed_file}")
            
            # Run trim_galore
            try:
                logger.info(f"Running command: {' '.join(cmd)}")
                result = subprocess.run(cmd, check=True, capture_output=True, text=True)
                logger.info(f"Trim Galore completed for sample: {sample_name}")
                logger.debug(f"Output:\n{result.stdout}")
                
                if result.stderr:
                    logger.warning(f"Warnings:\n{result.stderr}")
                    
                # Verify output file exists
                if trimmed_file.exists():
                    logger.info(f"Trimmed file created: {trimmed_file}")
                else:
                    logger.warning(f"Expected trimmed file not found: {trimmed_file}")

                # Optional: drop adapter-unsupported full-cycle reads (reviewer 1, comment 3)
                if self.drop_unsupported and trimmed_file.exists():
                    cycle_len = (self.cycle_length if self.cycle_length is not None
                                 else self.detect_cycle_length(file_path))
                    # Keep the original file name in a sibling directory: collapse
                    # derives sample names from it, so changing the stem would
                    # silently rename samples downstream.
                    filtered_file = self.filtered_dir / trimmed_file.name
                    stats = self.filter_unsupported_reads(trimmed_file, cycle_len,
                                                          filtered_file)
                    trimmed_files[sample_name] = filtered_file
                    self.unsupported_stats.append({
                        'Sample': sample_name,
                        'Raw_file': str(file_path),
                        **stats,
                    })
                    logger.info(
                        f"[{sample_name}] cycle length L={cycle_len}; removed "
                        f"{stats['removed_reads']}/{stats['total_trimmed_reads']} "
                        f"adapter-unsupported reads ({stats['removed_pct']}%) -> "
                        f"{filtered_file}"
                    )

            except subprocess.CalledProcessError as e:
                logger.error(f"Trim Galore failed: {str(e)}")
                logger.error(f"Error output:\n{e.stderr}")
                raise
            except Exception as e:
                logger.error(f"Error running Trim Galore: {str(e)}")
                raise
        
        logger.info("Trim Galore completed for all samples")
        return trimmed_files
    
    def run_multiqc(self, input_dir: Optional[str] = None) -> Path:
        """
        Run MultiQC to generate reports
        
        Args:
            input_dir: Input directory (default: trimmed_dir)
            
        Returns:
            MultiQC report path
        """
        logger.info("Starting MultiQC to generate reports...")
        
        if input_dir is None:
            input_dir = str(self.trimmed_dir)
        
        # Using Python API
        try:
            logger.info("Generating report using MultiQC Python API...")
            multiqc_output = self.reports_dir / "multiqc_report"
            
            report = multiqc.run(
                analysis_dir=input_dir,
                output_dir=str(self.reports_dir),
                filename="multiqc_report.html",
                template="default",
                force=True,
                quiet=False
            )
            
            logger.info(f"MultiQC report generation complete: {multiqc_output}.html")
            return multiqc_output
            
        except Exception as e:
            logger.warning(f"Failed using Python API, trying command line: {str(e)}")
            
            # Using command line
            try:
                cmd = [
                    self.multiqc_path,
                    input_dir,
                    '-o', str(self.reports_dir),
                    '-f',
                    '-v'
                ]
                
                logger.info(f"Running command: {' '.join(cmd)}")
                result = subprocess.run(cmd, check=True, capture_output=True, text=True)
                
                multiqc_output = self.reports_dir / "multiqc_report.html"
                logger.info(f"MultiQC report generation complete: {multiqc_output}")
                logger.debug(f"Output:\n{result.stdout}")
                
                return multiqc_output
                
            except subprocess.CalledProcessError as e:
                logger.error(f"MultiQC run failed: {str(e)}")
                logger.error(f"Error output:\n{e.stderr}")
                raise
    
    def get_qc_summary(self) -> pd.DataFrame:
        """
        Get QC summary statistics
        
        Returns:
            A DataFrame containing QC statistics
        """
        summary_data = []
        
        # Find report files generated by trim_galore
        for report_file in self.trimmed_dir.glob("*trimming_report.txt"):
            sample_name = report_file.stem.replace('_trimming_report', '')
            
            try:
                with open(report_file, 'r') as f:
                    content = f.read()
                    
                    stats = {
                        'Sample': sample_name,
                        'File': report_file.name
                    }
                    
                    lines = content.split('\n')
                    for line in lines:
                        if 'reads processed' in line:
                            stats['Total_reads'] = line.split(':')[1].strip()
                        elif 'reads with adapters' in line:
                            stats['Reads_with_adapters'] = line.split(':')[1].strip().split()[0]
                        elif 'reads written' in line:
                            stats['Reads_passed'] = line.split(':')[1].strip().split()[0]
                        elif 'Total basepairs processed' in line:
                            stats['Total_bp'] = line.split(':')[1].strip()
                        elif 'Quality-trimmed' in line:
                            stats['Quality_trimmed_bp'] = line.split(':')[1].strip()
                            
                    summary_data.append(stats)
                    
            except Exception as e:
                logger.warning(f"Failed to parse report file {report_file}: {str(e)}")
        
        if summary_data:
            return pd.DataFrame(summary_data)
        else:
            logger.warning("No QC report files found")
            return pd.DataFrame()
    
    def run_pipeline(self, fastq_list_file: str) -> Dict:
        """
        Run a complete single-end QC pipeline
        
        Args:
            fastq_list_file: List file containing fastq file paths
            
        Returns:
            A dictionary containing processing results
        """
        logger.info(f"Starting single-end QC pipeline, input file: {fastq_list_file}")
        
        try:
            fastq_files = self.read_fastq_list(fastq_list_file)
            trimmed_files = self.run_trim_galore(fastq_files)
            multiqc_report = self.run_multiqc()
            qc_summary = self.get_qc_summary()
            
            if not qc_summary.empty:
                summary_file = self.reports_dir / "qc_summary.csv"
                qc_summary.to_csv(summary_file, index=False)
                logger.info(f"CSV summary saved: {summary_file}")

            unsupported_summary = pd.DataFrame()
            if self.drop_unsupported and self.unsupported_stats:
                unsupported_summary = pd.DataFrame(self.unsupported_stats)
                unsupported_file = self.reports_dir / "unsupported_reads_summary.csv"
                unsupported_summary.to_csv(unsupported_file, index=False)
                logger.info(f"Unsupported-endpoint filter summary saved: {unsupported_file}")

            result = {
                'status': 'success',
                'trimmed_files': trimmed_files,
                'multiqc_report': multiqc_report,
                'qc_summary': qc_summary,
                'unsupported_summary': unsupported_summary,
                'output_dir': str(self.output_dir)
            }
            
            logger.info("Single-end QC pipeline completed!")
            return result
            
        except Exception as e:
            logger.error(f"Single-end QC pipeline failed: {str(e)}")
            return {
                'status': 'failed',
                'error': str(e)
            }


def prepare_input_list(input_path: str) -> str:
    """
    Helper function: Process user input. If it's a directory or a single file, 
    convert it to a temporary list file.
    """
    input_path = Path(input_path)
    tmp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False)
    tmp_file_path = tmp_file.name

    if input_path.is_file():
        if input_path.suffix in ['.txt', '.list']:
            # Already a list file, return original path
            tmp_file.close()
            os.remove(tmp_file_path)
            return str(input_path)
        elif input_path.suffix in ['.gz', '.fastq', '.fq']:
            # Is a single fastq file
            tmp_file.write(str(input_path.absolute()) + '\n')
            tmp_file.close()
            return tmp_file_path
    elif input_path.is_dir():
        # Is a directory, scan for fastq files inside
        fastq_extensions = ['.fastq.gz', '.fq.gz', '.fastq', '.fq']
        found_files = False
        for ext in fastq_extensions:
            for fastq_file in input_path.glob(f"*{ext}"):
                tmp_file.write(str(fastq_file.absolute()) + '\n')
                found_files = True
        tmp_file.close()
        
        if not found_files:
            os.remove(tmp_file_path)
            raise FileNotFoundError(f"No fastq files found in directory {input_path}")
            
        return tmp_file_path
    else:
        tmp_file.close()
        os.remove(tmp_file_path)
        raise FileNotFoundError(f"Input path does not exist: {input_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Single-end sequencing data QC processing tool (Trim Galore + MultiQC)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        '-i', '--input',
        type=str,
        required=True,
        help="Input path. Can be a .txt list file containing fastq paths, a directory containing fastq files, or a single fastq file."
    )
    parser.add_argument(
        '-o', '--output_dir',
        type=str,
        default="qc_results",
        help="Output directory path."
    )
    parser.add_argument(
        '-t', '--threads',
        type=int,
        default=4,
        help="Number of threads for Trim Galore."
    )
    parser.add_argument(
        '-q', '--quality',
        type=int,
        default=20,
        help="Quality trimming threshold."
    )
    parser.add_argument(
        '-l', '--length',
        type=int,
        default=18,
        help="Minimum read retention length after trimming."
    )
    parser.add_argument(
        '--trim_galore_path',
        type=str,
        default="trim_galore",
        help="Path to trim_galore executable."
    )
    parser.add_argument(
        '--multiqc_path',
        type=str,
        default="multiqc",
        help="Path to multiqc executable."
    )
    parser.add_argument(
        '--adapter',
        type=str,
        default="",
        help="Adapter sequence for trimming. If empty, Trim Galore auto-detects."
    )
    parser.add_argument(
        '--drop-unsupported',
        action='store_true',
        default=False,
        help="Optional (default OFF): remove trimmed reads that still have full "
             "cycle length L, i.e. reads with no detectable 3' adapter. Per-reviewer "
             "1 comment 3. Per-library L is auto-detected unless --cycle-length is set."
    )
    parser.add_argument(
        '--cycle-length',
        type=int,
        default=None,
        help="Fixed sequencing cycle count L for all libraries. Default: infer L "
             "per library from the raw FASTQ. Only used with --drop-unsupported."
    )

    args = parser.parse_args()

    try:
        # 1. Preprocess input path, convert to list file
        list_file = prepare_input_list(args.input)
        
        # 2. Initialize processor
        qc = QCProcessor(
            output_dir=args.output_dir,
            trim_galore_path=args.trim_galore_path,
            multiqc_path=args.multiqc_path,
            threads=args.threads,
            quality=args.quality,
            length=args.length,
            adapter=args.adapter,
            drop_unsupported=args.drop_unsupported,
            cycle_length=args.cycle_length
        )
        
        # 3. Run pipeline
        result = qc.run_pipeline(list_file)
        
        # 4. Output results
        if result['status'] == 'success':
            print("\n" + "="*40)
            print("QC processing completed successfully!")
            print(f"Output directory: {result['output_dir']}")
            print(f"MultiQC Report: {result['multiqc_report']}.html")
            
            if not result['qc_summary'].empty:
                print("\nQC Summary Statistics:")
                print(result['qc_summary'].to_string())

            unsupported = result.get('unsupported_summary')
            if unsupported is not None and not unsupported.empty:
                print("\nAdapter-unsupported reads removed:")
                print(unsupported[['Sample', 'cycle_length', 'total_trimmed_reads',
                                   'removed_reads', 'removed_pct']].to_string(index=False))
            print("="*40)
        else:
            print(f"\nQC processing failed: {result['error']}")
            exit(1)
            
    except Exception as e:
        print(f"\nAn error occurred: {str(e)}")
        exit(1)
    finally:
        if 'list_file' in locals() and os.path.abspath(list_file) != os.path.abspath(args.input) and os.path.exists(list_file):
            os.remove(list_file)


if __name__ == "__main__":
    main()
