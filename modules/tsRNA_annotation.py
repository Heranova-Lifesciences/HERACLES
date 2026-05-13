#!/usr/bin/env python3
"""
tsRNA Analysis Module
"""

import os
import sys
import logging
import subprocess
import argparse
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
from collections import defaultdict
import pandas as pd
import re

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class SimpleTsRNAAnalyzer:
    
    def __init__(
        self,
        index_dir: str,
        output_dir: str = "tsRNA_results",
        min_length: int = 15,
        max_length: int = 50,
        mismatch: int = 0,
        threads: int = 4,
        bowtie_path: str = "bowtie",
        keep_temp: bool = False
    ):
        self.index_dir = Path(index_dir)
        self.output_dir = Path(output_dir)
        self.min_length = min_length
        self.max_length = max_length
        self.mismatch = mismatch
        self.threads = threads
        self.bowtie_path = bowtie_path
        self.keep_temp = keep_temp
        
        # Find index files
        self.mature_index = self._find_bowtie_index()
        self.mature_fasta = self._find_fasta_file()
        
        self._validate_files()
        
        self.tRNA_seqs, self.tRNA_lengths = self._load_tRNA_seqs(self.mature_fasta)
        
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"SimpleTsRNAAnalyzer initialized")
        logger.info(f"Loaded {len(self.tRNA_seqs)} tRNA sequences")
        logger.info(f"Output directory: {self.output_dir}")
    
    def _find_bowtie_index(self) -> Path:
        """Find bowtie index file"""
        patterns = ["mature", "*.1.ebwt", "hg38-tRNA"]
        
        for pattern in patterns:
            if pattern.endswith(".ebwt"):
                ebwt_files = list(self.index_dir.glob(pattern))
                if ebwt_files:
                    index_base = ebwt_files[0].name.replace(".1.ebwt", "")
                    return self.index_dir / index_base
            else:
                if (self.index_dir / f"{pattern}.1.ebwt").exists():
                    return self.index_dir / pattern
        return self.index_dir / "mature"
    
    def _find_fasta_file(self) -> Path:
        """Find tRNA FASTA file"""
        # FIX: Be more specific to avoid picking up random .fa files like collapsed.fa
        fasta_patterns = ["hg38-tRNA.fa", "tRNA.fa", "mature.fa"]
        
        for pattern in fasta_patterns:
            if (self.index_dir / pattern).exists():
                return self.index_dir / pattern
                
        logger.warning("Could not find specific tRNA fasta, falling back to first .fa file. Ensure this is correct!")
        fasta_files = [f for f in self.index_dir.glob("*.fa") if 'collapsed' not in f.name]
        if fasta_files:
            return fasta_files[0]
            
        return self.index_dir / "hg38-tRNA.fa"
    
    def _validate_files(self):
        """Validate required files"""
        index_extensions = ['.1.ebwt', '.2.ebwt', '.3.ebwt', '.4.ebwt', '.rev.1.ebwt', '.rev.2.ebwt']
        index_found = any((self.mature_index.parent / f"{self.mature_index.name}{ext}").exists() for ext in index_extensions)
        
        if not index_found:
            logger.error(f"Bowtie index files not found for {self.mature_index}")
            raise FileNotFoundError(f"Bowtie index not found in {self.index_dir}")
        
        if not self.mature_fasta.exists():
            logger.error(f"tRNA FASTA file not found: {self.mature_fasta}")
            raise FileNotFoundError(f"tRNA FASTA file not found in {self.index_dir}")
        
        logger.info(f"Found bowtie index: {self.mature_index}")
        logger.info(f"Found tRNA FASTA: {self.mature_fasta}")

    def _load_tRNA_seqs(self, fasta_path: Path) -> Tuple[Dict[str, str], Dict[str, int]]:
        """Load tRNA sequences and lengths"""
        sequences, lengths = {}, {}
        try:
            with open(fasta_path, 'r') as f:
                current_id, current_seq = None, []
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        if current_id:
                            full_seq = ''.join(current_seq)
                            sequences[current_id], lengths[current_id] = full_seq, len(full_seq)
                        current_id = line[1:].split()[0]
                        current_seq = []
                    else:
                        current_seq.append(line)
                if current_id:
                    full_seq = ''.join(current_seq)
                    sequences[current_id], lengths[current_id] = full_seq, len(full_seq)
            return sequences, lengths
        except Exception as e:
            logger.error(f"Failed to load tRNA sequences: {e}")
            return {}, {}
    
    def _parse_read_count_from_header(self, header: str) -> Tuple[str, int]:
        """Parse read count from FASTA header"""
        if header.startswith('>'): header = header[1:]
        
        parts = header.split('_')
        if len(parts) >= 2 and parts[-1].isdigit():
            return '_'.join(parts[:-1]), int(parts[-1])
            
        if '_x' in header:
            x_parts = header.split('_x')
            if len(x_parts) == 2 and x_parts[1].isdigit():
                return x_parts[0], int(x_parts[1])
                
        return header, 1
    
    def _filter_fasta_by_length(self, input_fasta: Path, output_fasta: Path) -> int:
        """Filter FASTA file by length. FIX: DO NOT EXPAND SEQUENCES"""
        kept_seqs = 0
        try:
            with open(input_fasta, 'r') as fin, open(output_fasta, 'w') as fout:
                current_header, current_seq = None, []
                for line in fin:
                    line = line.strip()
                    if line.startswith('>'):
                        if current_header is not None:
                            seq = ''.join(current_seq)
                            if self.min_length <= len(seq) <= self.max_length:
                                fout.write(f">{current_header}\n{seq}\n")
                                kept_seqs += 1
                        current_header = line[1:]
                        current_seq = []
                    else:
                        current_seq.append(line)
                
                if current_header is not None:
                    seq = ''.join(current_seq)
                    if self.min_length <= len(seq) <= self.max_length:
                        fout.write(f">{current_header}\n{seq}\n")
                        kept_seqs += 1
            
            logger.info(f"Length filtering: {kept_seqs} unique sequences kept (unexpanded)")
            return kept_seqs
        except Exception as e:
            logger.error(f"Length filtering failed: {e}")
            return 0
    
    def _run_bowtie_simple(self, input_fasta: Path, sam_output: Path) -> Tuple[int, Path]:
        """Run bowtie alignment"""
        temp_fasta = sam_output.with_suffix('.temp.fasta')
        
        filtered_count = self._filter_fasta_by_length(input_fasta, temp_fasta)
        if filtered_count == 0:
            logger.warning("No sequences passed length filtering")
            if not self.keep_temp and temp_fasta.exists(): temp_fasta.unlink() # FIX: cleanup on early return
            return 0, sam_output
            
        cmd = [
            self.bowtie_path, "-p", str(self.threads), str(self.mature_index),
            "-f", str(temp_fasta), "-v", str(self.mismatch),
            "-a", "--best", "--strata", "-S", "--norc", str(sam_output)
        ]
        
        logger.info(f"Running bowtie: {' '.join(cmd)}")
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True, encoding='utf-8', timeout=3600)
            mapped_reads = 0
            for line in result.stderr.split('\n'):
                if '# reads with at least one alignment:' in line:
                    mapped_reads = int(line.split(':')[1].strip().split()[0])
            logger.info(f"Bowtie alignment: {mapped_reads}/{filtered_count} unique sequences aligned")
            
            if not self.keep_temp and temp_fasta.exists(): temp_fasta.unlink()
            return mapped_reads, sam_output
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Bowtie failed: {e}\n{e.stderr[:500]}")
            if not self.keep_temp and temp_fasta.exists(): temp_fasta.unlink()
            return 0, sam_output
    
    def _parse_sam_file(self, sam_file: Path) -> Dict[str, Dict[str, Any]]:
        """
        Parse SAM file with proper multi-mapping resolution per read.
        """
        read_alignments = defaultdict(list)

        if not sam_file.exists():
            return {}

        total_lines = 0
        skipped_unmapped = 0
        skipped_unknown_ref = 0
        skipped_short_cigar = 0

        try:
            with open(sam_file, 'r') as f:
                for line in f:
                    if line.startswith('@'):
                        continue
                    total_lines += 1
                    parts = line.strip().split('\t')
                    if len(parts) < 10:
                        continue

                    read_id, flag, ref_name, pos, cigar = parts[0], int(parts[1]), parts[2], int(parts[3]), parts[5]
                    if flag & 0x4 or ref_name == "*":
                        skipped_unmapped += 1
                        continue
                    if ref_name not in self.tRNA_lengths:
                        skipped_unknown_ref += 1
                        continue

                    # Calculate alignment length (read length on reference)
                    align_len = 0
                    if cigar != "*":
                        i = 0
                        while i < len(cigar):
                            j = i
                            while j < len(cigar) and cigar[j].isdigit():
                                j += 1
                            if j > i:
                                op = cigar[j]
                                if op in 'M=X':
                                    align_len += int(cigar[i:j])
                                i = j + 1
                            else:
                                i += 1

                    if align_len == 0:
                        skipped_short_cigar += 1
                        continue

                    end_pos = pos + align_len - 1
                    tRNA_len = self.tRNA_lengths[ref_name]
                    tsRNA_type = self._classify_tsRNA(ref_name, pos, end_pos, align_len, tRNA_len)

                    seq_fragment = ""
                    if pos >= 1 and end_pos <= tRNA_len:
                        seq_fragment = self.tRNA_seqs[ref_name][pos - 1:end_pos]

                    read_alignments[read_id].append({
                        'tRNA': ref_name,
                        'start': pos,
                        'end': end_pos,
                        'length': align_len,
                        'type': tsRNA_type,
                        'sequence': seq_fragment,
                        'cigar': cigar
                    })
            
            logger.info(
                f"SAM stats: {total_lines} alignments, "
                f"{skipped_unmapped} unmapped, "
                f"{skipped_unknown_ref} unknown ref, "
                f"{skipped_short_cigar} zero-length, "
                f"{sum(len(v) for v in read_alignments.values())} kept "
                f"(from {len(read_alignments)} unique reads)"
            )

            if skipped_unknown_ref > 0 and skipped_unknown_ref > total_lines * 0.5:
                logger.warning(
                    f"{skipped_unknown_ref}/{total_lines} alignments have unknown reference names. "
                    f"Check that the FASTA file matches the Bowtie index. "
                    f"Known refs: {len(self.tRNA_lengths)}, first few: "
                    f"{list(self.tRNA_lengths.keys())[:5]}"
                )

            # 2. Resolve multi-mappings and aggregate counts
            tsRNA_dict = {}
            
            for read_id, alignments in read_alignments.items():
                # Get the real count from the header (since we didn't expand)
                _, count = self._parse_read_count_from_header(read_id)
                
                if len(alignments) == 1:
                    best_hit = alignments[0]
                else:
                    # Multi-mapping! Apply XY Minimum Principle to pick ONE best hit
                    best_hit = min(alignments, key=lambda x: self._get_xy_sort_key(x['tRNA']))
                
                # Create unique key for the final dictionary
                key = f"{best_hit['tRNA']}:{best_hit['start']}-{best_hit['end']}:{best_hit['type']}"
                
                if key not in tsRNA_dict:
                    tsRNA_dict[key] = {
                        'tRNA': best_hit['tRNA'],
                        'start': best_hit['start'],
                        'end': best_hit['end'],
                        'length': best_hit['length'],
                        'type': best_hit['type'],
                        'sequence': best_hit['sequence'],
                        'cigar': best_hit['cigar'],
                        'count': 0
                    }
                
                # Add the REAL weight (count) ONLY to the chosen best hit
                tsRNA_dict[key]['count'] += count
                
            logger.info(f"Parsed SAM: {len(tsRNA_dict)} unique tsRNA loci identified (multi-mappings resolved)")
            return tsRNA_dict
            
        except Exception as e:
            logger.error(f"Failed to parse SAM file: {e}")
            return {}
    
    def _classify_tsRNA(self, tRNA_name, start_pos, end_pos, fragment_length, tRNA_length):
        starts_at_5prime = (start_pos == 1)
        ends_at_3prime = (end_pos == tRNA_length)
        
        if starts_at_5prime:
            if fragment_length <= 30: return "tRF-5"
            elif 31 <= fragment_length <= 40: return "tiRNA-5"
            else: return "tiRNA-5L"
        elif ends_at_3prime:
            if fragment_length <= 30: return "tRF-3"
            elif 31 <= fragment_length <= 40: return "tiRNA-3"
            else: return "tiRNA-3L"
        return "tRF-i"
    
    def _get_xy_sort_key(self, tRNA_name: str):
        try:
            parts = tRNA_name.split('-')
            if len(parts) >= 5:
                return (int(parts[-2]), int(parts[-1]))
        except ValueError: pass
        return (float('inf'), float('inf'))

    def analyze_sample(self, sample_name, fasta_file, output_prefix=None):
        logger.info(f"Starting analysis for sample: {sample_name}")
        if not fasta_file.exists():
            return {'status': 'error', 'message': f'FASTA file does not exist: {fasta_file}'}
        
        if output_prefix is None: output_prefix = sample_name
        
        sam_file = self.output_dir / f"{output_prefix}_aligned.sam"
        annotation_file = self.output_dir / f"{output_prefix}_annotation.tsv"
        count_file = self.output_dir / f"{output_prefix}_counts.tsv"
        
        try:
            mapped_count, sam_path = self._run_bowtie_simple(fasta_file, sam_file)
            
            if mapped_count == 0:
                logger.warning("No reads aligned successfully, trying with 1 mismatch...")
                original_mismatch = self.mismatch
                self.mismatch = 1
                mapped_count, sam_path = self._run_bowtie_simple(fasta_file, sam_file)
                self.mismatch = original_mismatch
                if mapped_count == 0:
                    return {'status': 'warning', 'message': 'No reads aligned successfully'}
            
            tsRNA_results = self._parse_sam_file(sam_path)
            if not tsRNA_results:
                return {'status': 'warning', 'message': 'No tsRNAs identified'}
            
            # 1. Generate Annotation
            annotation_data = []
            for key, info in tsRNA_results.items():
                annotation_data.append({
                    'tsRNA_id': key, 'tRNA': info['tRNA'], 'start': info['start'],
                    'end': info['end'], 'length': info['length'], 'type': info['type'],
                    'sequence': info['sequence'], 'count': info['count'], 'cigar': info['cigar']
                })
            
            annotation_df = pd.DataFrame(annotation_data)[['tsRNA_id', 'tRNA', 'start', 'end', 'length', 'type', 'sequence', 'count', 'cigar']]
            annotation_df.to_csv(annotation_file, sep='\t', index=False)
            
            # 2. XY Deduplication: same type+sequence mapping to different tRNA loci
            #    (tRNA gene duplication produces identical sequence fragments at different genomic loci)
            groups = defaultdict(list)
            empty_seq_count = 0
            for key, info in tsRNA_results.items():
                seq = info['sequence']
                if not seq:
                    empty_seq_count += 1
                    group_key = key
                else:
                    group_key = f"{info['type']}:{seq}"
                groups[group_key].append((key, info))

            final_counts = {}
            dedup_groups = 0
            for group_key, items in groups.items():
                if len(items) == 1:
                    orig_key, info = items[0]
                    final_counts[f"{orig_key}:{info['sequence']}"] = info['count']
                else:
                    dedup_groups += 1
                    winner = min(items, key=lambda x: self._get_xy_sort_key(x[1]['tRNA']))
                    total_count = sum(item[1]['count'] for item in items)
                    final_counts[f"{winner[0]}:{winner[1]['sequence']}"] = total_count

            logger.info(
                f"XY dedup: {len(tsRNA_results)} loci → {len(final_counts)} features "
                f"({dedup_groups} groups collapsed, {empty_seq_count} empty-seq skipped)"
            )
            
            count_df = pd.DataFrame.from_dict(final_counts, orient='index', columns=['count'])
            count_df.index.name = 'tsRNA_id'
            count_df = count_df.reset_index().sort_values(by=['count'], ascending=False).set_index('tsRNA_id')
            count_df.to_csv(count_file, sep='\t')
            
            # Cleanup
            if not self.keep_temp:
                for temp_file in self.output_dir.glob(f"{output_prefix}*.temp.*"):
                    if temp_file.exists(): temp_file.unlink()
            
            # Stats
            total_tsRNAs = len(count_df)
            total_reads = count_df['count'].sum()
            type_counts = {}
            for tsRNA_id in count_df.index:
                parts = tsRNA_id.split(':')
                if len(parts) >= 4:
                    ttype = parts[2]
                    type_counts[ttype] = type_counts.get(ttype, 0) + count_df.loc[tsRNA_id, 'count']
            
            logger.info(f"Analysis completed! Unique seqs: {total_tsRNAs}, Total reads: {total_reads}")
            for ttype, count in sorted(type_counts.items()):
                logger.info(f"  {ttype}: {count} reads")
            
            return {
                'status': 'success', 'sample_name': sample_name, 'tsRNA_count': total_tsRNAs,
                'total_reads': total_reads, 'type_counts': type_counts,
                'annotation_file': str(annotation_file), 'count_file': str(count_file)
            }
            
        except Exception as e:
            logger.error(f"Failed to analyze sample {sample_name}: {e}")
            import traceback; traceback.print_exc()
            return {'status': 'error', 'sample_name': sample_name, 'error': str(e)}


def main():
    parser = argparse.ArgumentParser(description='Simplified tsRNA analysis tool')
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('-i', '--input', help='Input FASTA file')
    input_group.add_argument('-s', '--samples', help='Sample list file (format: sample_name\\tfasta_path)')
    
    parser.add_argument('-d', '--index-dir', required=True, help='tRNA index directory')
    parser.add_argument('-n', '--sample-name', help='Sample name (used with -i)')
    parser.add_argument('-o', '--output-dir', default='tsRNA_results', help='Output directory')
    parser.add_argument('--min-len', type=int, default=18, help='Minimum length')
    parser.add_argument('--max-len', type=int, default=45, help='Maximum length')
    parser.add_argument('-v', '--mismatch', type=int, default=0, help='Mismatches allowed')
    parser.add_argument('-t', '--threads', type=int, default=4, help='Number of threads')
    parser.add_argument('--bowtie-path', default='bowtie', help='Path to bowtie')
    parser.add_argument('--keep-temp', action='store_true', help='Keep temporary files')
    parser.add_argument('--debug', action='store_true', help='Enable debug mode')
    
    args = parser.parse_args()
    if args.debug: logging.getLogger().setLevel(logging.DEBUG)
    
    try:
        analyzer = SimpleTsRNAAnalyzer(
            index_dir=args.index_dir, output_dir=args.output_dir,
            min_length=args.min_len, max_length=args.max_len,
            mismatch=args.mismatch, threads=args.threads,
            bowtie_path=args.bowtie_path, keep_temp=args.keep_temp
        )
        
        if args.input:
            sample_name = args.sample_name or Path(args.input).stem.replace('_collapsed', '')
            result = analyzer.analyze_sample(sample_name=sample_name, fasta_file=Path(args.input))
            if result['status'] == 'success':
                logger.info(f"Success! Annotation: {result['annotation_file']}, Counts: {result['count_file']}")
            else: sys.exit(1)
        
        elif args.samples:
            with open(args.samples, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line or line.startswith('#'): continue
                    parts = line.split()
                    s_name = parts[0] if len(parts) >= 2 else Path(parts[0]).stem.replace('_collapsed', '')
                    s_path = parts[1] if len(parts) >= 2 else parts[0]
                    
                    if not Path(s_path).exists():
                        logger.error(f"FASTA not found: {s_path}"); continue
                    
                    result = analyzer.analyze_sample(sample_name=s_name, fasta_file=Path(s_path))
            logger.info("All samples processing completed.")
    except Exception as e:
        logger.error(f"Program failed: {e}"); sys.exit(1)

if __name__ == "__main__":
    main()
