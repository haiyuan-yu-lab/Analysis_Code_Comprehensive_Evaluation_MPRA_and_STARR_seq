#!/usr/bin/env python3
"""
Script to generate count matrices from STARR-seq BAM files.

Supports two modes:
1. Full mode (default): Process BAM → generate alignment BED → count.
2. Count-only mode (--count_only): Skip BAM processing and directly generate counts from existing processed alignments.
"""

# --------------------------------------------------
# Standard Library Imports
# --------------------------------------------------
import os
import sys
import glob
import pandas as pd
import argparse
import pysam
import pybedtools
import shutil

from subprocess import call, PIPE, run, Popen, STDOUT
from multiprocessing import Pool, cpu_count

# --------------------------------------------------
# Local Module Imports
# --------------------------------------------------
project_root = os.path.abspath(os.path.join(os.getcwd(), '..', '..'))
src_path = os.path.join(project_root, 'src')
sys.path.insert(0, src_path)

import utils
import starrseq_bam_processing

# --------------------------------------------------
# Utility Functions
# --------------------------------------------------

def check_and_index_bam(bam_file):
    """
    Ensure BAM index (.bai) file exists. Create one if missing.
    """
    bai_file = bam_file.split('.')[0] + '.bai'
    if not os.path.exists(bai_file):
        print(f"Creating index for {bam_file}")
        pysam.index(bam_file)
    else:
        print(f"Index exists: {bai_file}")

import os
import pybedtools
import utils  # assumes your utils has safe_remove and safe_bedsort_optimized

def make_bin(window_size,
             step_size,
             chrom_sizes_path,
             out_dir,
             tmp_dir=None,
             n_cpu=1,
             mem_per_thread="5G"):
    """
    Generate genome-wide sliding windows and write a sorted BED.

    Args:
        window_size (int): window length (e.g., 100)
        step_size (int): step size (e.g., 10)
        chrom_sizes_path (str): path to UCSC-style chrom sizes
        out_dir (str): directory to write outputs
        tmp_dir (str|None): tmp dir for sorting (defaults to $TMPDIR or /tmp inside sorter)
        n_cpu (int): CPUs available; sorter will use about half
        mem_per_thread (str): GNU sort -S per-thread memory (e.g., "5G")

    Returns:
        str: path to the sorted windows BED
    """
    os.makedirs(out_dir, exist_ok=True)

    # 1) Make windows via pybedtools
    windows_bt = pybedtools.BedTool().window_maker(g=chrom_sizes_path,
                                                   w=int(window_size),
                                                   s=int(step_size))
    windows_df = windows_bt.to_dataframe(disable_auto_names=True, header=None)

    # 2) Write unsorted temp BED
    tmp_bed = os.path.join(out_dir, f"{window_size}_{step_size}_tmp.windows.bed")
    windows_df.to_csv(tmp_bed, sep="\t", index=False, header=False)
    del windows_df

    # 3) Sort to final BED
    sorted_bed = os.path.join(out_dir, f"{window_size}_{step_size}_tmp.windows.sorted.bed")
    parallel = max(1, int(n_cpu // 2))  # use ~half of available CPUs

    # Prefer the optimized sorter if available
    if hasattr(utils, "safe_bedsort_optimized"):
        utils.safe_bedsort_optimized(
            tmp_bed,
            sorted_bed,
            tmpdir=tmp_dir,
            parallel=parallel,
            mem=mem_per_thread
        )
    else:
        # Fallback to the original simple sorter
        utils.safe_bedsort(tmp_bed, sorted_bed)

    # 4) Clean up temp
    utils.safe_remove(tmp_bed)


def process_single_bam(args_tuple):
    """
    Process BAM file to generate alignment BED.

    Args:
        args_tuple: (bam_file, args, lib_type, tmp_bedtools_dir, use_multiprocessing)

    Returns:
        str: Processing status message.
    """
    bam_file, args, lib_type, tmp_bedtools_dir, use_multiprocessing = args_tuple

    sample_name = os.path.basename(bam_file).replace('.bam', '')
    sample_output_dir = os.path.join(args.output_dir, lib_type, sample_name)
    os.makedirs(sample_output_dir, exist_ok=True)

    processor = starrseq_bam_processing.ProcessBamSTARR(
        bam_dir=bam_file,
        min_frag_size=args.min_frag_size,
        max_frag_size=args.max_frag_size,
        prefix=sample_output_dir + '/',
        tmp_dir=tmp_bedtools_dir,
        n_cpu=args.n_cpu
    )

    processor.process_bam(
        chrom_size_dir=args.chrom_sizes,
        remove_duplicate=True,
        use_multiprocessing=use_multiprocessing
    )

    return f"Processed {sample_name} ({lib_type})"

def generate_single_count(args_tuple):
    """
    Generate fragment or bin-level count from pre-processed alignments.

    Args:
        args_tuple: (bam_file, args, lib_type, tmp_bedtools_dir)

    Returns:
        str: Counting status message.
    """
    bam_file, args, lib_type, tmp_bedtools_dir = args_tuple

    sample_name = os.path.basename(bam_file).replace('.bam', '')
    sample_output_dir = os.path.join(args.output_dir, lib_type, sample_name)

    processor = starrseq_bam_processing.ProcessBamSTARR(
        bam_dir=bam_file,
        min_frag_size=args.min_frag_size,
        max_frag_size=args.max_frag_size,
        prefix=sample_output_dir + '/',
        tmp_dir=tmp_bedtools_dir,
        n_cpu=args.n_cpu
    )

    if args.count_type == 'bin':
        # Determine whether to collapse based on UMI usage
        if args.umi_usage == 'Both':
            collapse = False
        elif args.umi_usage == 'DNA':
            collapse = lib_type != 'DNA'
        elif args.umi_usage == 'RNA':
            collapse = lib_type != 'RNA'
        else:  # 'None'
            collapse = True

        processor.get_genomic_bin_count(
            windowSize=args.bin_size,
            stepSize=args.step_size,
            collapse=collapse,
            chrom_size_dir=args.chrom_sizes,
            tmp_bin_path=args.output_dir
        )

    elif args.count_type == 'fragment':
        processor.get_frag_count(collapse=False)

    return f"Generated count for {sample_name} ({lib_type})"

# --------------------------------------------------
# Main
# --------------------------------------------------

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-d', '--data_dir',
        help="Directory containing 'DNA' and 'RNA' subfolders with BAM files."
    )
    parser.add_argument(
        '-o', '--output_dir',
        help="Directory to store intermediate and final outputs."
    )
    parser.add_argument('--min_frag_size', type=int, default=150, help="Minimum fragment size to include (default: 150).")
    parser.add_argument('--max_frag_size', type=int, default=800, help="Maximum fragment size to include (default: 800).")
    parser.add_argument('--bin_size', type=int, default=100, help="Genomic bin size for counting (default: 100 bp).")
    parser.add_argument('--step_size', type=int, default=10, help="Step size for sliding window (default: 10 bp).")
    parser.add_argument('--count_type', choices=['bin', 'fragment'], default='bin', help="Count by genomic bins or by fragments (default: bin).")
    parser.add_argument('--chrom_sizes', help="Path to chromosome size file.")
    parser.add_argument('--umi_usage', choices=['Both', 'DNA', 'RNA', 'None'], default='None', help="Specify UMI usage (default: None).")
    parser.add_argument('--n_cpu', type=int, default=4, help="Number of CPUs to use (default: 4).")
    parser.add_argument('--count_only', action='store_true', default=False, help="Skip BAM processing and only generate counts (default: False).")

    args = parser.parse_args()

    # Validate input directory
    if not os.path.isdir(args.data_dir):
        sys.exit(f"Error: Data directory not found: {args.data_dir}")

    dna_dir = os.path.join(args.data_dir, 'DNA')
    rna_dir = os.path.join(args.data_dir, 'RNA')

    if not os.path.isdir(dna_dir) or not os.path.isdir(rna_dir):
        sys.exit(f"Error: Missing 'DNA' or 'RNA' subdirectory in {args.data_dir}")

    print("DNA and RNA directories found. Proceeding...")

    # Gather all BAM files
    dna_bam_files = sorted(glob.glob(os.path.join(dna_dir, '*.bam')))
    rna_bam_files = sorted(glob.glob(os.path.join(rna_dir, '*.bam')))

    print(f"Found {len(dna_bam_files)} DNA BAM files.")
    print(f"Found {len(rna_bam_files)} RNA BAM files.")

    for bam_file in dna_bam_files + rna_bam_files:
        check_and_index_bam(bam_file)

    # Prepare output and temp directories
    utils.set_dir(args.output_dir)
    print(f"Output directory: {args.output_dir}")

    tmp_bedtools_dir = os.path.join(args.output_dir, 'tmp_bedtools')
    utils.set_dir(tmp_bedtools_dir)
    pybedtools.helpers.set_tempdir(tmp_bedtools_dir)
    print(f"pybedtools temp dir: {tmp_bedtools_dir}")

    # Prepare job arguments
    if not args.count_only:
        print("Running full BAM processing ...")
        args_list = [(bam, args, 'DNA', tmp_bedtools_dir, False) for bam in dna_bam_files]
        args_list += [(bam, args, 'RNA', tmp_bedtools_dir, False) for bam in rna_bam_files]

        with Pool(processes=min(args.n_cpu, len(args_list))) as pool:
            results = pool.map(process_single_bam, args_list)

        print("Running count generation ... ")
        # Make bin
        print("Make bin bed file first ...")
        make_bin(window_size=args.bin_size,
                 step_size=args.step_size,
                 chrom_sizes_path=args.chrom_sizes,
                 out_dir=args.output_dir,
                 tmp_dir=tmp_bedtools_dir,
                 n_cpu=args.n_cpu,
                 mem_per_thread="5G")
        
        # Run to obtain count matrix
        args_list = [(bam, args, 'DNA', tmp_bedtools_dir) for bam in dna_bam_files]
        args_list += [(bam, args, 'RNA', tmp_bedtools_dir) for bam in rna_bam_files]

        with Pool(processes=min(args.n_cpu, len(args_list))) as pool:
            results = pool.map(generate_single_count, args_list)

    else:
        print("Running in count-only mode...")
        args_list = [(bam, args, 'DNA', tmp_bedtools_dir) for bam in dna_bam_files]
        args_list += [(bam, args, 'RNA', tmp_bedtools_dir) for bam in rna_bam_files]

        with Pool(processes=min(args.n_cpu, len(args_list))) as pool:
            results = pool.map(generate_single_count, args_list)

    for res in results:
        print(res)

    print("All tasks completed.")

    # --------------------------------------------------
    # Remove pybedtools temporary directory
    # --------------------------------------------------
    try:
        shutil.rmtree(tmp_bedtools_dir)
        print(f"Temporary directory removed: {tmp_bedtools_dir}")
    except Exception as e:
        print(f"Warning: Failed to remove temporary directory: {tmp_bedtools_dir}")
        print(f"Reason: {e}")

if __name__ == '__main__':
    main()