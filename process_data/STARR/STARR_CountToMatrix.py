#!/usr/bin/env python3
"""
Merge per-sample count files into a final count matrix.

Inputs:
  --dna_dir         directory containing folders DNA1, DNA2, ...
  --rna_dir         directory containing folders RNA1, RNA2, ...
  --count_filename  the per-sample count file inside each folder (e.g., counts.bed)
  --output_dir
  --count_type      fragment | bin
  --chrom_sizes     required for fragment mode
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
import re

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

class CountMatrixMerger:
    def __init__(self, dna_files, rna_files, sample_list, count_type, chrom_sizes=None):
        self.dna_files = dna_files
        self.rna_files = rna_files
        self.sample_list = sample_list  # e.g., ["DNA1","DNA2",...,"RNA1","RNA2",...]
        self.count_type = count_type
        self.chrom_sizes = chrom_sizes

    def get_chr_list(self):
        chr_ref = [x.split('\t')[0] for x in open(self.chrom_sizes, 'r').readlines()]
        all_other_chr = [x for x in chr_ref if '_' in x]
        chr_list = [[x] for x in chr_ref if x not in all_other_chr]
        chr_list = chr_list + [all_other_chr]
        chr_list = [y for x in chr_list for y in x]
        return chr_list

    def get_count_per_sample(self, file_path, chrom, sample):
        with open(file_path, 'r') as input_file:
            content = [line.strip().split('\t') for line in input_file if line.startswith(chrom)]

        if content:
            df = pd.DataFrame(content)
            # enforce numeric for start/end
            df = df.astype({1: int, 2: int})
            df = df.set_index([0, 1, 2, 3])
            df.columns = [sample]
            return df
        else:
            # empty frame with right index shape
            df = pd.DataFrame(columns=[sample])
            df.index = pd.MultiIndex.from_arrays([[], [], [], []], names=[0, 1, 2, 3])
            return df

    def write_fragment_matrix(self, output_dir):
        utils.set_dir(output_dir)
        chr_list = self.get_chr_list()
        all_files = self.dna_files + self.rna_files

        out_path = os.path.join(output_dir, 'count_mat.bed')
        if os.path.exists(out_path):
            os.remove(out_path)

        for chrom in chr_list:
            count_mat = None
            for i, file_path in enumerate(all_files):
                sample = self.sample_list[i]
                df = self.get_count_per_sample(file_path, chrom, sample)
                count_mat = df if count_mat is None else count_mat.join(df, how='outer')

            if count_mat is not None and not count_mat.empty:
                count_mat = count_mat[self.sample_list].fillna(0).reset_index()
                count_mat.to_csv(out_path, sep='\t', mode='a', index=False, header=False)

    def write_bin_matrix(self, output_dir):
        utils.set_dir(output_dir)
        file_list = self.dna_files + self.rna_files
        if not file_list:
            raise ValueError("No input files found for bin mode.")

        num_rows = utils.get_row_number(file_list[0])

        open_files = [open(f, 'r') for f in file_list]
        out_path = os.path.join(output_dir, 'count_mat_all_bins.bed')
        with open(out_path, 'w') as out_file:
            for _ in range(num_rows):
                line_parts = [f.readline().strip().split('\t') for f in open_files]
                base = line_parts[0][:-1]
                counts = [line[-1] for line in line_parts]
                # write only if at least one count > 0
                if any(float(c) > 0 for c in counts):
                    out_file.write('\t'.join(base + counts) + '\n')

        for f in open_files:
            f.close()

    def merge(self, output_dir):
        if self.count_type == 'fragment':
            assert self.chrom_sizes is not None, "chrom_sizes must be specified for fragment mode."
            self.write_fragment_matrix(output_dir)
        elif self.count_type == 'bin':
            self.write_bin_matrix(output_dir)
        else:
            raise ValueError("count_type must be either 'fragment' or 'bin'")

def _natural_index(name: str, expected_prefix: str):
    """
    Returns (number, name) for sorting things like DNA1, DNA2, ..., RNA1, ...
    Unmatched names go to the end but stay deterministic.
    """
    m = re.fullmatch(rf"{re.escape(expected_prefix)}(\d+)", name)
    if m:
        return (int(m.group(1)), name)
    return (10**9, name)

def collect_samples(root_dir: str, expected_prefix: str, count_filename: str):
    """
    Find subfolders like DNA1, DNA2, ... (or RNA1, RNA2, ...) and return:
        samples: [(sample_name, file_path), ...] sorted naturally by the trailing number.
    Ensures the count file exists inside each folder.
    """
    samples = []
    if not root_dir:
        return samples

    for entry in os.listdir(root_dir):
        subdir = os.path.join(root_dir, entry)
        if not os.path.isdir(subdir):
            continue
        if not entry.startswith(expected_prefix):
            continue

        candidate = os.path.join(subdir, count_filename)
        if os.path.isfile(candidate):
            samples.append((entry, candidate))
        else:
            # tolerate nested structures: search recursively for the exact basename
            matches = glob.glob(os.path.join(subdir, "**", count_filename), recursive=True)
            if matches:
                samples.append((entry, matches[0]))

    samples.sort(key=lambda x: _natural_index(x[0], expected_prefix))
    return samples

def main():
    parser = argparse.ArgumentParser(description="Merge per-sample STARR-seq counts into a matrix.")
    parser.add_argument('--dna_dir', required=True, help="Directory containing DNA1, DNA2, ... subfolders.")
    parser.add_argument('--rna_dir', required=True, help="Directory containing RNA1, RNA2, ... subfolders.")
    parser.add_argument('--count_filename', required=True, help="Basename of per-sample count file (inside each subfolder).")
    parser.add_argument('--output_dir', required=True, help="Directory to save merged count matrix.")
    parser.add_argument('--count_type', choices=['fragment', 'bin'], required=True, help="Type of count matrix.")
    parser.add_argument('--chrom_sizes', help="Required for fragment count mode.")
    args = parser.parse_args()

    dna_samples = collect_samples(args.dna_dir, "DNA", args.count_filename)
    rna_samples = collect_samples(args.rna_dir, "RNA", args.count_filename)

    if not dna_samples:
        raise FileNotFoundError(f"No samples found under DNA dir: {args.dna_dir}")
    if not rna_samples:
        raise FileNotFoundError(f"No samples found under RNA dir: {args.rna_dir}")

    # Build lists in deterministic order and use folder names as sample names
    dna_files = [p for (name, p) in dna_samples]
    rna_files = [p for (name, p) in rna_samples]
    sample_list = [name for (name, _) in dna_samples] + [name for (name, _) in rna_samples]

    utils.set_dir(args.output_dir)

    merger = CountMatrixMerger(
        dna_files=dna_files,
        rna_files=rna_files,
        sample_list=sample_list,
        count_type=args.count_type,
        chrom_sizes=args.chrom_sizes
    )
    merger.merge(args.output_dir)

if __name__ == '__main__':
    main()