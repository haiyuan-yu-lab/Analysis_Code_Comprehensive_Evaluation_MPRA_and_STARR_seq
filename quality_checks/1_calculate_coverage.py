#!/usr/bin/env python3
"""
Script for processing STARR-seq count matrices and performing coverage analysis.

If raw (unfiltered) data is provided, the script filters regions based on 
(1) total DNA read count thresholds,
(2) DNA read count thresholds across each DNA sample. 
The filtered data will be saved and used to generate statistics of genome-wide coverage and user-specified region coverage. 

Need to specify the directory to import functions in utils.py, and quality_control.py. 

"""

# --------------------------------------------------
# Standard Library Imports
# --------------------------------------------------
import os
import sys
import glob
import math
from subprocess import call, PIPE, run, Popen, STDOUT
from multiprocessing import Pool, cpu_count
from random import sample

# --------------------------------------------------
# Third-Party Library Imports
# --------------------------------------------------
import numpy as np
import pandas as pd
import argparse
import pybedtools
import pysam
import seaborn as sns
from matplotlib import pyplot as plt
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.axes_grid1 import make_axes_locatable  # Used for advanced subplot layouts
import matplotlib.patches as mpatches
import matplotlib.ticker as mtick
from matplotlib.offsetbox import AnchoredText
import matplotlib.colors as clr
from matplotlib import cm
from matplotlib.colors import Normalize 
from scipy.interpolate import interpn
import scipy
import statsmodels.stats.multitest as smm
from Bio import SeqIO

# --------------------------------------------------
# Matplotlib Configuration
# --------------------------------------------------
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Helvetica Neue'

# --------------------------------------------------
# Local Module Imports and Configuration
# --------------------------------------------------
# Dynamically add src directory to sys.path
project_root = os.path.abspath(os.path.join(os.getcwd(), '..'))
src_path = os.path.join(project_root, 'src')
sys.path.insert(0, src_path)

import utils
from quality_control import CheckCoverage


# Set the project root directory (assuming current working directory is within a subdirectory)
project_root = os.path.abspath(os.path.join(os.getcwd(), '..'))

# Set the temporary directory for pybedtools.
pybedtools.helpers.set_tempdir('/fs/cbsuhy02/storage/jz855/tmp/')

# --------------------------------------------------
# Command-Line Argument Parsing
# --------------------------------------------------
parser = argparse.ArgumentParser(
    description="Filter STARR-seq count matrices and perform coverage analysis."
)
parser.add_argument(
    '-c', "--count_mat", 
    help="Input count matrix path(s). Use space-separated list if multiple.",
    nargs='+'
)
parser.add_argument(
    '-d', '--num_rep_dna', 
    type=int, default=3, 
    help="Number of DNA replicates (default: 3)."
)
parser.add_argument(
    '-r', '--num_rep_rna', 
    type=int, default=3, 
    help="Number of RNA replicates (default: 3)."
)
parser.add_argument(
    '-o', "--out_path", 
    help="Output directory or file path for results."
)
parser.add_argument(
    '-t', '--tmp_path', 
    help="Temporary directory path to use during processing."
)
parser.add_argument(
    '--region', 
    help="Specific regions to process (space-separated list).",
    nargs='+'
)
parser.add_argument(
    '--prefiltered', 
    action='store_true', 
    help="Flag indicating that prefiltering has been performed."
)
parser.add_argument(
    '--filtered_dir',
    help="Directory containing prefiltered files. Use this when prefiltering has been done already."
)
parser.add_argument(
    '--region_only', 
    action='store_true', 
    help="Flag to check region coverage only, without further processing."
)

# Parse the command-line arguments.
args = parser.parse_args()

# ------------------------------
# Filter Count Matrix Based on Different Thresholds
# ------------------------------
# If the data has not been prefiltered, perform filtering.
if not args.prefiltered:
    # Get the total number of rows in the first count matrix file.
    num_rows = get_row_number(args.count_mat[0])
    
    # Construct the base directory for filtered data output.
    filtered_data_dir = os.path.join(args.out_path, 'filtered_data')
    set_dir(filtered_data_dir)
    
    # Extract the input filename (without extension) for naming output files.
    input_base = args.count_mat[0].split('/')[-1].split('.')[0]
    
    # Construct the path for the sorted BED file.
    sorted_bed_path = os.path.join(filtered_data_dir, f"{input_base}_sorted.bed")
    
    # Sort the BED file using safe_bedsort and save it to the filtered_data directory.
    safe_bedsort(args.count_mat[0], sorted_bed_path)
    
    # Open the sorted file for reading.
    sorted_file = open(sorted_bed_path, 'r')
    
    # Open output files for filtering based on the sum of DNA counts.
    out_sum_files = {
        5: open(os.path.join(filtered_data_dir, '5_sum.bed'), 'w'),
        10: open(os.path.join(filtered_data_dir, '10_sum.bed'), 'w'),
        20: open(os.path.join(filtered_data_dir, '20_sum.bed'), 'w'),
        50: open(os.path.join(filtered_data_dir, '50_sum.bed'), 'w'),
        100: open(os.path.join(filtered_data_dir, '100_sum.bed'), 'w')
    }
    
    # Open output files for filtering based on the number of samples passing a raw count threshold.
    out_raw_files = {
        1: open(os.path.join(filtered_data_dir, '1.bed'), 'w'),
        5: open(os.path.join(filtered_data_dir, '5.bed'), 'w'),
        10: open(os.path.join(filtered_data_dir, '10.bed'), 'w'),
        20: open(os.path.join(filtered_data_dir, '20.bed'), 'w'),
        50: open(os.path.join(filtered_data_dir, '50.bed'), 'w'),
        100: open(os.path.join(filtered_data_dir, '100.bed'), 'w')
    }
    
    # Set the minimum number of DNA replicates that must pass the raw count threshold.
    # Here, we assume that the DNA replicate columns are the first part of the replicate columns.
    min_dna_samples = args.num_rep_dna

    # Process each row in the sorted file.
    for row_idx in range(num_rows):
        # Read the current line and split into columns.
        line_fields = sorted_file.readline().strip().split('\t')
        
        # Extract DNA replicate counts from the last (num_rep_dna + num_rep_rna) columns.
        # Assume that the first args.num_rep_dna columns among these are the DNA counts.
        dna_counts = line_fields[-(args.num_rep_dna + args.num_rep_rna): -args.num_rep_rna]
        
        # ------------------------------
        # Filtering Based on Sum of DNA Counts
        # ------------------------------
        # Convert DNA counts to integers and compute their total.
        total_dna = np.sum([int(float(x)) for x in dna_counts])
        # For each sum threshold, write the line if the total DNA count meets or exceeds the threshold.
        for threshold, out_file in out_sum_files.items():
            if total_dna >= threshold:
                out_file.write('\t'.join(line_fields) + '\n')
        
        # ------------------------------
        # Filtering Based on the Number of Samples with a Minimum Raw Count
        # ------------------------------
        # For each threshold, check if at least min_dna_samples have counts above or equal to the threshold.
        for threshold, out_file in out_raw_files.items():
            passing_samples = [x for x in dna_counts if int(float(x)) >= threshold]
            if len(passing_samples) >= min_dna_samples:
                out_file.write('\t'.join(line_fields) + '\n')
    
    # Close all file handles.
    sorted_file.close()
    for f in out_sum_files.values():
        f.close()
    for f in out_raw_files.values():
        f.close()
    
    # Build the list of filtered file paths for downstream processing.
    file_list = [
        sorted_bed_path,
        os.path.join(filtered_data_dir, '5_sum.bed'),
        os.path.join(filtered_data_dir, '10_sum.bed'),
        os.path.join(filtered_data_dir, '20_sum.bed'),
        os.path.join(filtered_data_dir, '50_sum.bed'),
        os.path.join(filtered_data_dir, '100_sum.bed'),
        os.path.join(filtered_data_dir, '1.bed'),
        os.path.join(filtered_data_dir, '5.bed'),
        os.path.join(filtered_data_dir, '10.bed'),
        os.path.join(filtered_data_dir, '20.bed'),
        os.path.join(filtered_data_dir, '50.bed'),
        os.path.join(filtered_data_dir, '100.bed')
    ]
else:
    # If files are already prefiltered, use the provided filtered directory.
    input_base = args.count_mat[0].split('/')[-1].split('.')[0]
    file_list = [
        os.path.join(args.filtered_dir, f"{input_base}_sorted.bed"),
        os.path.join(args.filtered_dir, '5_sum.bed'),
        os.path.join(args.filtered_dir, '10_sum.bed'),
        os.path.join(args.filtered_dir, '20_sum.bed'),
        os.path.join(args.filtered_dir, '50_sum.bed'),
        os.path.join(args.filtered_dir, '100_sum.bed'),
        os.path.join(args.filtered_dir, '1.bed'),
        os.path.join(args.filtered_dir, '5.bed'),
        os.path.join(args.filtered_dir, '10.bed'),
        os.path.join(args.filtered_dir, '20.bed'),
        os.path.join(args.filtered_dir, '50.bed'),
        os.path.join(args.filtered_dir, '100.bed')
    ]
    
# Print the list of filtered file paths for verification.
print(file_list)


# Instantiate a CheckCoverage object to calculate genome-wide coverage.
coverage = CheckCoverage()

# Define the effective genome size for hg38.
hg38_effective_genome_size = 2913022398

# Only perform genome-wide coverage analysis if the 'region_only' flag is not set.
if not args.region_only:
    # Initialize a list to store coverage statistics for each file.
    coverage_stats = []
    
    # Loop through each file in the list of filtered files.
    for file_path in file_list:
        # Calculate the genome-wide coverage for the current file.
        # The method calc_genome_wide_coverage returns a tuple, and we use the second element (index 1)
        # which represents the number of base pairs covered.
        base_pairs_covered = coverage.calc_genome_wide_coverage(file_path)[1]
        
        # Get the total number of fragments (rows) in the file.
        num_fragments = get_row_number(file_path)
        
        # Extract the file name from the file path.
        file_name = file_path.split('/')[-1]
        
        # Append the collected statistics as a tuple:
        # (base pairs covered, number of fragments, file name, effective genome size)
        coverage_stats.append((base_pairs_covered, num_fragments, file_name, hg38_effective_genome_size))
    
    # Create a DataFrame from the list of coverage statistics.
    genome_coverage_df = pd.DataFrame(coverage_stats, columns=['num_bp', 'num_frag', 'query_file', 'size'])
    
    # Construct the output file path and save the DataFrame as a tab-separated text file.
    output_coverage_path = os.path.join(args.out_path, 'genome_wide_coverage_df.txt')
    genome_coverage_df.to_csv(output_coverage_path, sep='\t', index=False, header=True)
    
    print('Finished genome-wide coverage calculation.')
    
# ------------------------------
# Calculate Region-Specific Coverage
# ------------------------------
# Check if the user has provided any regions to analyze.
if len(args.region) > 0:
    # Initialize a list to collect coverage DataFrames for each specified region.
    region_coverage_list = []
    
    # Iterate over each region file provided in the command-line arguments.
    for region_file in args.region:
        # Print a message indicating the start of processing for the current region.
        print('Start processing region:', region_file)
        
        # Calculate the coverage for the current region using the pre-defined coverage object.
        # This function takes a region file and a list of filtered files as inputs.
        current_region_coverage = coverage.calc_region_coverage(region_file, file_list)
        
        # Append the resulting DataFrame to the list.
        region_coverage_list.append(current_region_coverage)
        
        # Print a message indicating that processing for the current region is complete.
        print('Finished processing region:', region_file)
    
    # Concatenate all region-specific coverage DataFrames into one master DataFrame.
    region_coverage_df = pd.concat(region_coverage_list, ignore_index=True, axis=0)
    
    # Construct the output file path for the region coverage DataFrame.
    output_region_coverage_path = os.path.join(args.out_path, 'region_coverage_df.txt')
    
    # Save the concatenated DataFrame to a tab-separated text file.
    region_coverage_df.to_csv(output_region_coverage_path, sep='\t', index=False, header=True)
    
    # Print a final message indicating that region coverage calculation is complete.
    print('Finished region coverage calculation.')

