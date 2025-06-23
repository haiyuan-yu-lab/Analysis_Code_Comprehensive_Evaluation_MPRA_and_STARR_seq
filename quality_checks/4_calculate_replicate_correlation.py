#!/usr/bin/env python3
"""
replicate_correlation_analysis.py

Description:
  This script performs replicate correlation analysis on count data.
  It calculates pairwise correlations for replicate counts using:
    - Raw Counts Per Million (CPM)
    - Log-transformed CPM
    - Log2 ratio (derived from the log-transformed CPM)

Inputs:
  --count_mat : Path to the input count matrix file (tab-separated, without header).
                The count matrix is expected to have many columns, where the last columns correspond
                to the replicates for DNA and RNA.
  --num_DNA_rep : Number of DNA replicates present in the count matrix (default: 3).
  --num_RNA_rep : Number of RNA replicates present in the count matrix (default: 3).
  --dname       : Dataset name, which will be appended to the output for clarity.
  --out_path    : Path to the output directory where result files will be saved. Ensure the path ends with a '/'.

Outputs:
  The script generates three output files in the specified output directory:
    1. cpm_count_corr.txt     - Pairwise correlation results using raw CPM values.
    2. logcpm_count_corr.txt  - Pairwise correlation results using log-transformed CPM values.
    3. log2ratio_corr.txt     - Pairwise correlation results using log2 ratio values.
  
Usage:
  python replicate_correlation_analysis.py --count_mat path/to/count_matrix.txt \
    --out_path path/to/output_dir/ --num_DNA_rep 3 --num_RNA_rep 3 --dname dataset_name
"""

import argparse
import os
import pandas as pd
import gzip
import sys

# Dynamically add the 'src' directory to sys.path so that modules from the project can be imported.
project_root = os.path.abspath(os.path.join(os.getcwd(), '..'))
src_path = os.path.join(project_root, 'src')
sys.path.insert(0, src_path)

from quality_control import CheckRepCorr

# -------------------------
# Parse command-line arguments
# -------------------------
parser = argparse.ArgumentParser(
    description="Perform replicate correlation analysis on count data using CPM, logCPM, and log2 ratios."
)
parser.add_argument('-c', '--count_mat',
                    help='Path to input count matrix (tab-separated file without header)')
parser.add_argument('-o', '--out_path',
                    help='Path to the output directory where results will be saved (e.g., /path/to/output/)')
parser.add_argument('-d', '--num_DNA_rep', type=int, default=3,
                    help='Number of DNA replicates in the count matrix (default: 3)')
parser.add_argument('-r', '--num_RNA_rep', type=int, default=3,
                    help='Number of RNA replicates in the count matrix (default: 3)')
parser.add_argument('-n', '--dname',
                    help='Dataset name to be added to output files for reference')
args = parser.parse_args()

# -------------------------
# Prepare column names based on replication counts
# -------------------------
# Create descriptive column names for the DNA and RNA replicate columns.
# It is assumed that the last (num_DNA_rep + num_RNA_rep) columns in the count matrix correspond to these replicates.
dna_columns = [f'DNA{i}' for i in range(1, args.num_DNA_rep + 1)]
rna_columns = [f'RNA{i}' for i in range(1, args.num_RNA_rep + 1)]

# -------------------------
# Read and concatenate the count matrix
# -------------------------
# Read the input count matrix in manageable chunks to handle large files without exhausting memory.
chunk_size = 1000000  # Number of rows per chunk
data_chunks = []     # List to collect each chunk

# Iterate through the count matrix file chunk by chunk.
for chunk in pd.read_csv(args.count_mat, sep='\t', header=None, chunksize=chunk_size):
    # Extract only the last columns corresponding to the replicates.
    count_columns = chunk.columns.tolist()[-(args.num_DNA_rep + args.num_RNA_rep):]
    chunk = chunk[count_columns]
    # Rename these columns with the descriptive DNA and RNA replicate names.
    chunk.columns = dna_columns + rna_columns
    data_chunks.append(chunk)

# Combine all chunks into a single DataFrame.
data = pd.concat(data_chunks, axis=0, ignore_index=True)
print("Total rows in data:", len(data))

# -------------------------
# Run replicate correlation analysis
# -------------------------
# Initialize the CheckRepCorr object, which handles the correlation calculations.
qc = CheckRepCorr(num_rep_dna=args.num_DNA_rep, num_rep_rna=args.num_RNA_rep)

# (1) Calculate raw Counts Per Million (CPM) values (without log transformation).
raw_cpm = qc.calc_cpm(data, log=False)

# (2) Compute pairwise correlations for DNA replicates using raw CPM values.
dna_cpm_corr = qc.compute_pairwise_corr(raw_cpm[dna_columns], dna_columns)
dna_cpm_corr['count_type'] = 'cpm'  # Label the result with the count type

# (3) Compute pairwise correlations for RNA replicates using raw CPM values.
rna_cpm_corr = qc.compute_pairwise_corr(raw_cpm[rna_columns], rna_columns)
rna_cpm_corr['count_type'] = 'cpm'

# Combine DNA and RNA correlation results into one DataFrame.
cpm_corr_df = pd.concat([dna_cpm_corr, rna_cpm_corr], ignore_index=True)
print("Total pairwise correlations (CPM):", len(cpm_corr_df))

# Append the dataset name to the results and save to output file.
cpm_corr_df['Dataset'] = args.dname
cpm_corr_df.to_csv(args.out_path + 'cpm_count_corr.txt', sep='\t', header=True, index=False)
print('Finished CPM count correlation')

# -------------------------
# Calculate and analyze log-transformed CPM
# -------------------------
# Compute log-transformed CPM values.
log_cpm = qc.calc_cpm(data, log=True)

# Compute pairwise correlations for DNA replicates using log-transformed CPM.
dna_logcpm_corr = qc.compute_pairwise_corr(log_cpm[dna_columns], dna_columns)
dna_logcpm_corr['count_type'] = 'logcpm'

# Compute pairwise correlations for RNA replicates using log-transformed CPM.
rna_logcpm_corr = qc.compute_pairwise_corr(log_cpm[rna_columns], rna_columns)
rna_logcpm_corr['count_type'] = 'logcpm'

# Combine DNA and RNA logCPM correlation results.
logcpm_corr_df = pd.concat([dna_logcpm_corr, rna_logcpm_corr], ignore_index=True)
print("Total pairwise correlations (logCPM):", len(logcpm_corr_df))

# Append the dataset name to the results and save to output file.
logcpm_corr_df['Dataset'] = args.dname
logcpm_corr_df.to_csv(args.out_path + 'logcpm_count_corr.txt', sep='\t', header=True, index=False)
print('Finished logCPM count correlation')

# -------------------------
# Calculate and analyze log2 ratio
# -------------------------
# Compute the log2 ratio values using the log-transformed CPM data.
log2_ratio = qc.calc_log2Ratio(log_cpm)
replicates_list = log2_ratio.columns.tolist()

# Compute pairwise correlations for the log2 ratio values.
log2_ratio_corr = qc.compute_pairwise_corr(log2_ratio, replicates_list)

# Save the log2 ratio correlation results to an output file.
log2_ratio_corr.to_csv(args.out_path + 'log2ratio_corr.txt', sep='\t', header=True, index=False)
print('Finished log2 ratio correlation')