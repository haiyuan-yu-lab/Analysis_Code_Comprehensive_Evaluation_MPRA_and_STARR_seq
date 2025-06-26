import numpy as np
import pandas as pd
import pybedtools
import os

# -------------------- PyBedTools Temp Dir Setup -------------------
# Set temporary directory for PyBedTools operations
pybedtools.helpers.set_tempdir('/fs/cbsuhy02/storage/jz855/tmp/')

# -------------------- Project Root -------------------
# Get the current working directory as the root for relative paths
project_root = os.getcwd()

# ==============================================================================
# Load GRO-cap Element Files (Bidirectional, Divergent, Unidirectional)
# ==============================================================================

# Load bidirectional elements (60bp resolution)
data_path = os.path.join(project_root, 'data', 'Reference', 'K562_GRO_cap', 'GRO_cap_element',
                         'original_files', 'GROcap-K562-02-Calls_1_bidirectional_peaks_element_60bp.bed')
bidirectional = pd.read_csv(data_path, sep='\t', header=None)
print(len(bidirectional))

# Load divergent elements
data_path = os.path.join(project_root, 'data', 'Reference', 'K562_GRO_cap', 'GRO_cap_element',
                         'original_files', 'GROcap-K562-02-Calls_1_divergent_peaks_element_60bp.bed')
divergent = pd.read_csv(data_path, sep='\t', header=None)
print(len(divergent))

# Load unidirectional elements
data_path = os.path.join(project_root, 'data', 'Reference', 'K562_GRO_cap', 'GRO_cap_element',
                         'original_files', 'GROcap-K562-02-Calls_1_unidirectional_peaks_element_60bp.bed')
unidirectional = pd.read_csv(data_path, sep='\t', header=None)
print(len(unidirectional))

# ==============================================================================
# Merge Bidirectional and Divergent GRO-cap Elements
# ==============================================================================

# Convert to index for set-based comparison
tmp1 = bidirectional.set_index([0,1,2])
tmp2 = divergent.set_index([0,1,2])

# Identify unique and shared entries between bidirectional and divergent sets
bidirectional_unique = set(tmp1.index.tolist()).difference(tmp2.index.tolist())
divergent_unique = set(tmp2.index.tolist()).difference(tmp1.index.tolist())
bidirectional_divergent_both = set(tmp1.index.tolist()).intersection(tmp2.index.tolist())

print(len(bidirectional_unique))
print(len(divergent_unique))
print(len(bidirectional_divergent_both))

# Combine both into a single DataFrame and annotate as 'divergent' by default
gro_cap_element_all = pd.concat([bidirectional[[0,1,2]], divergent[[0,1,2]]], ignore_index=True)
gro_cap_element_all['class'] = 'divergent'

# Overwrite to 'bidirectional' for entries uniquely in bidirectional set
gro_cap_element_all = gro_cap_element_all.set_index([0,1,2])
gro_cap_element_all.loc[list(bidirectional_unique), 'class'] = 'bidirectional'

# Reset index and remove duplicates
gro_cap_element_all = gro_cap_element_all.reset_index().drop_duplicates([0,1,2])
print(len(gro_cap_element_all))

# ==============================================================================
# Append Unidirectional Elements
# ==============================================================================

tmp = unidirectional[[0,1,2]]
tmp['class'] = 'unidirectional'

gro_cap_element_all = pd.concat([gro_cap_element_all, tmp], ignore_index=True)
print(len(gro_cap_element_all))
print(gro_cap_element_all.groupby('class').size())  # View breakdown by class

# Add size column (end - start)
gro_cap_element_all['size'] = gro_cap_element_all[2] - gro_cap_element_all[1]

# Save merged GRO-cap elements
data_path = os.path.join(project_root, 'data', 'Reference', 'K562_GRO_cap', 'GRO_cap_element',
                         'processed_files', 'GROcap_elements_all.bed.gz')
gro_cap_element_all.to_csv(data_path, sep='\t', index=False, header=False, compression='gzip')

# ==============================================================================
# Process ENCODE cCRE v4 BED File for K562
# ==============================================================================

# Load gzipped cCRE BED
data_path = os.path.join(project_root, 'data', 'Reference', 'cCRE_v4', 'ENCFF286VQG.bed.gz')
ccre_df = pd.read_csv(data_path, sep='\t', header=None)
print(len(ccre_df))

# Compute element size
ccre_df['size'] = ccre_df[2] - ccre_df[1]

# Keep [chr, start, end, ID, class (col 9), size]
ccre_df = ccre_df[[0, 1, 2, 3, 9, 'size']]

# Save processed cCRE BED
data_path = os.path.join(project_root, 'data', 'Reference', 'cCRE_v4', 'processed_files', 'K562_cCRE_v4_processed.bed.gz')
ccre_df.to_csv(data_path, sep='\t', index=False, header=False, compression='gzip')

# ==============================================================================
# Load and Process GENCODE v45 for Promoter Annotation
# ==============================================================================
import pyranges as pr

# Load GTF with PyRanges
data_path = os.path.join(project_root, 'data', 'Reference', 'hg38', 'GENCODE_v45', 'original_files', 'gencode.v45.annotation.gtf.gz')
v45 = pr.read_gtf(data_path)
print(len(v45))
# Convert to pandas DataFrame
v45 = v45.df

print(v45.groupby(['Feature']).size())

# Filter for transcripts
v45_protein_coding = v45[v45['Feature'] == 'transcript']
print(len(v45_protein_coding))

# Keep only protein-coding genes
v45_protein_coding = v45_protein_coding[v45_protein_coding['gene_type'] == 'protein_coding']
print(len(v45_protein_coding))

# Select relevant columns
v45_protein_coding = v45_protein_coding[['Chromosome', 'Start', 'End', 'Score', 'Strand',
                                         'Source', 'Feature', 'gene_id', 'gene_type', 'gene_name']]
print(len(v45_protein_coding))

# Remove duplicates based on location
v45_protein_coding = v45_protein_coding.drop_duplicates(['Chromosome', 'Start', 'End'])
print(len(v45_protein_coding))

# Sort and export full protein-coding transcript BED
v45_protein_coding = v45_protein_coding.sort_values(['Chromosome', 'Start'])
data_path = os.path.join(project_root, 'data', 'Reference', 'hg38', 'GENCODE_v45', 'processed_files',
                         'gencode_v45_protein_coding_annotation.bed.gz')
v45_protein_coding.to_csv(data_path, sep='\t', header=None, index=None, compression='gzip')
print(len(v45_protein_coding))

# ==============================================================================
# Define Promoter Regions (±500bp around TSS)
# ==============================================================================

# Separate by strand
v45_protein_coding_plus = v45_protein_coding[v45_protein_coding['Strand'] == '+']
print(len(v45_protein_coding_plus))
v45_protein_coding_minus = v45_protein_coding[v45_protein_coding['Strand'] == '-']
print(len(v45_protein_coding_minus))

# Compute proximal regions
distance_list = [200, 500, 2000]
filename_str_list = ['400b', '1kb', '4kb']

for d, f in zip(distance_list, filename_str_list):
    print(d)
    print(f)
    promoter_plus = v45_protein_coding_plus.copy()
    promoter_plus['Start_p'] = promoter_plus['Start'] - d
    promoter_plus['End_p'] = promoter_plus['Start'] + d

    promoter_minus = v45_protein_coding_minus.copy()
    promoter_minus['Start_p'] = promoter_minus['End'] - d
    promoter_minus['End_p'] = promoter_minus['End'] + d

    # Merge promoters from both strands
    promoter_both = pd.concat([promoter_plus, promoter_minus], ignore_index=True)
    print(len(promoter_both))

    # Rearrange and rename columns for BED-like format
    promoter_both.columns = [0,1,2,3,4,5,6,7,8,9,10,11]
    promoter_both = promoter_both[[0,10,11,3,4,1,2,5,6,7,8,9]]
    promoter_both = promoter_both.sort_values([0, 1])
    promoter_both = promoter_both.drop_duplicates()
    print(len(promoter_both))

    # Remove duplicated TSS-centered regions
    promoter_both = promoter_both.drop_duplicates([0, 10, 11])
    print(len(promoter_both))

    # Save to BED file
    fname = 'promoters_'+f+'_GENCODE_v45_protein_coding_tss_centered.bed.gz'
    data_path = os.path.join(project_root, 'data', 'Reference', 'hg38', 'GENCODE_v45', 'processed_files',fname)
    promoter_both.to_csv(data_path, sep='\t', header=False, index=False)
    print('-------')