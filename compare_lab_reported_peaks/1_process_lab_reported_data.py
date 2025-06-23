import numpy as np
import pandas as pd
import pybedtools
import os


# -------------------- PyBedTools Temp Dir Setup -------------------
pybedtools.helpers.set_tempdir('/fs/cbsuhy02/storage/jz855/tmp/')

def merge_bins(data):
    
    data = data[['seqnames', 'start', 'end']]
    data = data.sort_values(['seqnames', 'start'])
    
    bed = pybedtools.BedTool.from_dataframe(data)
    merged = bed.merge()
    merged = merged.to_dataframe(disable_auto_names=True, header=None)
    
    merged = merged.sort_values([0,1])
    
    return(merged)

##############################################################################

# Specify root directory
project_root = os.path.abspath(os.path.join(os.getcwd(), '..'))

# LentiMPRA

# Read in measured activity values for each element
# Build full path to the file
data_path = os.path.join(project_root, 'data', 'lab_reported_data', 'original_files', 'LentiMPRA', 'K562_data_Vikram_2024_0521.txt')
lenti_activity = pd.read_csv(data_path, sep='\t')


# Read in metadata for each element (e.g., location, name, strand, category)
data_path = os.path.join(project_root, 'data', 'lab_reported_data', 'original_files', 'LentiMPRA', 'lentiMPRA_vikarm_2024_0521.txt')
lenti_reference = pd.read_csv(data_path, sep='\t')
print(len(lenti_reference))

# Merge activity data with metadata on 'name'
lenti_reference = lenti_reference.set_index('name')
lenti_activity = lenti_activity.set_index('name')
lenti_data = lenti_reference.join(lenti_activity, how='outer')
print(len(lenti_data))

# Reset index for downstream analysis
lenti_data = lenti_data.reset_index()

# Subset to relevant columns
lenti_data = lenti_data[['chr', 'start', 'end', 'name', 'strand', 'category', 
                         'replicate 1 [log2(rna/dna)]', 
                         'replicate 2 [log2(rna/dna)]', 
                         'replicate 3 [log2(rna/dna)]', 'mean']]

# Replace any missing values with 'NA'
lenti_data = lenti_data.fillna('NA')

# Identify negative controls from merged dataset
neg_ctrl = lenti_data[lenti_data['name'].map(lambda x: 'shuffled' in x)]

# Drop rows without mean logFC value
neg_ctrl = neg_ctrl[neg_ctrl['mean'] != 'NA']

lenti_data = lenti_data[lenti_data['mean'] != 'NA']

# Compute 95th percentile of negative control activity as a threshold for calling 'active'
neg_ctrl_activity = neg_ctrl['mean'].values
neg_ctrl_95pct = np.percentile(neg_ctrl_activity, 95)

# Label each element as 'active' or 'inactive' based on threshold
lenti_data['class'] = 'inactive'
lenti_data.loc[lenti_data[lenti_data['mean'] >= neg_ctrl_95pct].index.tolist(), 'class'] = 'active'
print(lenti_data.groupby(['class']).size())

# ------------------ Subset Genomic Elements ------------------

# Keep only elements with valid start coordinate (i.e., genomic elements)
lenti_genomic_element = lenti_data[(lenti_data['start'] != 'NA')]
lenti_genomic_element = lenti_genomic_element.astype({'start':float, 'end':float})
lenti_genomic_element = lenti_genomic_element.astype({'start':int, 'end':int})

# Show breakdown by category and active/inactive class
print(lenti_genomic_element.groupby(['category', 'class']).size())

# ------------------ Merge Inactive Elements ------------------

lenti_inactive = lenti_genomic_element[(lenti_genomic_element['class'] == 'inactive')]
print(len(lenti_inactive))

# Format BED-style columns
lenti_inactive = lenti_inactive[['chr', 'start', 'end', 'strand']]
lenti_inactive.columns = ['seqnames', 'start', 'end', 'strand']

# Split by strand and merge overlapping elements
forward = lenti_inactive[lenti_inactive['strand'] == '+']
reverse = lenti_inactive[lenti_inactive['strand'] == '-']

lenti_inactive_forward_merged = merge_bins(forward)
lenti_inactive_forward_merged[3] = '+'

lenti_inactive_reverse_merged = merge_bins(reverse)
lenti_inactive_reverse_merged[3] = '-'

# Concatenate merged elements from both strands
lenti_inactive_merged = pd.concat([lenti_inactive_forward_merged, lenti_inactive_reverse_merged], ignore_index=True, axis=0)
print(len(lenti_inactive_merged))

# Calculate size of each merged region
lenti_inactive_merged[4] = lenti_inactive_merged[2]-lenti_inactive_merged[1]

# ------------------ Merge Active Elements ------------------

lenti_active = lenti_genomic_element[(lenti_genomic_element['class'] == 'active')]
print(len(lenti_active))

# Format BED-style columns
lenti_active = lenti_active[['chr', 'start', 'end', 'strand']]
lenti_active.columns = ['seqnames', 'start', 'end', 'strand']

# Split by strand and merge overlapping elements
forward = lenti_active[lenti_active['strand'] == '+']
reverse = lenti_active[lenti_active['strand'] == '-']

lenti_active_forward_merged = merge_bins(forward)
lenti_active_forward_merged[3] = '+'

lenti_active_reverse_merged = merge_bins(reverse)
lenti_active_reverse_merged[3] = '-'

# Combine merged active elements
lenti_active_merged = pd.concat([lenti_active_forward_merged, lenti_active_reverse_merged], ignore_index=True, axis=0)
print(len(lenti_active_merged))

# Calculate region size
lenti_active_merged[4] = lenti_active_merged[2]-lenti_active_merged[1]

# ------------------ Deduplicate Active Elements ------------------

lenti_active_merged_unique = lenti_active_merged.drop_duplicates([0,1,2])
print(len(lenti_active_merged_unique))

# ------------------ Save Outputs ------------------

data_path = os.path.join(project_root, 'data', 'lab_reported_data', 'processed_files', 'LentiMPRA', 'inactive_merged_all.bed')
lenti_inactive_merged.to_csv(data_path, sep='\t', header=None, index=None)

data_path = os.path.join(project_root, 'data', 'lab_reported_data', 'processed_files', 'LentiMPRA', 'active_merged_all.bed')
lenti_active_merged.to_csv(data_path, sep='\t', header=None, index=None)

data_path = os.path.join(project_root, 'data', 'lab_reported_data', 'processed_files', 'LentiMPRA', 'active_merged_all_unique_either_orientation.bed')
lenti_active_merged_unique.to_csv(data_path, sep='\t', header=False, index=False)
print(len(lenti_active_merged_unique))
print('finished processing LentiMPRA')

##############################################################################

# Read original tiling MPRA calls (hg19) from ENCODE file
data_path = os.path.join(project_root, 'data', 'lab_reported_data', 'original_files', 'TilingMPRA', 'OL13_ENCFF436CRS.bed.gz')
ol13_original = pd.read_csv(data_path, sep='\t', compression='gzip', header=None)
print(len(ol13_original))

# Convert 1-based start to 0-based for BED format
ol13_original[1] = ol13_original[1] - 1

# Convert -log10(p) and -log10(adj.p) to raw p-values
ol13_original['raw_p'] = 10 ** (-ol13_original[9].values)
ol13_original['adj_p'] = 10 ** (-ol13_original[10].values)

# Read hg38-lifted version and merge info
data_path = os.path.join(project_root, 'data', 'lab_reported_data', 'original_files', 'TilingMPRA', 'OL13_hg38.bed')
ol13_hg38 = pd.read_csv(data_path, sep='\t', header=None)
print(len(ol13_hg38))

# Reformat hg38 columns for joining
ol13_hg38 = ol13_hg38.set_index(3)
ol13_hg38 = ol13_hg38[[0, 1, 2]]
ol13_hg38.columns = ['seqnames', 'start', 'end']

# Join original data (hg19) with lifted coordinates (hg38)
ol13_original = ol13_original.set_index(3)
ol13_original = ol13_original.join(ol13_hg38, how='outer')
ol13_original = ol13_original.reset_index()

# -----------------------------------------------------
# Call enhancer/inactive using original criteria:
# log2FC ≥ 1 & adj.p < 0.01 → enhancer
# -----------------------------------------------------
ol13_original['original_call'] = 'inactive'
ol13_original.loc[ol13_original[(ol13_original[6] >= 1) & 
                                (ol13_original['adj_p'] < 0.01)].index.tolist(), 'original_call'] = 'enhancer'

# --------------------------------------------
# Split into enhancer vs inactive sets
# --------------------------------------------
ol13_original_enhancer = ol13_original[ol13_original['original_call'] == 'enhancer']
print(len(ol13_original_enhancer))

ol13_original_inactive = ol13_original[ol13_original['original_call'] == 'inactive']
print(len(ol13_original_inactive))

# --------------------------------------------
# Merge inactive regions across genome
# --------------------------------------------
ol13_original_inactive = ol13_original_inactive[[0, 1, 2, 5]]
ol13_original_inactive.columns = ['seqnames', 'start', 'end', 'strand']

# Merge overlapping or adjacent bins
ol13_original_inactive_merged = merge_bins(ol13_original_inactive)
print(len(ol13_original_inactive_merged))

# Add region size
ol13_original_inactive_merged[4] = ol13_original_inactive_merged[2] - ol13_original_inactive_merged[1]
print(ol13_original_inactive_merged[4].describe())

# Save merged inactive regions
data_path = os.path.join(project_root, 'data', 'lab_reported_data', 'processed_files', 'TilingMPRA', 'OL13_merged_inactive_either_orientation.bed')
ol13_original_inactive_merged.to_csv(data_path, sep='\t', header=None, index=None)

# --------------------------------------------
# Reformat enhancer calls to match expected schema
# --------------------------------------------
ol13_original_enhancer = ol13_original_enhancer[[3, 'seqnames', 'start', 'end', 6, 'raw_p', 'adj_p', 'original_call']]
ol13_original_enhancer = ol13_original_enhancer.astype({'start': int, 'end': int})
ol13_original_enhancer.columns = ['ID', 'chr', 'start', 'stop', 'log2FoldChange', 'pvalue', 'padj', 'original_call']

# --------------------------------------------
# Prepare enhancer elements for peak merging
# --------------------------------------------
enhancer = ol13_original_enhancer[ol13_original_enhancer['original_call'] == 'enhancer']
enhancer['strand'] = '+'
enhancer['z_score'] = 0
enhancer['neg_ctrl'] = '.'

# Rename columns to match internal format
enhancer = enhancer[['ID', 'chr', 'start', 'stop', 'log2FoldChange', 'strand', 'pvalue', 'padj', 'z_score', 'neg_ctrl', 'original_call']]
enhancer.columns = ['ID', 'seqnames', 'start', 'end', 'logFC', 'strand', 'P.Value', 'adj.P.Val', 'z_score', 'neg_ctrl', 'call']

merged = merge_bins(enhancer)
merged['strand'] = '.'
merged['size'] = merged[2].values - merged[1].values

data_path = os.path.join(project_root, 'data', 'lab_reported_data', 'processed_files', 'TilingMPRA', 'OL13_merged_enhancer_peaks_in_either_orientation.bed')
merged.to_csv(data_path, sep='\t', index=False, header=False)
print('finished saving TilingMPRA OL13')

##############################################################################

# ------------------------------------------------------
# Load associated BED results (with -log10(p) scores)
# ------------------------------------------------------
data_path = os.path.join(project_root, 'data', 'lab_reported_data', 'original_files', 'TilingMPRA', 'OL45_ENCFF558RYR.bed.gz')
ol45_bed = pd.read_csv(data_path, sep='\t', compression='gzip', header=None)
print(len(ol45_bed))

# ------------------------------------------------------
# Convert -log10(p) and -log10(adj.p) to raw values
# ------------------------------------------------------
ol45_bed['raw_p'] = 10 ** (-ol45_bed[9].values)
ol45_bed['adj_p'] = 10 ** (-ol45_bed[10].values)

# ------------------------------------------------------
# Classify fragments as enhancers if logFC ≥ 1 and adj.p < 0.01
# ------------------------------------------------------
ol45_bed['original_call'] = 'inactive'
ol45_bed.loc[ol45_bed[(ol45_bed[6] >= 1) & 
                      (ol45_bed['adj_p'] < 0.01)].index.tolist(), 'original_call'] = 'enhancer'

print(ol45_bed.groupby(['original_call']).size())

# ------------------------------------------------------
# Separate enhancer vs inactive fragments
# ------------------------------------------------------
ol45_original_enhancer = ol45_bed[ol45_bed['original_call'] == 'enhancer']
print(len(ol45_original_enhancer))

ol45_original_inactive = ol45_bed[ol45_bed['original_call'] == 'inactive']
print(len(ol45_original_inactive))

# ------------------------------------------------------
# Merge inactive fragments across genome
# ------------------------------------------------------
ol45_original_inactive = ol45_original_inactive[[0, 1, 2, 5]]  # [chrom, start, end, strand]
ol45_original_inactive.columns = ['seqnames', 'start', 'end', 'strand']

# Merge overlapping bins
ol45_original_inactive_merged = merge_bins(ol45_original_inactive)

# Add region size
ol45_original_inactive_merged[4] = ol45_original_inactive_merged[2] - ol45_original_inactive_merged[1]

# Save merged inactive regions
data_path = os.path.join(project_root, 'data', 'lab_reported_data', 'processed_files', 'TilingMPRA', 'OL45_merged_inactive_either_orientation.bed')
ol45_original_inactive_merged.to_csv(data_path, sep='\t', header=None, index=None)

# ------------------------------------------------------
# Format enhancer fragments for peak merging and export
# ------------------------------------------------------
ol45_original_enhancer = ol45_original_enhancer[[0, 1, 2, 3, 6, 'raw_p', 'adj_p', 'original_call']]
ol45_original_enhancer.columns = ['chr', 'start', 'end', 'ID', 'logFC', 'raw_p', 'adj_p', 'original_call']
ol45_original_enhancer.head(1)

# Reformat and annotate enhancer bins
enhancer = ol45_original_enhancer[ol45_original_enhancer['original_call'] == 'enhancer']
enhancer['strand'] = '+'
enhancer['z_score'] = 0
enhancer['neg_ctrl'] = '.'

# Rename columns to match internal format
enhancer = enhancer[['ID', 'chr', 'start', 'end', 'logFC', 'strand', 'raw_p', 'adj_p', 'z_score', 'neg_ctrl', 'original_call']]
enhancer.columns = ['ID', 'seqnames', 'start', 'end', 'logFC', 'strand', 'P.Value', 'adj.P.Val', 'z_score', 'neg_ctrl', 'call']

merged = merge_bins(enhancer)
merged['strand'] = '.'
merged['size'] = merged[2].values - merged[1].values

data_path = os.path.join(project_root, 'data', 'lab_reported_data', 'processed_files', 'TilingMPRA', 'OL45_merged_enhancer_peaks_in_either_orientation.bed')
merged.to_csv(data_path, sep='\t', index=False, header=False)
print('finished saving TilingMPRA OL45')

##############################################################################

# ------------------------------------------------------
# Load BED file with MPRA statistics (e.g., logFC, -log10(p), -log10(adj.p))
# ------------------------------------------------------
data_path = os.path.join(project_root, 'data', 'lab_reported_data', 'original_files', 'TilingMPRA', 'OL43_ENCFF879JQC.bed.gz')
ol43_bed = pd.read_csv(data_path, sep='\t', header=None, compression='gzip')
print(len(ol43_bed))

# ------------------------------------------------------
# Convert -log10(p-value) and -log10(adj.p) to raw p-values
# ------------------------------------------------------
ol43_bed['raw_p'] = 10 ** (-ol43_bed[9].values)
ol43_bed['adj_p'] = 10 ** (-ol43_bed[10].values)

# ------------------------------------------------------
# Assign enhancer label to bins based on logFC and adj.p thresholds
# ------------------------------------------------------
ol43_bed['original_call'] = 'inactive'
ol43_bed.loc[ol43_bed[(ol43_bed[6] >= 1) & 
                      (ol43_bed['adj_p'] < 0.01)].index.tolist(), 'original_call'] = 'enhancer'
print(ol43_bed.groupby(['original_call']).size())

# ------------------------------------------------------
# Split enhancer vs inactive bins
# ------------------------------------------------------
ol43_original_enhancer = ol43_bed[ol43_bed['original_call'] == 'enhancer']
print(len(ol43_original_enhancer))

ol43_original_inactive = ol43_bed[ol43_bed['original_call'] == 'inactive']
print(len(ol43_original_inactive))

# ------------------------------------------------------
# Prepare inactive bins for merging (format as BED)
# ------------------------------------------------------
ol43_original_inactive = ol43_original_inactive[[0, 1, 2, 5]]  # [chr, start, end, strand]
ol43_original_inactive.columns = ['seqnames', 'start', 'end', 'strand']

# ------------------------------------------------------
# Merge overlapping or adjacent inactive bins
# ------------------------------------------------------
ol43_original_inactive_merged = merge_bins(ol43_original_inactive)

# Compute size of each merged bin
ol43_original_inactive_merged[4] = ol43_original_inactive_merged[2] - ol43_original_inactive_merged[1]

# Save merged inactive bins to BED file
data_path = os.path.join(project_root, 'data', 'lab_reported_data', 'processed_files', 'TilingMPRA', 'OL43_merged_inactive_either_orientation.bed')
ol43_original_inactive_merged.to_csv(data_path,sep='\t', header=None, index=None)

# ------------------------------------------------------
# Prepare enhancer bins with relevant metadata for merging
# ------------------------------------------------------
ol43_original_enhancer = ol43_original_enhancer[[0, 1, 2, 3, 6, 'raw_p', 'adj_p', 'original_call']]
ol43_original_enhancer.columns = ['chr', 'start', 'end', 'ID', 'logFC', 'raw_p', 'adj_p', 'original_call']

# Reformat and annotate enhancer bins
enhancer = ol43_original_enhancer[ol43_original_enhancer['original_call'] == 'enhancer']
enhancer['strand'] = '+'
enhancer['z_score'] = 0
enhancer['neg_ctrl'] = '.'

# Reorder and rename columns to match unified format
enhancer = enhancer[['ID', 'chr', 'start', 'end', 'logFC', 'strand', 'raw_p', 'adj_p', 'z_score', 'neg_ctrl', 'original_call']]
enhancer.columns = ['ID', 'seqnames', 'start', 'end', 'logFC', 'strand', 'P.Value', 'adj.P.Val', 'z_score', 'neg_ctrl', 'call']

merged = merge_bins(enhancer)
merged['strand'] = '.'
merged['size'] = merged[2] - merged[1]

data_path = os.path.join(project_root, 'data', 'lab_reported_data', 'processed_files', 'TilingMPRA', 'OL43_merged_enhancer_peaks_in_either_orientation.bed')
merged.to_csv(data_path, sep='\t', index=False, header=False)
print('finished saving TilingMPRA OL43')

##############################################################################

# No need to process WHG-STARR-seq

##############################################################################

# ATAC-STARR-seq

# Load ATAC-STARR-seq output with csaw results (q-values, logFC, etc.)
data_path = os.path.join(project_root, 'data', 'lab_reported_data', 'original_files', 'ATAC_STARR_seq', 'atacSTARR.ultra_deep.corrected.csaw.hg38.v10.common_file_formatted.txt.gz')
atac_starr_csaw = pd.read_csv(data_path, sep='\t')
print(len(atac_starr_csaw))  # Total number of bins/fragments

# --------------------------------------------
# Filter for active bins:
# - Q-value threshold: -log10(q) ≥ 3  =>  q ≤ 0.001
# - logFC > 0          =>  enhancer-like activity
# --------------------------------------------
atac_starr_csaw['original_call'] = 'inactive'
atac_starr_csaw.loc[atac_starr_csaw[(atac_starr_csaw['logFC'] > 0 ) & (atac_starr_csaw['minusLog10QValue'] >= 3)].index.tolist(), 'original_call'] = 'enhancer'
print(atac_starr_csaw.groupby(['original_call']).size())

# ------------------------------------------------------
# Split enhancer vs inactive bins
# ------------------------------------------------------
atac_starr_original_enhancer = atac_starr_csaw[atac_starr_csaw['original_call'] == 'enhancer']
print(len(ol43_original_enhancer))

atac_starr_original_inactive = atac_starr_csaw[atac_starr_csaw['original_call'] == 'inactive']
print(len(atac_starr_original_inactive))

# ------------------------------------------------------
# Prepare inactive bins for merging (format as BED)
# ------------------------------------------------------
atac_starr_original_inactive = atac_starr_original_inactive[['seqnames', 'start' ,'end', 'name', 'score', 'strand']]

# ------------------------------------------------------
# Merge overlapping or adjacent inactive bins
# ------------------------------------------------------
atac_starr_original_inactive_merged = merge_bins(atac_starr_original_inactive)

# Compute size of each merged bin
atac_starr_original_inactive_merged[4] = atac_starr_original_inactive_merged[2] - atac_starr_original_inactive_merged[1]

# Save merged inactive bins to BED file
data_path = os.path.join(project_root, 'data', 'lab_reported_data', 'processed_files', 'ATAC_STARR_seq', 'ATAC_STARR_merged_inactive_either_orientation.bed')
atac_starr_original_inactive_merged.to_csv(data_path, sep='\t', header=None, index=None)

# Enhancers
# --------------------------------------------
# Format as BED-like DataFrame for downstream use
# --------------------------------------------
atac_starr_original_enhancer = atac_starr_original_enhancer[['seqnames', 'start' ,'end', 'name', 'score', 'strand']]

atac_starr_original_enhancer_merged = merge_bins(atac_starr_original_enhancer)
atac_starr_original_enhancer_merged[4] = atac_starr_original_enhancer_merged[2]-atac_starr_original_enhancer_merged[1]

# --------------------------------------------
# Save to BED file (e.g., for visualization or overlap analysis)
# --------------------------------------------
data_path = os.path.join(project_root, 'data', 'lab_reported_data', 'processed_files', 'ATAC_STARR_seq', 'ATAC_STARR_merged_enhancer_peaks_in_either_orientation.bed')
atac_starr_original_enhancer_merged.to_csv(data_path, sep='\t', header=False, index=False)
print('finished saving ATAC-STARR-seq')

##############################################################################

# Merge lab-reported into one single bed file

# -------------------- LentiMPRA --------------------
data_path = os.path.join(project_root, 'data', 'lab_reported_data', 'processed_files', 'LentiMPRA', 'active_merged_all_unique_either_orientation.bed')
lenti_merged_unique_active = pd.read_csv(data_path, sep='\t', header=None)
print(len(lenti_merged_unique_active))

lenti_merged_unique_active = lenti_merged_unique_active[[0,1,2]]
lenti_merged_unique_active[3] = 'LentiMPRA_ENCSR382BVV_active'

data_path = os.path.join(project_root, 'data', 'lab_reported_data', 'processed_files', 'LentiMPRA', 'inactive_merged_all.bed')
lenti_inactive_merged = pd.read_csv(data_path, sep='\t', header=None)
print(len(lenti_inactive_merged))

lenti_inactive_merged = lenti_inactive_merged.drop_duplicates([0,1,2])
print(len(lenti_inactive_merged))

lenti_inactive_merged = lenti_inactive_merged[[0,1,2]]
lenti_inactive_merged[3] = 'LentiMPRA_ENCSR382BVV_inactive'

lenti = pd.concat([lenti_merged_unique_active, lenti_inactive_merged], axis=0, ignore_index=True)
print(len(lenti))

print('LentiMPRA')
print(lenti.groupby([3]).size())

# -------------------- TilingMPRA OL13 --------------------

data_path = os.path.join(project_root, 'data', 'lab_reported_data', 'processed_files', 'TilingMPRA', 'OL13_merged_enhancer_peaks_in_either_orientation.bed')
ol13_merged_active = pd.read_csv(data_path, sep='\t', header=None)
print(len(ol13_merged_active))

ol13_merged_active = ol13_merged_active[[0,1,2]]
ol13_merged_active[3] = 'TilingMPRA_ENCSR394HXI_active'

data_path = os.path.join(project_root, 'data', 'lab_reported_data', 'processed_files', 'TilingMPRA', 'OL13_merged_inactive_either_orientation.bed')
ol13_merged_inactive = pd.read_csv(data_path, sep='\t', header=None)
print(len(ol13_merged_inactive))

ol13_merged_inactive = ol13_merged_inactive[[0,1,2]]
ol13_merged_inactive[3] = 'TilingMPRA_ENCSR394HXI_inactive'

ol13_merged_inactive.head(3)

tilingMPRA_ENCSR394HXI = pd.concat([ol13_merged_active, ol13_merged_inactive], axis=0, ignore_index=True)
print(len(tilingMPRA_ENCSR394HXI))

print('TilingMPRA OL13')
print(tilingMPRA_ENCSR394HXI.groupby([3]).size())

# -------------------- TilingMPRA OL45 --------------------

data_path = os.path.join(project_root, 'data', 'lab_reported_data', 'processed_files', 'TilingMPRA', 'OL45_merged_enhancer_peaks_in_either_orientation.bed')
ol45_merged_active = pd.read_csv(data_path, sep='\t', header=None)
print(len(ol45_merged_active))

ol45_merged_active = ol45_merged_active[[0,1,2]]
ol45_merged_active[3] = 'TilingMPRA_ENCSR363XER_active'

data_path = os.path.join(project_root, 'data', 'lab_reported_data', 'processed_files', 'TilingMPRA', 'OL45_merged_inactive_either_orientation.bed')
ol45_merged_inactive = pd.read_csv(data_path, sep='\t', header=None)
print(len(ol45_merged_inactive))

ol45_merged_inactive = ol45_merged_inactive[[0,1,2]]
ol45_merged_inactive[3] = 'TilingMPRA_ENCSR363XER_inactive'

tilingMPRA_ENCSR363XER = pd.concat([ol45_merged_active, ol45_merged_inactive], axis=0, ignore_index=True)
print(len(tilingMPRA_ENCSR363XER))

print('TilingMPRA OL45')
print(tilingMPRA_ENCSR363XER.groupby([3]).size())

# -------------------- TilingMPRA OL43 --------------------

data_path = os.path.join(project_root, 'data', 'lab_reported_data', 'processed_files', 'TilingMPRA', 'OL43_merged_enhancer_peaks_in_either_orientation.bed')
ol43_merged_active = pd.read_csv(data_path, sep='\t', header=None)
print(len(ol43_merged_active))

ol43_merged_active = ol43_merged_active[[0,1,2]]
ol43_merged_active[3] = 'TilingMPRA_ENCSR917SFD_active'

data_path = os.path.join(project_root, 'data', 'lab_reported_data', 'processed_files', 'TilingMPRA', 'OL43_merged_inactive_either_orientation.bed')
ol43_merged_inactive = pd.read_csv(data_path, sep='\t', header=None)
print(len(ol43_merged_inactive))

ol43_merged_inactive = ol43_merged_inactive[[0,1,2]]
ol43_merged_inactive[3] = 'TilingMPRA_ENCSR917SFD_inactive'

tilingMPRA_ENCSR917SFD = pd.concat([ol43_merged_active, ol43_merged_inactive], axis=0, ignore_index=True)
print(len(tilingMPRA_ENCSR917SFD))

print('TilingMPRA OL43')
print(tilingMPRA_ENCSR917SFD.groupby([3]).size())

# -------------------- WHG-STARR-seq --------------------
data_path = os.path.join(project_root, 'data', 'lab_reported_data', 'original_files', 'WHG_STARR_seq', 'ENCFF908UFR.bed.gz')
whg_starr_peaks = pd.read_csv(data_path, sep='\t', header=None, compression='gzip')
print(len(whg_starr_peaks))

whg_starr_peaks = whg_starr_peaks[[0,1,2]]
whg_starr_peaks[3] = 'WHG_STARR_ENCSR661FOW_active'

print('WHG-STARR-seq')
print(len(whg_starr_peaks))


# -------------------- ATAC-STARR-seq --------------------

data_path = os.path.join(project_root, 'data', 'lab_reported_data', 'processed_files', 'ATAC_STARR_seq', 'ATAC_STARR_merged_enhancer_peaks_in_either_orientation.bed')
atac_merged_active = pd.read_csv(data_path, sep='\t', header=None)
print(len(atac_merged_active))

atac_merged_active = atac_merged_active[[0,1,2]]
atac_merged_active[3] = 'ATAC_STARR_ENCSR312UQM_active'

data_path = os.path.join(project_root, 'data', 'lab_reported_data', 'processed_files', 'ATAC_STARR_seq', 'ATAC_STARR_merged_inactive_either_orientation.bed')
atac_merged_inactive = pd.read_csv(data_path, sep='\t', header=None)
print(len(atac_merged_inactive))

atac_merged_inactive = atac_merged_inactive[[0,1,2]]
atac_merged_inactive[3] = 'ATAC_STARR_ENCSR312UQM_inactive'

atac_starr = pd.concat([atac_merged_active, atac_merged_inactive], axis=0, ignore_index=True)
print(len(atac_starr))

print('ATAC-STARR-seq')
print(atac_starr.groupby([3]).size())

# -------------------- Merge All --------------------

# Combine all datasets into one large BED
all_results = pd.concat([lenti, tilingMPRA_ENCSR394HXI, tilingMPRA_ENCSR917SFD, tilingMPRA_ENCSR363XER, 
                         whg_starr_peaks, atac_starr], axis=0, ignore_index=True)
print(len(all_results))

# Save to BED file
output_path = os.path.join(project_root, 'data', 'lab_reported_data', 'processed_files', 'all_merged_results_from_lab_reported_data.bed')
all_results.to_csv(output_path, sep='\t', header=False, index=False)
print('finished merging all results')

# -------------------- Merge TilingMPRA --------------------

# Combine all TilingMPRA results into one large BED
all_tilingMPRA_results = pd.concat([tilingMPRA_ENCSR394HXI, tilingMPRA_ENCSR917SFD, tilingMPRA_ENCSR363XER], axis=0, ignore_index=True)
print(len(all_tilingMPRA_results))

# Save to BED file
output_path = os.path.join(project_root, 'data', 'lab_reported_data', 'processed_files', 'TilingMPRA', 'all_TilingMPRA_results_from_lab_reported_data.bed')
all_tilingMPRA_results.to_csv(output_path, sep='\t', header=False, index=False)
print('finished merging all results')

