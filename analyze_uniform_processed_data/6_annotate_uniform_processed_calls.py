# -------------------- Standard Library Imports --------------------
import argparse
import glob
import os
import sys
from multiprocessing import Pool, cpu_count
from subprocess import PIPE, Popen, STDOUT, call, run

# -------------------- Scientific Computing ------------------------
import numpy as np
import pandas as pd
import scipy

# -------------------- Bioinformatics Libraries --------------------
import pybedtools
import pyBigWig
import pysam

# -------------------- PyBedTools Temp Dir Setup -------------------
# Set temporary directory for PyBedTools operations
pybedtools.helpers.set_tempdir('/fs/cbsuhy02/storage/jz855/tmp/')

# -------------------- Project Root -------------------
# Get the current working directory as the root for relative paths
project_root = os.path.abspath(os.path.join(os.getcwd(), '..'))
sys.path.insert(0, project_root)
from src.utils import *


def get_active_inactive_regions(merged_peak_path, bin_data_path):
    """
    Load and process active and inactive MPRA regions into a single annotated DataFrame.

    Parameters:
    - merged_peak_path: path to BED file of merged active peaks
    - bin_data_path: path to gzipped BED file of all tested bins

    Returns:
    - combined_df: DataFrame containing active and inactive regions with a 'call' column
    """

    # --------------------------------------------
    # Load merged enhancer peaks (active regions)
    # --------------------------------------------
    peak = pd.read_csv(merged_peak_path, sep='\t', header=None)
    print(f"Loaded {len(peak)} merged enhancer peaks")

    active = peak[[0, 1, 2, 3, 4, 5, 6, 13]].copy()
    active.columns = ['chr', 'start', 'end', 'name', 'logFC', 'strand', 'z_score', 'size']
    active = active[~active['chr'].str.contains('chrM|_')]  # Filter chrM and alt contigs
    active['call'] = 'active'
    print(f"Filtered to {len(active)} active peaks")

    # --------------------------------------------
    # Identify inactive regions
    # --------------------------------------------

    # Re-load for clean columns
    peak = peak[[0, 1, 2, 3, 4, 5, 6]].copy()
    peak.columns = ['chr', 'start', 'end', 'name', 'logFC', 'strand', 'z_score']
    peak = peak.sort_values(['chr', 'start'])

    # Load all bins tested
    all_bins = pd.read_csv(bin_data_path, sep='\t', header=None, compression='gzip')
    all_bins = all_bins[[0, 1, 2, 3, 4, 5, 8, 10]]
    all_bins.columns = ['chr', 'start', 'end', 'name', 'logFC', 'strand', 'call_raw', 'FDR']
    all_bins = all_bins.sort_values(['chr', 'start'])

    # BedTool conversion
    all_bins_bed = pybedtools.BedTool.from_dataframe(all_bins)
    peak_bed = pybedtools.BedTool.from_dataframe(peak)

    # Get bins not overlapping any active peak
    inactive_bins_bed = all_bins_bed.intersect(peak_bed, v=True)
    print(f"Found {len(inactive_bins_bed)} inactive bins (non-overlapping)")

    # Merge adjacent bins to create regions
    inactive = inactive_bins_bed.merge().to_dataframe(disable_auto_names=True, header=None)
    inactive.columns = ['chr', 'start', 'end']
    inactive['name'] = [f'inactive_region_{i}' for i in range(1, len(inactive) + 1)]
    inactive['logFC'] = np.nan
    inactive['strand'] = '.'
    inactive['z_score'] = np.nan
    inactive['size'] = inactive['end'] - inactive['start']
    inactive['call'] = 'inactive'

    inactive = inactive[['chr', 'start', 'end', 'name', 'logFC', 'strand', 'z_score', 'size', 'call']]
    inactive = inactive.sort_values(['chr', 'start'])
    inactive = inactive[~inactive['chr'].str.contains('chrM|_')]
    print(f"Final inactive regions: {len(inactive)}")

    # --------------------------------------------
    # Combine active and inactive into one DataFrame
    # --------------------------------------------

    active = active[['chr', 'start', 'end', 'name', 'logFC', 'strand', 'z_score', 'size', 'call']]
    combined_df = pd.concat([active, inactive], ignore_index=True)
    print(f"Combined DataFrame has {len(combined_df)} rows")

    return combined_df

def compute_GROcap_signal(region_info):
    """
    Computes the GRO-cap signal (forward, reverse, and total) for a given genomic region.
    Signals are then normalized by the region size (length in bp).

    Parameters:
    -----------
    region_info : tuple
        A tuple containing (chrom, start, end, name) representing a genomic region.

    Returns:
    --------
    list
        A list containing:
        [name,
         raw_forward_count, raw_reverse_count, raw_total_count,
         normalized_forward_signal, normalized_reverse_signal, normalized_total_signal]
    """

    chrom, start, end, name = region_info

    # Paths to forward and reverse GRO-cap bigWig signal files
    gro_cap_forward_path = os.path.join(project_root, 'data', 'reference', 'K562_GRO_cap', 'K562_GROcap_hg38_aligned_pl.bw')
    gro_cap_reverse_path = os.path.join(project_root, 'data', 'reference', 'K562_GRO_cap', 'K562_GROcap_hg38_aligned_mn.bw')

    # Load bigWig files
    gro_cap_bw_forward = pyBigWig.open(gro_cap_forward_path)
    gro_cap_bw_reverse = pyBigWig.open(gro_cap_reverse_path)

    # Sum of GRO-cap signals across the region
    forward_count = np.nansum(gro_cap_bw_forward.values(chrom, start, end))
    reverse_count = -np.nansum(gro_cap_bw_reverse.values(chrom, start, end))  # reverse is negative strand
    total_count = forward_count + reverse_count

    # Compute region length
    size = end - start

    # Normalize signals by region size
    forward_signal = forward_count / size
    reverse_signal = reverse_count / size
    total_signal = total_count / size

    return [
        name,
        forward_count, reverse_count, total_count,
        forward_signal, reverse_signal, total_signal
    ]

def annotate_GROcap_signal(tested_df):
    """
    Annotate a DataFrame of genomic regions with GRO-cap signal information.

    Parameters:
    -----------
    tested_df : pandas.DataFrame
        Input DataFrame containing columns:
        ['chr', 'start', 'end', 'name', 'logFC', 'strand', 'z_score', 'size', 'call']

    Returns:
    --------
    pandas.DataFrame
        Updated DataFrame with additional GRO-cap columns:
        ['forward_raw', 'reverse_raw', 'total_raw', 
         'forward_signal', 'reverse_signal', 'total_signal']
    """

    # -------------------------------------------
    # Step 1: Prepare list of region info for signal extraction
    # -------------------------------------------
    region_info_list = tested_df[['chr', 'start', 'end', 'name']] \
        .apply(lambda x: [x[0], x[1], x[2], x[3]], axis=1).tolist()

    # -------------------------------------------
    # Step 2: Compute GRO-cap signal in parallel
    # -------------------------------------------
    with Pool(100) as pool:
        signal_list = pool.map(compute_GROcap_signal, region_info_list)

    # -------------------------------------------
    # Step 3: Convert signal list to DataFrame and index by region name
    # -------------------------------------------
    signal_df = pd.DataFrame(
        signal_list,
        columns=[
            'name',
            'forward_raw', 'reverse_raw', 'total_raw',
            'forward_signal', 'reverse_signal', 'total_signal'
        ]
    ).set_index('name')

    # -------------------------------------------
    # Step 4: Join signal data to original region data by 'name'
    # -------------------------------------------
    tested_df = tested_df.set_index('name')
    tested_df = tested_df.join(signal_df, how='left')

    # -------------------------------------------
    # Step 5: Reorder and reset index for final output
    # -------------------------------------------
    tested_df = tested_df.reset_index()
    tested_df = tested_df[
        ['chr', 'start', 'end', 'name', 'logFC', 'strand', 'z_score', 'size', 'call',
         'forward_raw', 'reverse_raw', 'total_raw',
         'forward_signal', 'reverse_signal', 'total_signal']
    ]

    return tested_df

def map_GROcap_elements(tested_df):
    """
    Intersect GRO-cap annotated regions with GRO-cap elements and compute overlap percentages.

    Parameters
    ----------
    tested_df : pd.DataFrame
        DataFrame containing tested genomic regions with GRO-cap signal already annotated.
        Required columns include:
        ['chr', 'start', 'end', 'name', 'logFC', 'strand', 'z_score', 'size', 'call',
         'forward_raw', 'reverse_raw', 'total_raw',
         'forward_signal', 'reverse_signal', 'total_signal']

    Returns
    -------
    pd.DataFrame
        Annotated DataFrame with GRO-cap element overlaps and percent overlap columns.
        Adds:
        ['gro_chr', 'gro_start', 'gro_end', 'gro_type', 'gro_size', 'gro_overlap_bp',
         'pct_region', 'pct_GROcap_element']
    """
    
    # Convert to BedTool
    tested_bed = pybedtools.BedTool.from_dataframe(tested_df)

    # Intersect with GRO-cap elements; keep all tested regions (`woa=True`)
    intersect = tested_bed.intersect(gro_cap_element_bed, wao=True).to_dataframe(disable_auto_names=True, header=None)
    
    # Rename columns (tested_df has 15 columns; GRO-cap adds 5 more + 1 for overlap)
    intersect.columns = (
        ['chr', 'start', 'end', 'name', 'logFC', 'strand', 'z_score', 'size', 'call',
         'forward_raw', 'reverse_raw', 'total_raw',
         'forward_signal', 'reverse_signal', 'total_signal'] +  # existing
        ['gro_chr', 'gro_start', 'gro_end', 'gro_type', 'gro_size',  # from GRO-cap elements
         'gro_overlap_bp']  # number of base pairs overlapping
    )

    intersect['gro_size'] = intersect['gro_size'].replace('.', -1)
    
    # Calculate percent overlap metrics
    intersect['gro_overlap_pct_region'] = intersect['gro_overlap_bp'] / intersect['size'] * 100
    intersect['gro_overlap_pct_GROcap_element'] = intersect['gro_overlap_bp'] / intersect['gro_size'] * 100

    # Keep only the best overlapping GRO-cap element per region
    intersect = intersect.sort_values('gro_overlap_bp', ascending=False)
    intersect = intersect.drop_duplicates(subset='name')

    return intersect

def map_cCRE(tested_df):
    """
    Annotate genomic regions with cCRE (candidate cis-Regulatory Element) overlap information.

    This function takes a DataFrame of tested genomic regions (e.g., active/inactive enhancers),
    assumes GRO-cap signal and element annotations have already been applied, and maps each region
    to overlapping cCRE elements from ENCODE.

    For each region, overlap percentage with cCRE elements is calculated both from the perspective
    of the tested region and the cCRE element. The best overlapping cCRE is kept per region.

    Parameters:
    ----------
    tested_df : pandas.DataFrame
        A DataFrame containing tested genomic regions with prior GRO-cap annotations.
        Expected to have at least the first 15 columns including region info, signal, and size.

    Returns:
    -------
    intersect : pandas.DataFrame
        Annotated DataFrame with overlap info from cCREs, including:
        - Base pair overlap
        - Percent overlap with region and with cCRE
        - Best match per region
    """

    # Convert input DataFrame to BedTool for overlap operations
    tested_bed = pybedtools.BedTool.from_dataframe(tested_df)

    # Intersect with cCRE elements (retain all tested regions even if no overlap: `wao=True`)
    intersect = tested_bed.intersect(cCRE_bed, wao=True).to_dataframe(disable_auto_names=True, header=None)

    # Rename columns:
    # - First 15 columns come from `tested_df`
    # - Next 5 columns come from the cCRE BED file
    # - One column is the number of overlapping base pairs (from `woa=True`)
    intersect.columns = (
        ['chr', 'start', 'end', 'name', 'logFC', 'strand', 'z_score', 'size', 'call',
         'forward_raw', 'reverse_raw', 'total_raw',
         'forward_signal', 'reverse_signal', 'total_signal', 
         'gro_chr', 'gro_start', 'gro_end', 'gro_type', 'gro_size',  # from GRO-cap elements
         'gro_overlap_bp', 
         'gro_overlap_pct_region', 'gro_overlap_pct_GROcap_element'] +  # calculated earlier
        ['ccre_chr', 'ccre_start', 'ccre_end', 'ccre_name', 'ccre_type', 'ccre_size',  # from cCRE BED
         'ccre_overlap_bp']  # base pair overlap with cCRE
    )
    
    # Handle missing cCRE size values returned as '.' and convert to numeric
    intersect['ccre_size'] = intersect['ccre_size'].replace('.', -1)
    intersect['ccre_size'] = pd.to_numeric(intersect['ccre_size'], errors='coerce')
    
    # Calculate percent of the tested region and cCRE element that overlaps
    intersect['ccre_overlap_pct_region'] = intersect['ccre_overlap_bp'] / intersect['size'] * 100
    intersect['ccre_overlap_pct_ccre'] = intersect['ccre_overlap_bp'] / intersect['ccre_size'] * 100

    # Sort by highest cCRE overlap and keep only the top match per tested region
    intersect = intersect.sort_values('ccre_overlap_bp', ascending=False)
    intersect = intersect.drop_duplicates(subset='name')

    return intersect

def annotate_promoter(tested_df, column_name, promoter_region_bed):
    """
    Annotate tested genomic regions with promoter proximity.

    For each region in the input DataFrame, this function determines whether it overlaps
    a provided promoter region (e.g., ±200bp, ±1kb, ±4kb from TSS), using a minimum of 90%
    reciprocal overlap with the region (`f=0.9`). It adds a binary column indicating promoter overlap.

    Parameters:
    ----------
    tested_df : pandas.DataFrame
        DataFrame containing genomic regions with at least ['chr', 'start', 'end', 'name'] columns.
    
    column_name : str
        Name of the new column to be added, indicating promoter overlap ('Y' or 'N').
    
    promoter_region_bed : pybedtools.BedTool
        A BedTool object containing promoter regions to be used for overlap.

    Returns:
    -------
    pandas.DataFrame
        Input DataFrame with an additional binary column named `column_name` marking promoter overlap.
    """

    # Convert the tested regions to BedTool format for overlap analysis
    tested_bed = pybedtools.BedTool.from_dataframe(tested_df[['chr', 'start', 'end', 'name']])
    
    # Preserve the original column order for clean reconstruction later
    original_columns = tested_df.columns.tolist()
    
    # Intersect tested regions with promoter regions using 90% reciprocal overlap (f=0.9)
    intersect = tested_bed.intersect(promoter_region_bed, wa=True, f=0.9).to_dataframe(disable_auto_names=True, header=None)
    intersect = intersect.drop_duplicates()  # Remove duplicate overlaps if any

    # Create a copy of the input DataFrame to annotate
    df = tested_df.copy()
    
    # Initialize new promoter annotation column to 'N' (not proximal)
    df[column_name] = 'N'
    
    # Set index to 'name' for easy row access
    df = df.set_index('name')
    
    # For overlapping regions, set promoter annotation to 'Y'
    if len(intersect) > 0:
        df.loc[intersect[3].tolist(), column_name] = 'Y'
        
    # Restore original indexing and column order + new annotation
    df = df.reset_index()
    df = df[original_columns + [column_name]]
    
    return df


def classify_GROcap_signal_levels(tested_df, output_col, combined=True, low=0.01, high=0.08):
    """
    Annotate GRO-cap signal strength into transcriptional activity classes.

    Parameters:
    ----------
    df : pandas.DataFrame
        Input DataFrame that contains a 'total_signal' column.

    output_col : str
        Name of the column to be created or overwritten to store transcriptional class labels.

    combined : bool (default=True)
        If True, all regions with zero signal are assigned to 'low_transcription'.
        If False, zero-signal regions are labeled as 'not_transcribed'.

    low : float (default=0.01)
        Threshold for minimum signal to be considered 'low_transcription' (if not combined).

    high : float (default=0.08)
        Threshold above which regions are considered 'high_transcription'.

    Returns:
    -------
    pandas.DataFrame
        The same DataFrame with an added or updated `output_col` containing:
        ['none_transcription', 'low_transcription', 'medium_transcription', 'high_transcription']
    """
    signal_col = 'total_signal'
    
    df = tested_df.copy()
    
    if combined:
        # All regions start as 'low_transcription' if combined mode is enabled
        df[output_col] = 'low_transcription'
    else:
        # If not combined, initialize with 'not_transcribed' for zero or missing signal
        df[output_col] = 'none_transcription'
        df.loc[df[(df[signal_col].map(lambda x: x <= low and x > 0))].index.tolist(), output_col] = 'low_transcription'
    
    # Classify medium and high transcription regions
    df.loc[df[(df[signal_col].map(lambda x: x > low and x <= high))].index.tolist(), output_col] = 'medium_transcription'
    df.loc[df[(df[signal_col].map(lambda x: x > high))].index.tolist(), output_col] = 'high_transcription'

    return df


def classify_binary_GROcap_signal(tested_df, output_col):
    """
    Classify regions based on whether GRO-cap signal is present or absent.

    Parameters:
    ----------
    df : pandas.DataFrame
        Input DataFrame with a 'total_signal' column.

    output_col : str
        Name of the column to store binary transcription call.

    Returns:
    -------
    pandas.DataFrame
        The same DataFrame with a new column `output_col` containing:
        ['not_transcribed', 'transcribed']
    """
    signal_col = 'total_signal'
    
    df = tested_df.copy()
    
    # Default to not transcribed
    df[output_col] = 'not_transcribed'
    
    # Update regions with any non-zero signal to 'transcribed'
    df.loc[df[(df[signal_col].map(lambda x: x > 0))].index.tolist(), output_col] = 'transcribed'

    return df

def classify_overlap_extent_with_elements(tested_df, element_type, element_col, element):
    """
    Classify the extent of overlap between tested genomic regions and a specified regulatory element
    (e.g., GRO-cap element types or cCRE types), based on reciprocal overlap criteria.

    Parameters:
    -----------
    tested_df : pandas.DataFrame
        DataFrame of tested genomic regions that have already been annotated with overlap metrics.

    element_type : str
        The specific type of element to classify overlap for (e.g., 'bidirectional', 'PLS').
        Special types like 'ELS' or 'ELS_PLS' are merged categories (e.g., dELS + pELS).

    element_col : str
        The column in tested_df that stores the type of the overlapped element 
        (e.g., 'gro_type' or 'ccre_type').

    element : str
        Which regulatory element set to use for overlap scoring.
        Should be one of:
            - 'GROcap': uses gro_overlap_pct_region and gro_overlap_pct_GROcap_element
            - 'cCRE': uses ccre_overlap_pct_region and ccre_overlap_pct_ccre

    Returns:
    --------
    pandas.DataFrame
        A copy of the input DataFrame with a new column named `element_type`, 
        where values represent:
            - 'high': ≥80% reciprocal overlap
            - 'moderate': ≥50% reciprocal overlap
            - 'low': any remaining overlap
            - '.': no overlap
            - or 'not_overlap' (only when element_type includes 'not_overlap')
    """
    
    df = tested_df.copy()
    df = df.set_index('name')
    original_columns = tested_df.columns.tolist()

    # Handle special case: annotate regions NOT overlapping the element set
    if 'not_overlap' in element_type:
        df[element_type] = element_type
        df.loc[df[df[element_col].map(lambda x: x != '.')].index.tolist(), element_type] = '.'
    
    else:
        # -----------------------------------------------
        # Determine subset of rows with the given element_type
        # For merged classes like ELS or ELS_PLS, combine multiple element labels
        # -----------------------------------------------
        if element_type == 'ELS':
            overlap_with_element = df[df[element_col].map(lambda x: x in ['dELS', 'pELS'])]
        elif element_type == 'ELS_PLS':
            overlap_with_element = df[df[element_col].map(lambda x: x in ['dELS', 'pELS', 'PLS'])]
        else:
            overlap_with_element = df[df[element_col] == element_type]
        
        # -----------------------------------------------
        # Apply element-specific reciprocal overlap thresholds
        # -----------------------------------------------
        if element == 'GROcap':
            reciprocal_50pct = overlap_with_element[
                (overlap_with_element['gro_overlap_pct_region'] >= 50) &
                (overlap_with_element['gro_overlap_pct_GROcap_element'] >= 50)
            ]
            reciprocal_80pct = overlap_with_element[
                (overlap_with_element['gro_overlap_pct_region'] >= 80) &
                (overlap_with_element['gro_overlap_pct_GROcap_element'] >= 80)
            ]

        elif element == 'cCRE':
            reciprocal_50pct = overlap_with_element[
                (overlap_with_element['ccre_overlap_pct_region'] >= 50) &
                (overlap_with_element['ccre_overlap_pct_ccre'] >= 50)
            ]
            reciprocal_80pct = overlap_with_element[
                (overlap_with_element['ccre_overlap_pct_region'] >= 80) &
                (overlap_with_element['ccre_overlap_pct_ccre'] >= 80)
            ]

        # -----------------------------------------------
        # Annotate overlap strength by defaulting to '.', then upgrading
        # -----------------------------------------------
        element_type = '_'.join(element_type.split('-'))
        df[element_type] = '.'
        df.loc[overlap_with_element.index.tolist(), element_type] = 'low'
        df.loc[reciprocal_50pct.index.tolist(), element_type] = 'moderate'
        df.loc[reciprocal_80pct.index.tolist(), element_type] = 'high'

    # Restore original ordering
    df = df.reset_index()
    df = df[original_columns + [element_type]]

    return df

def run_annotation_pipeline(merged_peak_path, bin_data_path):
    """
    Annotate active and inactive genomic regions with:
    - GRO-cap signal
    - GRO-cap element overlap
    - cCRE overlap
    - Promoter proximity
    - Transcription activity (signal-based)
    - Regulatory element classification (GRO-cap & cCRE)

    Parameters:
    -----------
    merged_peak_path : str
        Path to the merged enhancer peak file (BED format).

    bin_data_path : str
        Path to the all-bin tested result file (gzipped BED-like format).

    Returns:
    --------
    pandas.DataFrame
        Annotated dataframe of tested genomic regions.
    """

    # Step 1: Get tested regions (active + inactive)
    tested_df = get_active_inactive_regions(merged_peak_path, bin_data_path)

    # Step 2: Annotate GRO-cap signal from bigWig
    tested_df = annotate_GROcap_signal(tested_df)

    # Step 3: Map GRO-cap elements
    tested_df = map_GROcap_elements(tested_df)

    # Step 4: Map cCRE elements
    tested_df = map_cCRE(tested_df)

    # Step 5: Annotate promoter proximity
    tested_df = annotate_promoter(tested_df, 'promoter_400bp_tss', promoter_400bp_bed)
    tested_df = annotate_promoter(tested_df, 'promoter_1kb_tss', promoter_1kb_bed)
    tested_df = annotate_promoter(tested_df, 'promoter_4kb_tss', promoter_4kb_bed)

    # Step 6: Transcription classification (GRO-cap signal)
    tested_df = classify_GROcap_signal_levels(tested_df, 'separate_GROcap_signal_levels', combined=False)
    tested_df = classify_binary_GROcap_signal(tested_df, 'binary_transcription_class')

    # Step 7: GRO-cap element overlap classification
    for element_type in ['divergent', 'unidirectional', 'bidirectional', 'not_overlap_GROcap_elements']:
        tested_df = classify_overlap_extent_with_elements(
            tested_df,
            element_type=element_type,
            element_col='gro_type',
            element='GROcap'
        )

    # Step 8: cCRE element overlap classification
    cCRE_type_list = [
        'CA-CTCF', 'CA-H3K4me3', 'CA-TF', 'dELS', 'pELS',
        'Low-DNase', 'PLS', 'CA-only', 'not_overlap_cCRE',
        'ELS', 'ELS_PLS'
    ]

    for element_type in cCRE_type_list:
        tested_df = classify_overlap_extent_with_elements(
            tested_df,
            element_type=element_type,
            element_col='ccre_type',
            element='cCRE'
        )

    return tested_df

# Specify root directory
project_root = os.path.abspath(os.path.join(os.getcwd(), '..'))
print(project_root)

# -------------------- Load Reference Files -------------------------
print("Loading reference annotation files...")

gro_cap_element_bed = pybedtools.BedTool(os.path.join(project_root, 'data', 'reference', 'K562_GRO_cap', 'GRO_cap_element', 'processed_files', 'GROcap_elements_all.bed'))
cCRE_bed = pybedtools.BedTool(os.path.join(project_root, 'data', 'reference', 'cCRE_v4', 'processed_files', 'K562_cCRE_v4_processed.bed'))
promoter_400bp_bed = pybedtools.BedTool(os.path.join(project_root, 'data', 'reference', 'hg38', 'GENCODE_v45', 'processed_files', 'promoters_400b_GENCODE_v45_protein_coding_tss_centered.bed.gz'))
promoter_1kb_bed = pybedtools.BedTool(os.path.join(project_root, 'data', 'reference', 'hg38', 'GENCODE_v45', 'processed_files', 'promoters_1kb_GENCODE_v45_protein_coding_tss_centered.bed.gz'))
promoter_4kb_bed = pybedtools.BedTool(os.path.join(project_root, 'data', 'reference', 'hg38', 'GENCODE_v45', 'processed_files', 'promoters_4kb_GENCODE_v45_protein_coding_tss_centered.bed.gz'))

# -------------------- Define Datasets for Annotation -------------------------

lenti_datasets = [
    ('LentiMPRA', 'element_level_all_result.bed.gz', 'merged_enhancer_peak_in_either_orientation.bed.gz', 'annotated_tested_regions_either_orientation.bed.gz'),
    ('LentiMPRA', 'all_element_tested_in_both_orientations.bed.gz', 'merged_enhancer_peak_orientation_independent.bed.gz', 'annotated_tested_regions_both_orientations.bed.gz')
]

tiling_datasets = [
    'OL13_ENCSR394HXI', 'OL43_ENCSR917SFD', 'OL45_ENCSR363XER'
]

starr_datasets = [
    'ATAC_STARR_seq', 'WHG_STARR_seq'
]

# -------------------- Annotation for LentiMPRA -------------------------
print("\nRunning annotation for LentiMPRA...")
for i, (assay, bin_file, peak_file, output_file) in enumerate(lenti_datasets):
    print(f"\n[Step {i+1}] Annotating {output_file.split('_')[-2]} orientation for {assay}...")

    set_dir(os.path.join(project_root, 'data', 'uniform_processed_data', assay, 'annotated_results'))
    bin_data_path = os.path.join(project_root, 'data', 'uniform_processed_data', assay, 'element_level', bin_file)
    merged_peak_path = os.path.join(project_root, 'data', 'uniform_processed_data', assay, 'merged_peak', peak_file)

    annotated = run_annotation_pipeline(merged_peak_path, bin_data_path)
    print(f"Annotated {len(annotated)} regions")

    output_path = os.path.join(project_root, 'data', 'uniform_processed_data', assay, 'annotated_results', output_file)
    annotated.to_csv(output_path, sep='\t', header=False, index=False, compression='gzip')

# -------------------- Annotation for TilingMPRA -------------------------
print("\nRunning annotation for TilingMPRA...")
for dataset in tiling_datasets:
    print(f"\nAnnotating dataset: {dataset} (either orientation)")

    set_dir(os.path.join(project_root, 'data', 'uniform_processed_data', 'TilingMPRA', dataset, 'annotated_results'))

    bin_data_path = os.path.join(project_root, 'data', 'uniform_processed_data', 'TilingMPRA', dataset, 'element_level', 'element_level_all_result.bed.gz')
    merged_peak_path = os.path.join(project_root, 'data', 'uniform_processed_data', 'TilingMPRA', dataset, 'merged_peak', 'merged_enhancer_peak_in_either_orientation.bed.gz')

    annotated_either = run_annotation_pipeline(merged_peak_path, bin_data_path)
    print(f"Annotated {len(annotated_either)} regions (either orientation)")

    output_path = os.path.join(project_root, 'data', 'uniform_processed_data', 'TilingMPRA', dataset, 'annotated_results', 'annotated_tested_regions_either_orientation.bed.gz')
    annotated_either.to_csv(output_path, sep='\t', header=False, index=False, compression='gzip')

    if 'OL13' in dataset:
        print(f"\nAnnotating dataset: {dataset} (both orientations)")

        bin_data_path = os.path.join(project_root, 'data', 'uniform_processed_data', 'TilingMPRA', dataset, 'element_level', 'all_element_tested_in_both_orientations.bed.gz')
        merged_peak_path = os.path.join(project_root, 'data', 'uniform_processed_data', 'TilingMPRA', dataset, 'merged_peak', 'merged_enhancer_peak_orientation_independent.bed.gz')

        annotated_both = run_annotation_pipeline(merged_peak_path, bin_data_path)
        print(f"Annotated {len(annotated_both)} regions (both orientations)")

        output_path = os.path.join(project_root, 'data', 'uniform_processed_data', 'TilingMPRA', dataset, 'annotated_results', 'annotated_tested_regions_both_orientations.bed.gz')
        annotated_both.to_csv(output_path, sep='\t', header=False, index=False, compression='gzip')

# -------------------- Annotation for STARR-seq datasets -------------------------
print("\nRunning annotation for STARR-seq datasets...")
for dataset in starr_datasets:
    print(f"\nAnnotating dataset: {dataset} (either orientation)")

    set_dir(os.path.join(project_root, 'data', 'uniform_processed_data', dataset, 'annotated_results'))

    bin_data_path = os.path.join(project_root, 'data', 'uniform_processed_data', dataset, 'bin_level', 'bin_level_all_result.bed.gz')
    merged_peak_path = os.path.join(project_root, 'data', 'uniform_processed_data', dataset, 'merged_peak', 'merged_enhancer_peak_in_either_orientation.bed.gz')

    annotated_either = run_annotation_pipeline(merged_peak_path, bin_data_path)
    print(f"Annotated {len(annotated_either)} regions (either orientation)")

    output_path = os.path.join(project_root, 'data', 'uniform_processed_data', dataset, 'annotated_results', 'annotated_tested_regions_either_orientation.bed.gz')
    annotated_either.to_csv(output_path, sep='\t', header=False, index=False, compression='gzip')

    print(f"\nAnnotating dataset: {dataset} (both orientations)")

    bin_data_path = os.path.join(project_root, 'data', 'uniform_processed_data', dataset, 'bin_level', 'all_bin_tested_in_both_orientations.bed.gz')
    merged_peak_path = os.path.join(project_root, 'data', 'uniform_processed_data', dataset, 'merged_peak', 'merged_enhancer_peak_orientation_independent.bed.gz')

    annotated_both = run_annotation_pipeline(merged_peak_path, bin_data_path)
    print(f"Annotated {len(annotated_both)} regions (both orientations)")

    output_path = os.path.join(project_root, 'data', 'uniform_processed_data', dataset, 'annotated_results', 'annotated_tested_regions_both_orientations.bed.gz')
    annotated_both.to_csv(output_path, sep='\t', header=False, index=False, compression='gzip')

print("\nAll annotation tasks completed successfully!")
