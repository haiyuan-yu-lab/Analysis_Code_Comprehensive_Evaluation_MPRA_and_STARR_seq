# --------------------------------------------------
# Standard Library Imports
# --------------------------------------------------
import os
import sys
import glob
import pandas as pd
import numpy as np
import argparse
import pysam
import pybedtools
import shutil
import statsmodels.stats.multitest as smm
import shutil

from subprocess import call, PIPE, run, Popen, STDOUT
from multiprocessing import Pool, cpu_count

# --------------------------------------------------
# Local Module Imports
# --------------------------------------------------
# project_root = os.path.abspath(os.path.join(os.getcwd(), '..', '..', '..'))
src_path = os.path.join('/fs/cbsuhy01/storage/jz855/STARR_seq_code/public_analysis_code/src/')
sys.path.insert(0, src_path)

import utils

def merge_bins(data):
    """
    Merge overlapping or adjacent bins from BED-like input.

    Args:
        data (pd.DataFrame): Input BED dataframe with at least ['seqnames', 'start', 'end']

    Returns:
        pd.DataFrame: Merged, sorted BED intervals (columns 0,1,2)
    """
    data = data[['seqnames', 'start', 'end']]
    data = data.sort_values(['seqnames', 'start'])
    bed = pybedtools.BedTool.from_dataframe(data)
    merged = bed.merge()  # merge overlapping genomic regions
    merged = merged.to_dataframe(disable_auto_names=True, header=None)
    merged = merged.sort_values([0, 1])
    return merged

def _build_subject_for_mode(raw, mode):
    """
    Construct the 'subject' dataframe with a known column order so that
    pybedtools intersection column indices are stable.

    Returns a dataframe with columns:
      - 'both'  : ['seqnames', 'start', 'end', 'name', 'activity', 'strand', 
                       'activity_z_score', 'call', 
                       'logFC_for', 'P.Value_for', 'q_for', 'z_score_for', 'neg_ctrl_for', 'call_for',  
                       'logFC_rev', 'P.Value_rev', 'q_rev', 'z_score_rev', 'neg_ctrl_rev', 'call_rev']
      - 'either': ['seqnames', 'start', 'end', 'name', 'logFC', 'strand', 
                       'P.Value', 'adj.P.Val', 'z_score', 'neg_ctrl', 'call']
    """

    if mode == 'both':
        subject = raw[['seqnames', 'start', 'end', 'name', 'activity', 'strand', 
                       'activity_z_score', 'call', 
                       'logFC_for', 'P.Value_for', 'q_for', 'z_score_for', 'neg_ctrl_for', 'call_for',  
                       'logFC_rev', 'P.Value_rev', 'q_rev', 'z_score_rev', 'neg_ctrl_rev', 'call_rev']].sort_values(['seqnames', 'start'])
    elif mode == 'either':
        subject = raw[['seqnames', 'start', 'end', 'name', 'logFC', 'strand', 
                       'P.Value', 'adj.P.Val', 'z_score', 'neg_ctrl', 'call']].sort_values(['seqnames', 'start'])
    else:
        raise ValueError("mode must be 'both' or 'either'")
    return subject

def _calc_activity_for_merged_peak(args):
    """
    Worker for get_activity_after_merging(...).

    Args:
        args (tuple): (idx, merged_df, raw_df)
            - idx (int): Row index of merged peak
            - merged_df (pd.DataFrame): Merged BED intervals
            - raw_df (pd.DataFrame): Raw input bin-level data
            - mode:   'both' or 'either'

    Returns:
      [idx, mean_logFC, z_score_for_mean_logFC, logFC_list,
       summit_seqnames, summit_start, summit_end, summit_logFC, summit_z_score]
    """

    idx, merged, raw, mode = args

    # Select the peak region
    data = merged.loc[idx, :].to_frame().transpose()
    c = data[0].tolist()[0]
    data = pybedtools.BedTool.from_dataframe(data)

    # Filter raw data for same chromosome
    raw = raw[raw['seqnames'] == c]
    raw = raw.sort_values(['seqnames', 'start'])
    raw = pybedtools.BedTool.from_dataframe(raw)

    # Intersect peak with raw bins (full overlap only)
    intersect = data.intersect(raw, wao=True, F=1)
    intersect = intersect.to_dataframe(disable_auto_names=True, header=None)
    intersect = intersect[intersect[intersect.columns.tolist()[-1]] > 0] # keep overlaps only
    intersect = intersect.drop_duplicates([6]) # remove duplicate hits

    if(mode == 'both'):
        # Combine forward and reverse logFC columns
        activity_list = intersect[11].tolist() + intersect[17].tolist()
        avg_activity = np.mean(activity_list)
        # Z-score using negative control distribution
        avg_activity_z_score = (avg_activity-mean_neg_ctrl)/std_neg_ctrl

        # Determine summit (bin with max average logFC)
        intersect = intersect.sort_values(7, ascending=False)
        highest_summit_idx = intersect.index.tolist()[0]
    
    elif(mode == 'either'):

        # Calculate activity
        avg_activity = np.mean(intersect[7].tolist())
        avg_activity_z_score = (avg_activity-mean_neg_ctrl)/std_neg_ctrl
        activity_list = intersect[7].tolist()

        # Get the summit (bin with highest logFC)
        intersect = intersect.sort_values(7, ascending=False)
        highest_summit_idx = intersect.index.tolist()[0]

    return([idx, avg_activity, avg_activity_z_score, activity_list, 
            intersect.loc[highest_summit_idx, 3], # summit seqname
            intersect.loc[highest_summit_idx, 4], # summit start
            intersect.loc[highest_summit_idx, 5], # summit end
            intersect.loc[highest_summit_idx, 7], # summit logFC
            intersect.loc[highest_summit_idx, 12] # summit z-score
            ])

def get_activity_after_merging(merged, raw, mode):

    raw = _build_subject_for_mode(raw, mode)
    arg_list = [[idx, merged, raw, mode] for idx in merged.index.tolist()]

    with Pool(10) as pool:
        results = pool.map(_calc_activity_for_merged_peak, arg_list)

    df = pd.DataFrame(results, columns=['idx', 'mean_logFC', 'z_score_for_mean_logFC', 'logFC_list', 
                                        'summit_seqnames', 'summit_start', 'summit_end', 'summit_logFC', 'summit_z_score'])
    df = df.set_index('idx')

    merged = merged.join(df, how='outer')
    return(merged)

def intersect_neg(query, subject, overlap_size=100):
    
    """
    Intersect two BED-like DataFrames and retain entries from `query` that 
    overlap `subject` by exactly `overlap_size` base pairs.

    Args:
        query (pd.DataFrame): BED-like DataFrame (with at least 3 columns: chrom, start, end)
        subject (pd.DataFrame): BED-like DataFrame to intersect with
        overlap_size (int): Required overlap size to retain matches (default: 100 bp)

    Returns:
        pd.DataFrame: Subset of intersected entries with exact overlap of `overlap_size`
    """
    
    # Convert both DataFrames to BedTool objects
    query = pybedtools.BedTool.from_dataframe(query)
    subject = pybedtools.BedTool.from_dataframe(subject)
    
    # Perform intersection with "write all overlaps" (wao)
    vs = query.intersect(subject, wao=True)
    # Convert back to DataFrame
    vs = vs.to_dataframe(disable_auto_names=True, header=None)
    # Filter rows with exact overlap of specified size
    vs = vs[vs[vs.columns.tolist()[-1]] == overlap_size]
    # Drop duplicate entries based on the query fragment (assumes first 5 columns describe query)
    vs = vs.drop_duplicates([0,1,2,3,4])
    
    return(vs)

# === Argument Parser ===
parser = argparse.ArgumentParser()
parser.add_argument('-i', "--input_path", required=True, help="Diectory to Limma-Voom output")
parser.add_argument('-o', "--out_path", required=True, help="Directory to save the output")
parser.add_argument('-d', "--num_rep_DNA", type=int, required=True, help="Number of DNA replicates")
parser.add_argument('-r', "--num_rep_RNA", type=int, required=True, help="Number of RNA replicates")
parser.add_argument('-c', '--filtered_raw_count', help="Directory to filtered raw count file")
parser.add_argument('-m', '--mode', choices=['both', 'either'], required=True, help="Processing mode: 'both' for orientation-independent call, 'either' for either orientation calls")
parser.add_argument('--neg_ctrl_ref', help="Path to file listing negative control references, required for genome-wide STARR-seq")
parser.add_argument('--data_type', choices=['STARR', 'MPRA'], required=True, help="Processing type: 'STARR' for STARR-seq, 'MPRA' for MPRA")

# Parse arguments
args = parser.parse_args()

# Set current working directory to input directory
os.chdir(args.input_path)

# Set output directory
utils.set_dir(args.out_path)
utils.set_dir(args.out_path+'bin_level/')
utils.set_dir(args.out_path+'peak_level/')

# Set temporary directory for pybedtools
utils.set_dir(args.out_path+'tmp_pybedtools/')

# Set temporary directory for pybedtools
pybedtools.helpers.set_tempdir(args.out_path+'tmp_pybedtools/')

# =============================================
# Load limma output table
limma_out = pd.read_csv('./limma_out.txt', sep='\t')
print('limma_out '+str(len(limma_out)))

if(args.data_type == 'STARR'):
    col_list = limma_out.columns.tolist()
    limma_out['name'] = ['bin_'+str(i) for i in range(1, len(limma_out)+1)]
    limma_out['bin_id'] = ['bin_'+str(i) for i in range(1, len(limma_out)+1)]
    
    limma_out = limma_out.set_index(['name'])

    # Annotate bins that fall into negative control reference regions
    
    # Load negative control regions
    neg_ctrl_ref = pd.read_csv(args.neg_ctrl_ref, sep='\t', header=None)
    print('Number of negative control regions '+ str(len(neg_ctrl_ref)))

    # Process forward strand bins
    forward = limma_out[limma_out['strand'] == '+'][['seqnames', 'start', 'end', 'bin_id', 'strand']]
    # neg = neg_ctrl_ref[neg_ctrl_ref[4] == '+']
    forward = intersect_neg(forward, neg_ctrl_ref)

    # Process reverse strand bins
    reverse = limma_out[limma_out['strand'] == '-'][['seqnames', 'start', 'end', 'bin_id', 'strand']]
    # neg = neg_ctrl_ref[neg_ctrl_ref[4] == '-']
    reverse = intersect_neg(reverse, neg_ctrl_ref)

    # Combine negative control overlaps
    neg_ctrl_region = pd.concat([forward, reverse], ignore_index=True, axis=0)
    print(len(neg_ctrl_region))

    # Annotate limma output with negative control flag
    limma_out['neg_ctrl'] = 'N'
    limma_out.loc[neg_ctrl_region[3].tolist(), 'neg_ctrl'] = 'Y'

    print('Check number of negative controls and regions of interests in limma_out')
    print(limma_out.groupby(['neg_ctrl']).size())

    # Save negative control region activity
    neg_ctrl_region = limma_out[limma_out['neg_ctrl'] == 'Y'].reset_index()
    neg_ctrl_region = neg_ctrl_region[['seqnames', 'start' , 'end', 'name', 'strand', 'neg_ctrl', 'logFC', 'P.Value', 'adj.P.Val']]
    neg_ctrl_region.to_csv(args.out_path+'neg_ctrl_region.txt', sep='\t', index=False, header=True)
    print("finished saving negative control regions")

    limma_out = limma_out.reset_index()
    limma_out = limma_out[col_list+['neg_ctrl', 'name']]

elif(args.data_type == 'MPRA'):
    # Change column names to be consistent with STARR-seq
    limma_out.columns = ['seqnames', 'start', 'end', 'name', 'strand', 'ori_annot', 'annot', 'logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val', 'B']

    # Annotate negative controls
    limma_out['neg_ctrl'] = 'N'
    limma_out.loc[limma_out[limma_out['annot'] == 'neg_ctrl'].index.tolist(), 'neg_ctrl'] = 'Y'

    print('Check number of negative controls and regions of interests in limma_out')
    print(limma_out.groupby(['neg_ctrl']).size())

    # Save negative control region activity
    neg_ctrl_region = limma_out[limma_out['neg_ctrl'] == 'Y']
    neg_ctrl_region = neg_ctrl_region[['seqnames', 'start', 'end', 'name', 'strand', 'neg_ctrl', 'logFC', 'P.Value', 'adj.P.Val']]
    neg_ctrl_region.to_csv(args.out_path+'neg_ctrl_region.txt', sep='\t', index=False, header=True)
    print("finished saving negative control regions")

# Compute Z-scores using negative control bins

mean_neg_ctrl = np.mean(neg_ctrl_region['logFC'].values)
print(mean_neg_ctrl)
std_neg_ctrl = np.std(neg_ctrl_region['logFC'].values)
print(std_neg_ctrl)

# Set threshold for logFC based on Z=1.96
logFC_cutoff = 1.96*std_neg_ctrl+mean_neg_ctrl
print('logFC cutoff')
print(logFC_cutoff)

# Focus only genomic sequences
limma_out = limma_out[(limma_out['seqnames'] != '.')]
limma_out = limma_out[(limma_out['start'] > 0)]

# Calculate z-score for all regions
limma_out['z_score'] = (limma_out['logFC']-mean_neg_ctrl)/std_neg_ctrl

# make enhancer or repressor call based on logFC and adjusted p-value for all regions tested
limma_out['call'] = 'inactive'
idx_list = limma_out[(limma_out['logFC'] >= logFC_cutoff) & 
                     (limma_out['adj.P.Val'] < 0.05)].index.tolist()
limma_out.loc[idx_list, 'call'] = 'enhancer'
idx_list = limma_out[(limma_out['logFC'] <= -logFC_cutoff) & 
                     (limma_out['adj.P.Val'] < 0.05)].index.tolist()
limma_out.loc[idx_list, 'call'] = 'repressor'

print('Check number of original calls')
print(limma_out.groupby(['call']).size())

# Save bin-level/element-level data
file_path = os.path.join(args.out_path, "bin_level", "all_results_tested_in_either_orientation.bed.gz")

if os.path.exists(file_path):
    print(f"File exists: {file_path}")
    print(f"Continue to obtain peak level results")
else:
    print(f"File NOT found: {file_path}")
    print(f"Save Bin-level results first")

    if(args.data_type == 'MPRA'):
        limma_out = limma_out.set_index('name')
        limma_out = limma_out[['seqnames', 'start', 'end', 'strand','logFC', 'P.Value', 'adj.P.Val', 'neg_ctrl', 'z_score', 'call']]

    elif(args.data_type == 'STARR'):
        # Keep only selected columns and reformat index
        limma_out = limma_out.set_index(['seqnames', 'start', 'end', 'strand'])
        limma_out = limma_out[['logFC', 'P.Value', 'adj.P.Val', 'name', 'neg_ctrl', 'z_score', 'call']]

    if(args.filtered_raw_count != None):
        count = pd.read_csv(args.filtered_raw_count, sep='\t')
    else:
        count = pd.read_csv('./filtered_raw_count.txt', sep='\t')

    if(args.data_type == 'MPRA'):
        count.columns = ['seqnames', 'start', 'end', 'name', 'strand', 'orig_annot', 'annot'] + ['DNA{0}'.format(i) for i in range(1, args.num_rep_DNA+1)] + ['RNA{0}'.format(i) for i in range(1, args.num_rep_RNA+1)]
        count = count.set_index('name')[['DNA{0}'.format(i) for i in range(1, args.num_rep_DNA+1)] + ['RNA{0}'.format(i) for i in range(1, args.num_rep_RNA+1)]]

    elif(args.data_type == 'STARR'):
        count = count.set_index(['seqnames', 'start', 'end', 'strand'])
        print('filtered raw count '+ str(len(count)))

    # Join raw counts to annotated bin-level results
    all_results = limma_out.join(count, how='left')
    all_results = all_results.reset_index()
    # Exclude non-genomic sequences
    all_results = all_results[all_results['seqnames'] != '.']
    print('number of genomic sequences tested')
    print(len(all_results))
    all_results = all_results[['seqnames', 'start', 'end', 'name', 'logFC', 'strand', 'P.Value', 'adj.P.Val', 'z_score', 'neg_ctrl', 'call']+['DNA{0}'.format(i) for i in range(1, args.num_rep_DNA+1)] + ['RNA{0}'.format(i) for i in range(1, args.num_rep_RNA+1)]]

    # Save
    all_results.to_csv(file_path, sep='\t', index=False, header=False, compression='gzip')
    print('finished saving bin-level results in each orientation')

if(args.mode == 'either'):
    all_enhancers = limma_out[limma_out['call'] == 'enhancer'].reset_index()
    all_repressors = limma_out[limma_out['call'] == 'repressor'].reset_index()

    for data, fname in zip([all_enhancers, all_repressors], 
                           ['enhancer_peak_in_either_orientation.bed.gz', 'repressor_peak_in_either_orientation.bed.gz']):
        
        # Merge adjacent/overlapping bins/elements into peak regions
        if(len(data) == 0):
            continue
        merged = merge_bins(data)
        merged = get_activity_after_merging(merged, data, args.mode)

        # Add required columns to adjust format
        merged['strand'] = '.'
        merged['size'] = merged[2].values-merged[1].values

        merged.columns = ['seqnames', 'start', 'end', 'mean_logFC', 'z_score_for_mean_logFC', 'logFC_list', 
                          'summit_seqnames', 'summit_start', 'summit_end', 'summit_logFC', 'summit_z_score', 'strand', 'size']
        merged = merged.sort_values(['seqnames', 'start'])

        # Add peak names and re-arrange columns
        merged['name'] = ['peak{0}'.format(i) for i in range(1, len(merged)+1)]
        merged = merged[['seqnames', 'start', 'end', 'name', 'mean_logFC', 'strand', 
                         'z_score_for_mean_logFC', 'logFC_list', 
                         'summit_seqnames', 'summit_start', 'summit_end', 'summit_logFC', 'summit_z_score', 'size']]
        
        print('Number of peaks in '+fname+' '+str(len(merged)))

        # Write merged peak file
        merged.to_csv(args.out_path+'peak_level/'+fname, sep='\t', index=False, header=False, compression='gzip')
        print('finished saving ' + fname)

elif(args.mode == 'both'):

    limma_out = limma_out.reset_index()

    # Separate forward and reverse calls
    forward = limma_out[limma_out['strand'] == '+']
    forward = forward.set_index(['seqnames', 'start', 'end'])
    forward = forward[['logFC', 'P.Value', 'z_score', 'neg_ctrl', 'call']]
    forward.columns = [x+'_for' for x in forward.columns.tolist()]
    
    reverse = limma_out[limma_out['strand'] == '-']
    reverse = reverse.set_index(['seqnames', 'start', 'end'])
    reverse = reverse[['logFC', 'P.Value', 'z_score', 'neg_ctrl', 'call']]
    reverse.columns = [x+'_rev' for x in reverse.columns.tolist()]

    # Join both orientations and keep only complete pairs
    both = forward.join(reverse, how='outer')
    print('number of regions after merging orientatoins')
    print(len(both))
    both = both.dropna()
    print('number of regions tested in both orientations')
    print(len(both))

    # Average activity and compute z-score
    both['activity'] = (both['logFC_for'].values + both['logFC_rev'].values)/2
    both['activity_z_score'] = (both['activity']-mean_neg_ctrl)/std_neg_ctrl
    
    # Set up bin name and strand
    both = both.reset_index()
    both = both.sort_values(['seqnames', 'start'])
    both['name'] = ['bin_tested_both_{0}'.format(i) for i in range(1, len(both)+1)]
    both['strand'] = '.'
    both = both.set_index('name')

    # Recalculate FDRs for both strands
    tmp1 = both[['P.Value_for']]
    tmp1.columns = ['p']
    tmp1['type'] = 'for'
    tmp1 = tmp1.reset_index()
    
    tmp2 = both[['P.Value_rev']]
    tmp2.columns = ['p']
    tmp2['type'] = 'rev'
    tmp2 = tmp2.reset_index()
    
    tmp = pd.concat([tmp1, tmp2], axis=0)
    tmp['adj.p'] = smm.fdrcorrection(tmp['p'].values, alpha=0.05, method='indep', is_sorted=False)[1]

    # Split corrected p-values
    tmp1 = tmp[tmp['type'] == 'for']
    tmp1 = tmp1.set_index('name')
    tmp1 = tmp1[['adj.p']]
    tmp1.columns = ['q_for']
    
    tmp2 = tmp[tmp['type'] == 'rev']
    tmp2 = tmp2.set_index('name')
    tmp2 = tmp2[['adj.p']]
    tmp2.columns = ['q_rev']

    # Join back to main
    both = both.join(tmp1, how='inner')
    print(len(both))
    both = both.join(tmp2, how='inner')
    print(len(both))
    
    both = both.reset_index()

    # Final enhancer/repressor calls (requires both strands to be significant)
    both['call'] = 'inactive'
    both.loc[both[(both['logFC_for'] >= logFC_cutoff) & (both['q_for'] < 0.05) & 
                  (both['logFC_rev'] >= logFC_cutoff) & (both['q_rev'] < 0.05)].index.tolist(), 'call'] = 'enhancer'
    both.loc[both[(both['logFC_for'] <= -logFC_cutoff) & (both['q_for'] < 0.05) & 
                  (both['logFC_rev'] <= -logFC_cutoff) & (both['q_rev'] < 0.05)].index.tolist(), 'call'] = 'repressor'
    print('number of calls for regions tested in both orientations')
    print(both.groupby(['call']).size())
    
    # Save
    both = both[['seqnames', 'start', 'end', 'name', 'activity', 'strand', 
             'activity_z_score', 'call', 
             'logFC_for', 'P.Value_for', 'q_for', 'z_score_for', 'neg_ctrl_for', 'call_for',  
             'logFC_rev', 'P.Value_rev', 'q_rev', 'z_score_rev', 'neg_ctrl_rev', 'call_rev']]
    both.to_csv(args.out_path+'bin_level/all_results_tested_in_both_orientations.bed.gz', sep='\t', index=False, header=False, compression='gzip')
    print('finished saving all bin-level results for bins tested in both orientations')

    all_enhancers = both[both['call'] == 'enhancer']
    all_repressors = both[both['call'] == 'repressor']

    for data, fname in zip([all_enhancers, all_repressors], 
                           ['enhancer_peak_in_both_orientation.bed.gz', 'repressor_peak_in_both_orientation.bed.gz']):
        
        if(len(data) == 0):
            continue
        # Merge adjacent/overlapping bins/elements into peak regions
        merged = merge_bins(data)
        merged = get_activity_after_merging(merged, data, args.mode)

        # Add required columns to adjust format
        merged['strand'] = '.'
        merged['size'] = merged[2].values-merged[1].values

        merged.columns = ['seqnames', 'start', 'end', 'mean_logFC', 'z_score_for_mean_logFC', 'logFC_list', 
                          'summit_seqnames', 'summit_start', 'summit_end', 'summit_logFC', 'summit_z_score', 'strand', 'size']
        merged = merged.sort_values(['seqnames', 'start'])

        # Add peak names and re-arrange columns
        merged['name'] = ['peak_both_{0}'.format(i) for i in range(1, len(merged)+1)]
        merged = merged[['seqnames', 'start', 'end', 'name', 'mean_logFC', 'strand', 
                         'z_score_for_mean_logFC', 'logFC_list', 
                         'summit_seqnames', 'summit_start', 'summit_end', 'summit_logFC', 'summit_z_score', 'size']]
        
        print('Number of peaks in '+fname+' '+str(len(merged)))

        # Write merged peak file
        merged.to_csv(args.out_path+'peak_level/'+fname, sep='\t', index=False, header=False, compression='gzip')
        print('finished saving ' + fname)
    
# Cleanup the temp directory for pybedtools
shutil.rmtree(args.out_path+'tmp_pybedtools/')
