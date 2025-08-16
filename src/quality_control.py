# === Standard Library ===
import os
import glob
import argparse
from subprocess import call, PIPE, run, Popen, STDOUT
from multiprocessing import Pool, cpu_count
from random import sample

# === Data Handling ===
import pandas as pd
import numpy as np
import scipy
import statsmodels.stats.multitest as smm

# === Plotting ===
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm
from matplotlib.colors import Normalize
import matplotlib.colors as clr
import matplotlib.patches as mpatches
import matplotlib.ticker as mtick
from matplotlib.offsetbox import AnchoredText
from mpl_toolkits.axes_grid1 import make_axes_locatable

# === Plotting Configuration ===
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Helvetica Neue'

# === Bioinformatics Tools ===
import pysam
import pybedtools
pybedtools.helpers.set_tempdir('./')

# === Project Utilities ===
from utils import *


class CheckRepCorr():
    
    """
    Class to calculate and output correlation between replicates
    for CPM, log-CPM, and log2Ratio-transformed data.
    """
    
    def __init__(self, num_rep_dna, num_rep_rna):
        self.num_rep_dna = num_rep_dna
        self.num_rep_rna = num_rep_rna
        self.sample_list = ['DNA{0}'.format(i) for i in range(1, self.num_rep_dna+1)]+['RNA{0}'.format(i) for i in range(1, self.num_rep_rna+1)]
        self.min_num_rep = np.min([self.num_rep_dna, self.num_rep_rna])
    
    def calc_cpm(self, data, prior_count = 0.5, log=True):
        
        """
        Calculate counts per million (CPM), optionally log2-transformed.

        Args:
            data (pd.DataFrame): Raw count matrix
            prior_count (float): Prior added before log transform
            log (bool): Whether to return log2 CPM

        Returns:
            pd.DataFrame: CPM-transformed matrix
        """
        
        # Total read count per sample (library size) 
        col_sum = data.sum(axis=0)
        if(log):
            # Average library size across samples
            avg_lib_size = np.mean(col_sum)
            
            # Adjust prior based on each sample's library size
            adj_prior = prior_count*col_sum/avg_lib_size
            
            # Adjusted library size accounts for added prior
            adj_lib_size = col_sum + 2*adj_prior
            
            # Compute log2 CPM (in base 2)
            logcpm = (np.log2(data+adj_prior) + np.log2(1e6) - np.log2(adj_lib_size))/(np.log2(2))
            return(logcpm)
        else:
            # Return CPM if log=False
            return(data.mul(1e6).div(col_sum))
                   
    def calc_log2Ratio(self, log2cpm):
        """
        Calculate log2 ratio of RNA to DNA from log2 CPM.

        Args:
            log2cpm (pd.DataFrame): Log2 CPM matrix with RNA and DNA columns

        Returns:
            pd.DataFrame: log2(RNA / DNA) values
        """
        # Slice DNA samples from sample list
        if(self.num_rep_dna == 1):
            dna = log2cpm[self.sample_list[:self.num_rep_dna]]
            dna.columns = ['Rep1'] # Standardize single column name
        else:
            # Use up to the minimum number of replicates to match RNA
            dna = log2cpm[self.sample_list[:self.min_num_rep]]
            dna.columns = ['Rep{0}'.format(x) for x in range(1, self.min_num_rep+1)]
        
        # Slice RNA samples from sample list
        if(self.num_rep_dna > self.num_rep_rna):
            # If more DNA replicates than RN, only match the minimum
            rna = log2cpm[self.sample_list[ self.num_rep_dna:(self.num_rep_dna+self.min_num_rep) ]]
            rna.columns = ['Rep{0}'.format(x) for x in range(1, self.min_num_rep+1)]
        else:
            # Else use all RNA replicates
            rna = log2cpm[self.sample_list[ self.num_rep_dna:(self.num_rep_dna+self.num_rep_rna) ]]
            rna.columns = ['Rep{0}'.format(x) for x in range(1, self.num_rep_rna+1)]
            
        # Compute log2 ratio of RNA over DNA
        if(self.num_rep_dna > 1):
            # Element-wise substrat DNA from RNA (Rep1, Rep2,...)
            log2Ratio = rna-dna
        else:
            # Substract DNA Rep1 from each RNA column
            log2Ratio = rna[['Rep{0}'.format(x) for x in range(1, self.num_rep_rna+1)]].subtract(dna.Rep1, axis=0)
        
        return log2Ratio

    def compute_pairwise_corr(self, df, col_list):
        """
        Compute Spearman and Pearson correlations between all column pairs.

        Args:
            df (pd.DataFrame): Input matrix
            col_list (List[str]): Columns to correlate

        Returns:
            pd.DataFrame: Pairwise correlation results
        """
        # Initialize result lists
        spearmanr_list, pearsonr_list, pair_1_list, pair_2_list = [], [], [], []

        # Iterate over all unique column pairs
        for i in range(len(col_list)):
            for j in range(i + 1, len(col_list)):
                col_1, col_2 = col_list[i], col_list[j]
                # Compute Spearman and Pearson correlations between the two columns and append results to lists
                spearmanr_list.append(scipy.stats.spearmanr(df[col_1], df[col_2])[0])
                pearsonr_list.append(scipy.stats.pearsonr(df[col_1], df[col_2])[0])
                # Append results to lists
                pair_1_list.append(col_1)
                pair_2_list.append(col_2)

        # Combine into result DataFrame
        return pd.DataFrame({
            'Sample_1': pair_1_list,
            'Sample_2': pair_2_list,
            'Spearman_r': spearmanr_list,
            'Pearson_r': pearsonr_list
        })

    def compute_cpm_correlation(self, data, DNA_cols, RNA_cols, output_path, dataset_name):
        """
        Write CPM correlation file for DNA and RNA.

        Args:
            data (pd.DataFrame): Raw count matrix
            DNA_cols (List[str])
            RNA_cols (List[str])
            output_path (str)
            dataset_name (str)
        """
        # Compute CPM
        cpm = self.calc_cpm(data, log=False)

        # Compute pairwise correlation between DNA replicates
        df_DNA = self.compute_pairwise_corr(cpm[DNA_cols], DNA_cols)
        df_DNA['count_type'] = 'cpm'
        
        # Compute pairwise correlation between RNA replicates
        df_RNA = self.compute_pairwise_corr(cpm[RNA_cols], RNA_cols)
        df_RNA['count_type'] = 'cpm'
        
        # Combine DNA and RNA correlation results
        df = pd.concat([df_DNA, df_RNA], ignore_index=True)
        df['Dataset'] = dataset_name
    
        # Save to file
        df.to_csv(os.path.join(output_path, 'cpm_count_corr.txt'), sep='\t', index=False)

    def compute_log_cpm_correlation(self, data, DNA_cols, RNA_cols, output_path, dataset_name):
        """
        Write log-CPM correlation file for DNA and RNA.

        Args:
            data (pd.DataFrame): Raw count matrix
            DNA_cols (List[str])
            RNA_cols (List[str])
            output_path (str)
            dataset_name (str)
        """
        # Calculate log_CPM
        log_cpm = self.calc_cpm(data, log=True)
                   
        # Compute pairwise correlation between DNA replicates
        df_DNA = self.compute_pairwise_corr(log_cpm[DNA_cols], DNA_cols)
        df_DNA['count_type'] = 'log_cpm'
        
        # Compute pairwise correlation between RNA replicates
        df_RNA = self.compute_pairwise_corr(log_cpm[RNA_cols], RNA_cols)
        df_RNA['count_type'] = 'log_cpm'

        # Combine DNA and RNA correlation results
        df = pd.concat([df_DNA, df_RNA], ignore_index=True)
        df['Dataset'] = dataset_name

        # Save to file
        df.to_csv(os.path.join(output_path, 'log_cpm_count_corr.txt'), sep='\t', index=False)

    def compute_log2ratio_correlation(self, data, output_path):
        """
        Write log2Ratio correlation file across all replicates.

        Args:
            data (pd.DataFrame): log2Ratio matrix
            output_path (str)
        """
        # Get list of all replicate column names
        rep_list = data.columns.tolist()
        
        # Compute pairwise correlation between all replicates
        df = self.compute_pairwise_corr(data, rep_list)
        
        # Save to file
        df.to_csv(os.path.join(output_path, 'log2ratio_corr.txt'), sep='\t', index=False)         
                   

class CheckCoverage():
    """
    Class to compute genome-wide or region-specific coverage statistics
    from BED files using pybedtools.
    """
    
    def calc_genome_wide_coverage(self, query_file):
        """
        Calculate total number of base pairs covered genome-wide by query BED file
        
        Args:
            query_file (str):Path to a BED file
        
        Returns:
            tuple: (filename, total base-pair coverage)
        """
        query = pybedtools.BedTool(query_file)
        return((query_file, query.total_coverage()))
        
    def calc_bp_intersect_region(self, args):
                   
        """
        Compute coverage statistics of a query BED file over a target region BED.

        Args:
            args (tuple): (query_file_path, pybedtools.BedTool of region file)

        Returns:
            tuple: (
                filename,
                total bp covered within region,
                number of fragments overlapping region,
                number of region intervals with at least one overlap,
                total bp in region
            )
        """
        
        file, region_bed = args

        bed = pybedtools.BedTool(file)
                   
        # Intersect: get fragments that overlap region BED
        intersect_bed = bed.intersect(region_bed, wa=True)
        intersect = intersect_bed.to_dataframe(disable_auto_names=True, header=None)
        intersect = intersect.drop_duplicates()

        # Compute coverage of query BED on region BED
        query_cover_on_region_bed = region_bed.coverage(intersect_bed)
        query_cover_on_region_df = query_cover_on_region_bed.to_dataframe(disable_auto_names=True, header=None)

        # Extract statistics
        num_bp = np.sum(query_cover_on_region_df[query_cover_on_region_df.columns.tolist()[-3]].values)
        num_intersect_query = len(intersect)
        num_intersect_region = len(query_cover_on_region_df[query_cover_on_region_df[query_cover_on_region_df.columns.tolist()[-1]]>0])
        total_bp = np.sum(query_cover_on_region_df[query_cover_on_region_df.columns.tolist()[-2]].values)

        fname = file.split('/')[-1]

        return((fname, num_bp, num_intersect_query, num_intersect_region, total_bp))


    def calc_region_coverage(self, region_file, query_file_list):
        
        """
        Compute coverage statistics across multiple query files for a single region file.

        Args:
            region_file (str): Path to region BED file
            query_file_list (List[str]): List of query BED file paths

        Returns:
            pd.DataFrame: Coverage stats per query file, including:
                - query filename
                - number of bp covered in region
                - number of overlapping fragments
                - number of regions with overlaps
                - total region bp
                - region file path
        """
        
        region_bed = pybedtools.BedTool(region_file)

        # obtain num_bp
        arg_list = [(query_file, region_bed) for query_file in query_file_list]
        with Pool(5) as pool:
            result_list = pool.map(self.calc_bp_intersect_region, arg_list)

        result_df = pd.DataFrame(result_list, columns = ['query_file','num_bp', 'num_frag', 'num_region', 'size'])
        result_df['region'] = region_file

        return(result_df)
    
###########################################################################################    

class CheckLibRecoveryRate():
    """
    Class to compute basic library recovery statistics from a count matrix.
    Includes:
      - total count per library
      - number of fragments with detectable counts
      - percent of library recovered
    """
    
    def __init__(self, num_rep_dna, num_rep_rna):
                   
        """
        Initialize with the number of DNA and RNA replicates.

        Args:
            num_rep_dna (int): Number of DNA replicates
            num_rep_rna (int): Number of RNA replicates
        """
        self.num_rep_dna = num_rep_dna
        self.num_rep_rna = num_rep_rna
        
    def get_stats_from_count_mat(self, path_to_count_mat, count_cutoff=0):
        """
        Extract total counts and the number of fragments sequenced per library from a large count matrix,
        processing the file in chunks.

        If the file ends with '.txt', it is assumed to have a header with proper column names (e.g., DNA1, DNA2, ..., RNA1, RNA2, ...),
        so the function simply subsets these replicate columns.

        If the file ends with '.bed', it is assumed to have no header and that the last (num_rep_dna + num_rep_rna) columns correspond 
        to replicate counts. In that case, the function subsets those columns and assigns new column names (DNA1, ..., RNA1, ...).

        Args:
            path_to_count_mat (str): Path to the count matrix file.
            count_cutoff (int): Threshold count; a fragment is considered "sequenced" if its count exceeds this value.

        Returns:
            tuple:
                - tot_count_total (np.array): Total counts per library (for DNA and RNA replicates).
                - num_sequenced_total (np.array): Number of fragments with counts > count_cutoff per library.
                - num_frag (int): Total number of fragments (rows) processed.
        """
        import numpy as np

        # Initialize storage for chunk-wise statistics.
        tot_count_list = []
        num_sequenced_list = []
        num_frag = 0
        chunksize = 10**6  # Use a chunk size of 1,000,000 rows.

        # Determine the expected replicate column names.
        expected_replicate_cols = (['DNA{}'.format(i) for i in range(1, self.num_rep_dna + 1)] +
                                   ['RNA{}'.format(i) for i in range(1, self.num_rep_rna + 1)])
        n_replicates = len(expected_replicate_cols)

        # Process file differently depending on its extension.
        if path_to_count_mat.endswith('.txt'):
            # File has a header with column names.
            for chunk in pd.read_csv(path_to_count_mat, sep='\t', chunksize=chunksize, header=0):
                # Subset the columns using the expected replicate column names.
                chunk = chunk[expected_replicate_cols]
                tot_count = []
                num_sequenced = []
                # Loop through each replicate column.
                for col in expected_replicate_cols:
                    tot_count.append(chunk[col].sum())
                    num_sequenced.append((chunk[col] > count_cutoff).sum())
                tot_count_list.append(tot_count)
                num_sequenced_list.append(num_sequenced)
                num_frag += len(chunk)
        elif path_to_count_mat.endswith('.bed'):
            # File does not have a header.
            for chunk in pd.read_csv(path_to_count_mat, sep='\t', chunksize=chunksize, header=None):
                # Subset the last n_replicates columns, which correspond to the replicate counts.
                replicate_columns = chunk.columns[-n_replicates:]
                chunk = chunk[replicate_columns]
                # Assign new column names: DNA replicates first, then RNA replicates.
                chunk.columns = expected_replicate_cols
                tot_count = []
                num_sequenced = []
                # Loop through each replicate column.
                for col in expected_replicate_cols:
                    tot_count.append(chunk[col].sum())
                    num_sequenced.append((chunk[col] > count_cutoff).sum())
                tot_count_list.append(tot_count)
                num_sequenced_list.append(num_sequenced)
                num_frag += len(chunk)
        else:
            raise ValueError("Unsupported file extension. Expected a '.txt' or '.bed' file.")

        # Combine chunk-wise statistics into overall totals.
        tot_count_array = np.array(tot_count_list)
        tot_count_total = tot_count_array.sum(axis=0)

        num_sequenced_array = np.array(num_sequenced_list)
        num_sequenced_total = num_sequenced_array.sum(axis=0)

        return tot_count_total, num_sequenced_total, num_frag
    
    def calc_sequencing_coverage(self, tot_num_frag, num_sequenced_list):
        
        """
        Calculate percentage of fragments recovered (sequenced) per library.

        Args:
            tot_num_frag (int): Total number of fragments in the library
            num_sequenced_list (List[int]): Number of sequenced fragments per library

        Returns:
            pd.DataFrame: Coverage table with columns:
                - library
                - total (fragments)
                - num_sequenced
                - pct (percent sequenced)
        """
        
        df = pd.DataFrame({'library':['DNA{0}'.format(i) for i in range(1, self.num_rep_dna+1)]+['RNA{0}'.format(i) for i in range(1, self.num_rep_rna+1)], 
                   'total':[tot_num_frag]*(self.num_rep_dna+self.num_rep_rna), 
                   'num_sequenced':list(num_sequenced_list)})
        df['pct'] = df['num_sequenced'].div(df['total'])*100
        
        return(df)
    
    
    

    
    