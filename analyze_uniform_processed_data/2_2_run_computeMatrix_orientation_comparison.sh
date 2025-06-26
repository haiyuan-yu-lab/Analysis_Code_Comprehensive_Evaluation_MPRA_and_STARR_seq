#!/bin/bash
# run_computeMatrix_all_orientation_comparison.sh
#
# Description:
#   This script computes signal matrices comparing orientation-independent peaks
#   versus either orientation peaks for three datasets: LentiMPRA, ATAC_STARR_seq, and WHG_STARR_seq.
#   For each dataset, computeMatrix is run twice:
#     1. For DNase-seq and ATAC-seq signal profiles.
#     2. For histone modification signal profiles (H3K4me3 and H3K27ac).
#
#   The computeMatrix options used are:
#       --referencePoint center
#       --beforeRegionStartLength 1000
#       --afterRegionStartLength 1000
#       --sortRegions descend
#       --sortUsing mean
#
# Directories and file paths are defined based on the Final_Code_Sharing root.
#
# -------------------------------------------------------------------

# Define the base directory for Final_Code_Sharing.
BASE_DIR="/fs/cbsuhy01/storage/jz855/STARR_seq_code/Final_Code_Sharing"

# -------------------------------
# Define region files for each dataset using associative arrays.
# These arrays hold the region file paths for:
#   - Orientation-independent peaks (ori)
#   - Either orientation peaks (either)
# -------------------------------
declare -A orientation_indep_files
declare -A either_orientation_files

# LentiMPRA:
orientation_indep_files[LentiMPRA]="${BASE_DIR}/data/uniform_processed_data/LentiMPRA/merged_peak/merged_enhancer_peak_orientation_independent.bed.gz"
either_orientation_files[LentiMPRA]="${BASE_DIR}/data/uniform_processed_data/LentiMPRA/merged_peak/merged_enhancer_peak_from_either_in_tested_both.bed.gz"

# ATAC_STARR_seq:
orientation_indep_files[ATAC_STARR]="${BASE_DIR}/data/uniform_processed_data/ATAC_STARR_seq/merged_peak/merged_enhancer_peak_orientation_independent.bed.gz"
either_orientation_files[ATAC_STARR]="${BASE_DIR}/data/uniform_processed_data/ATAC_STARR_seq/merged_peak/merged_enhancer_peak_from_either_in_tested_both.bed.gz"

# WHG_STARR_seq:
orientation_indep_files[WHG_STARR]="${BASE_DIR}/data/uniform_processed_data/WHG_STARR_seq/merged_peak/merged_enhancer_peak_orientation_independent.bed.gz"
either_orientation_files[WHG_STARR]="${BASE_DIR}/data/uniform_processed_data/WHG_STARR_seq/merged_peak/merged_enhancer_peak_from_either_in_tested_both.bed.gz"

# -------------------------------
# Define BigWig signal file paths.
# -------------------------------
DNase_seq="${BASE_DIR}/data/reference/K562_DNase_seq/ENCFF972GVB.bigWig"
ATAC_seq="${BASE_DIR}/data/reference/K562_ATAC_seq/ENCFF102ARJ.bigWig"
H3K4me3="${BASE_DIR}/data/reference/K562_H3K4me3/fc_over_ctrl_ENCFF911JVK.bigWig"
H3K27ac="${BASE_DIR}/data/reference/K562_H3K27ac/fc_over_ctrl_ENCFF381NDD.bigWig"

# -------------------------------
# Define output directories and computeMatrix parameters.
# -------------------------------
OUT_DIR_DNASE_ATAC_ORI="${BASE_DIR}/data/output/computeMatrix/orientation_comparison/matrix_DNase_ATAC"
OUT_DIR_HISTONE_ORI="${BASE_DIR}/data/output/computeMatrix/orientation_comparison/matrix_Histone"

mkdir -p "$OUT_DIR_DNASE_ATAC_ORI"
mkdir -p "$OUT_DIR_HISTONE_ORI"

REFERENCE_POINT="center"
BEFORE_LENGTH=1000
AFTER_LENGTH=1000
SORT_REGIONS="descend"
SORT_USING="mean"

# -------------------------------
# Loop over each dataset and run computeMatrix for both sets of signals.
# -------------------------------
# We restrict the loop here to the three datasets for which orientation files are provided.
for dataset in "LentiMPRA" "ATAC_STARR" "WHG_STARR"; do
    # Retrieve the region file paths for orientation-independent and either orientation peaks.
    ori="${orientation_indep_files[$dataset]}"
    either="${either_orientation_files[$dataset]}"
    
    # Create a lowercase version of the dataset name for output naming.
    base=$(echo "$dataset" | tr '[:upper:]' '[:lower:]')
    
    echo "Processing dataset: $dataset (Orientation-independent vs Either orientation)"
    
    # ----------------------------------------------------------------
    # 1. DNase-seq and ATAC-seq profile comparison.
    # ----------------------------------------------------------------
    OUTPUT_MATRIX_DA="${OUT_DIR_DNASE_ATAC_ORI}/${base}.gz"
    OUTPUT_MATRIX_DA_TAB="${OUT_DIR_DNASE_ATAC_ORI}/${base}.tab"
    OUTPUT_MATRIX_DA_SORTED="${OUT_DIR_DNASE_ATAC_ORI}/${base}.bed"
    
    echo "  Running computeMatrix for DNase-seq/ATAC-seq for $dataset..."
    computeMatrix reference-point \
        -R "$ori" "$either" \
        -S "$DNase_seq" "$ATAC_seq" \
        --referencePoint $REFERENCE_POINT \
        --beforeRegionStartLength $BEFORE_LENGTH \
        --afterRegionStartLength $AFTER_LENGTH \
        --sortRegions $SORT_REGIONS \
        --sortUsing $SORT_USING \
        --outFileName "$OUTPUT_MATRIX_DA" \
        --outFileNameMatrix "$OUTPUT_MATRIX_DA_TAB" \
        --outFileSortedRegions "$OUTPUT_MATRIX_DA_SORTED" &
    
    # ----------------------------------------------------------------
    # 2. Histone modification profile comparison (H3K4me3 and H3K27ac).
    # ----------------------------------------------------------------
    OUTPUT_MATRIX_H="${OUT_DIR_HISTONE_ORI}/${base}.gz"
    OUTPUT_MATRIX_H_TAB="${OUT_DIR_HISTONE_ORI}/${base}.tab"
    OUTPUT_MATRIX_H_SORTED="${OUT_DIR_HISTONE_ORI}/${base}.bed"
    
    echo "  Running computeMatrix for Histone (H3K4me3/H3K27ac) for $dataset..."
    computeMatrix reference-point \
        -R "$ori" "$either" \
        -S "$H3K4me3" "$H3K27ac" \
        --referencePoint $REFERENCE_POINT \
        --beforeRegionStartLength $BEFORE_LENGTH \
        --afterRegionStartLength $AFTER_LENGTH \
        --sortRegions $SORT_REGIONS \
        --sortUsing $SORT_USING \
        --outFileName "$OUTPUT_MATRIX_H" \
        --outFileNameMatrix "$OUTPUT_MATRIX_H_TAB" \
        --outFileSortedRegions "$OUTPUT_MATRIX_H_SORTED" &
done

# Wait for all background computeMatrix jobs to complete.
wait

echo "All computeMatrix jobs for orientation comparison have completed."