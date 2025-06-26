#!/bin/bash
# run_computeMatrix_all_original_vs_processed.sh
#
# Description:
#   This script compares signal profiles between lab-reported (original) peaks and uniformly 
#   processed peaks for several datasets. For each dataset, computeMatrix is run twice:
#     1. To compare DNase-seq and ATAC-seq signals.
#     2. To compare H3K4me3 and H3K27ac signals.
#
#   The same parameters are used as in your original commands:
#     --referencePoint center --beforeRegionStartLength 1000 --afterRegionStartLength 1000
#     --sortRegions descend --sortUsing mean
#
# Directories and file paths are defined based on the Final_Code_Sharing root.
#
# -------------------------------------------------------------------

# Define the base directory for Final_Code_Sharing.
BASE_DIR="/fs/cbsuhy01/storage/jz855/STARR_seq_code/Final_Code_Sharing"

# -------------------------------
# Define region files for each dataset using associative arrays.
# -------------------------------
declare -A original_files
declare -A processed_files

# LentiMPRA:
original_files[LentiMPRA]="${BASE_DIR}/data/lab_reported_data/processed_files/LentiMPRA/active_merged_all.bed"
processed_files[LentiMPRA]="${BASE_DIR}/data/uniform_processed_data/LentiMPRA/merged_peak/merged_enhancer_peak_orientation_independent.bed.gz"

# ATAC-STARR-seq:
original_files[ATAC_STARR]="${BASE_DIR}/data/lab_reported_data/processed_files/ATAC_STARR_seq/ATAC_STARR_merged_enhancer_peaks_in_either_orientation.bed"
processed_files[ATAC_STARR]="${BASE_DIR}/data/uniform_processed_data/ATAC_STARR_seq/merged_peak/merged_enhancer_peak_orientation_independent.bed.gz"

# WHG-STARR-seq:
original_files[WHG_STARR]="${BASE_DIR}/data/lab_reported_data/original_files/WHG_STARR_seq/ENCFF908UFR.bed.gz"
processed_files[WHG_STARR]="${BASE_DIR}/data/uniform_processed_data/WHG_STARR_seq/merged_peak/merged_enhancer_peak_orientation_independent.bed.gz"

# TilingMPRA:
original_files[TilingMPRA]="${BASE_DIR}/data/lab_reported_data/processed_files/TilingMPRA/all_TilingMPRA_active_from_lab_reported_data.bed.gz"
processed_files[TilingMPRA]="${BASE_DIR}/data/uniform_processed_data/TilingMPRA/all_TilingMPRA_merged_enhancer_peak_in_either_orientation.bed.gz"

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
OUT_DIR_DNASE_ATAC="${BASE_DIR}/data/output/computeMatrix/original_vs_processed/matrix_DNase_ATAC"
OUT_DIR_HISTONE="${BASE_DIR}/data/output/computeMatrix/original_vs_processed/matrix_Histone"

mkdir -p "$OUT_DIR_DNASE_ATAC"
mkdir -p "$OUT_DIR_HISTONE"

REFERENCE_POINT="center"
BEFORE_LENGTH=1000
AFTER_LENGTH=1000
SORT_REGIONS="descend"
SORT_USING="mean"

# -------------------------------
# Loop over each dataset and run computeMatrix for both sets of signals.
# -------------------------------
for dataset in "LentiMPRA" "ATAC_STARR" "WHG_STARR" "TilingMPRA"; do
    # Retrieve the original and processed region file paths.
    orig="${original_files[$dataset]}"
    proc="${processed_files[$dataset]}"
    
    # Create a lowercase base name for naming output files.
    base=$(echo "$dataset" | tr '[:upper:]' '[:lower:]')
    
    echo "Processing dataset: $dataset"
    
    # ----------------------------------------------------------------
    # 1. DNase-seq and ATAC-seq profile comparison.
    # ----------------------------------------------------------------
    OUTPUT_MATRIX_DA="${OUT_DIR_DNASE_ATAC}/${base}.gz"
    OUTPUT_MATRIX_DA_TAB="${OUT_DIR_DNASE_ATAC}/${base}.tab"
    OUTPUT_MATRIX_DA_SORTED="${OUT_DIR_DNASE_ATAC}/${base}.bed"
    
    echo "  Running computeMatrix for DNase-seq/ATAC-seq for $dataset..."
    computeMatrix reference-point \
        -R "$orig" "$proc" \
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
    OUTPUT_MATRIX_H="${OUT_DIR_HISTONE}/${base}.gz"
    OUTPUT_MATRIX_H_TAB="${OUT_DIR_HISTONE}/${base}.tab"
    OUTPUT_MATRIX_H_SORTED="${OUT_DIR_HISTONE}/${base}.bed"
    
    echo "  Running computeMatrix for Histone (H3K4me3/H3K27ac) for $dataset..."
    computeMatrix reference-point \
        -R "$orig" "$proc" \
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

echo "All computeMatrix jobs have completed."