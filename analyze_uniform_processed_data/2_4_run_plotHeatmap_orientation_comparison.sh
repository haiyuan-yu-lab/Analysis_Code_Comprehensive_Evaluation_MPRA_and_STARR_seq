#!/bin/bash
# run_plotHeatmap_all_orientation_comparison.sh
#
# Description:
#   This script generates heatmaps from computeMatrix outputs comparing 
#   orientation-independent peaks versus either orientation peaks for three datasets:
#     - LentiMPRA
#     - ATAC_STARR_seq
#     - WHG_STARR_seq
#
#   For each dataset, two heatmaps are generated:
#     1. DNase-seq & ATAC-seq signal profiles, using the GnBu colormap.
#     2. Histone modification profiles (H3K4me3 & H3K27ac), using the BuPu colormap.
#
#   The computeMatrix parameters used are as follows:
#       --referencePoint center --beforeRegionStartLength 1000 --afterRegionStartLength 1000
#       --sortRegions descend --sortUsing mean
#
#   Input matrix files are expected in:
#       ${BASE_DIR}/output/computeMatrix/orientation_comparison/matrix_DNase_ATAC
#       ${BASE_DIR}/output/computeMatrix/orientation_comparison/matrix_Histone
#
#   Heatmap PDFs will be saved in:
#       ${BASE_DIR}/plot/heatmap/DNase_ATAC/orientation_comparison
#       ${BASE_DIR}/plot/heatmap/Histone/orientation_comparison
#
# ----------------------------------------------------------------------

# Define the base directory for Final_Code_Sharing.
BASE_DIR="/fs/cbsuhy01/storage/jz855/STARR_seq_code/Final_Code_Sharing"

# -------------------------------
# Define input matrix directories for orientation comparisons.
# -------------------------------
MATRIX_DNASE_ATAC_ORI_DIR="${BASE_DIR}/data/output/computeMatrix/orientation_comparison/matrix_DNase_ATAC"
MATRIX_HISTONE_ORI_DIR="${BASE_DIR}/data/output/computeMatrix/orientation_comparison/matrix_Histone"

# -------------------------------
# Define output directories for heatmaps (orientation comparisons).
# -------------------------------
# OUT_HEATMAP_DNASE_ATAC_ORI="${BASE_DIR}/plot/heatmap/DNase_ATAC/orientation_comparison"
# OUT_HEATMAP_HISTONE_ORI="${BASE_DIR}/plot/heatmap/Histone/orientation_comparison"

OUT_HEATMAP_DNASE_ATAC_ORI="${BASE_DIR}/plot/heatmap_with_legend/DNase_ATAC/orientation_comparison"
OUT_HEATMAP_HISTONE_ORI="${BASE_DIR}/plot/heatmap_with_legend/Histone/orientation_comparison"

mkdir -p "$OUT_HEATMAP_DNASE_ATAC_ORI"
mkdir -p "$OUT_HEATMAP_HISTONE_ORI"

# -------------------------------
# Loop over each dataset and run plotHeatmap for both signal profiles.
# -------------------------------
for dataset in "LentiMPRA" "ATAC_STARR" "WHG_STARR"; do
    # Convert the dataset name to lowercase for consistent output file naming.
    base=$(echo "$dataset" | tr '[:upper:]' '[:lower:]')
    
    echo "Generating orientation comparison heatmaps for dataset: $dataset"
    
    # ----------------------------------------------------------------
    # 1. DNase-seq & ATAC-seq heatmap (using GnBu colormap).
    # ----------------------------------------------------------------
    MATRIX_FILE_DA="${MATRIX_DNASE_ATAC_ORI_DIR}/${base}.gz"
    OUTFILE_DA="${OUT_HEATMAP_DNASE_ATAC_ORI}/${base}.pdf"
    
    echo "  Plotting DNase-seq & ATAC-seq heatmap for $dataset using matrix file: ${MATRIX_FILE_DA}"
    plotHeatmap -m "$MATRIX_FILE_DA" \
        -o "$OUTFILE_DA" \
        --dpi 300 --colorMap GnBu --alpha 0.8 \
        --yMin 0 --yMax 9 --zMin 0 --zMax 9 \
        --heatmapHeight 10 --heatmapWidth 5 \
        --legendLocation best \
        --refPointLabel Center --samplesLabel "DNase-seq" "ATAC-seq" \
        --regionsLabel "Both" "Either" &
    
    # ----------------------------------------------------------------
    # 2. Histone modifications heatmap (using BuPu colormap).
    # ----------------------------------------------------------------
    MATRIX_FILE_H="${MATRIX_HISTONE_ORI_DIR}/${base}.gz"
    OUTFILE_H="${OUT_HEATMAP_HISTONE_ORI}/${base}.pdf"
    
    echo "  Plotting Histone modifications heatmap for $dataset using matrix file: ${MATRIX_FILE_H}"
    plotHeatmap -m "$MATRIX_FILE_H" \
        -o "$OUTFILE_H" \
        --dpi 300 --colorMap BuPu --alpha 0.8 \
        --yMin 0 --yMax 20 --zMin 0 --zMax 5 \
        --heatmapHeight 10 --heatmapWidth 5 \
        --legendLocation best \
        --refPointLabel Center --samplesLabel "H3K4me3" "H3K27ac" \
        --regionsLabel "Both" "Either" &
done

# Wait for all background plotHeatmap jobs to complete.
wait

echo "All orientation comparison heatmaps have been generated."