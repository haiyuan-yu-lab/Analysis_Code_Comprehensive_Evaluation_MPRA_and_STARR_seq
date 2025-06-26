#!/bin/bash
# run_plotHeatmap_all_original_vs_processed.sh
#
# Description:
#   This script generates heatmaps from computeMatrix output comparing lab-reported (original)
#   peaks versus uniformly processed peaks for several datasets. For each dataset, the script
#   runs plotHeatmap twice:
#     1. To compare DNase-seq and ATAC-seq profiles using the GnBu colormap.
#     2. To compare histone modification profiles (H3K4me3 and H3K27ac) using the BuPu colormap.
#
# The computeMatrix output matrix files are expected to be found in:
#   ${BASE_DIR}/output/computeMatrix/original_vs_processed/matrix_DNase_ATAC
#   ${BASE_DIR}/output/computeMatrix/original_vs_processed/matrix_Histone
#
# Heatmap PDF outputs will be stored in:
#   ${BASE_DIR}/plot/heatmap/DNase_ATAC/original_vs_processed
#   ${BASE_DIR}/plot/heatmap/Histone/original_vs_processed
#
# This script uses the following fixed plotting parameters:
#   --referencePoint center --beforeRegionStartLength 1000 --afterRegionStartLength 1000
#   --sortRegions descend --sortUsing mean
#
# ----------------------------------------------------------------------

# Define the base directory for Final_Code_Sharing.
BASE_DIR="/fs/cbsuhy01/storage/jz855/STARR_seq_code/Final_Code_Sharing"

# -------------------------------
# Define directories for computeMatrix outputs.
# -------------------------------
MATRIX_DNASE_ATAC_DIR="${BASE_DIR}/data/output/computeMatrix/original_vs_processed/matrix_DNase_ATAC"
MATRIX_HISTONE_DIR="${BASE_DIR}/data/output/computeMatrix/original_vs_processed/matrix_Histone"

# -------------------------------
# Define output directories for heatmaps.
# -------------------------------
# OUT_HEATMAP_DNASE_ATAC="${BASE_DIR}/plot/heatmap/DNase_ATAC/original_vs_processed"
# OUT_HEATMAP_HISTONE="${BASE_DIR}/plot/heatmap/Histone/original_vs_processed"

OUT_HEATMAP_DNASE_ATAC="${BASE_DIR}/plot/heatmap_with_legend/DNase_ATAC/original_vs_processed"
OUT_HEATMAP_HISTONE="${BASE_DIR}/plot/heatmap_with_legend/Histone/original_vs_processed"

mkdir -p "$OUT_HEATMAP_DNASE_ATAC"
mkdir -p "$OUT_HEATMAP_HISTONE"

# -------------------------------
# Common computeMatrix parameters (for reference):
#   --referencePoint center
#   --beforeRegionStartLength 1000
#   --afterRegionStartLength 1000
#   --sortRegions descend
#   --sortUsing mean
# -------------------------------

# Loop over each dataset.
for dataset in "LentiMPRA" "ATAC_STARR" "WHG_STARR" "TilingMPRA"; do
    # Convert the dataset name to lowercase for consistent output file names.
    base=$(echo "$dataset" | tr '[:upper:]' '[:lower:]')
    
    echo "Generating heatmaps for dataset: $dataset"
    
    # ----------------------------------------------------------------
    # 1. DNase-seq & ATAC-seq heatmap (using GnBu colormap).
    # ----------------------------------------------------------------
    MATRIX_FILE_DA="${MATRIX_DNASE_ATAC_DIR}/${base}.gz"
    OUTFILE_DA="${OUT_HEATMAP_DNASE_ATAC}/${base}.pdf"
    
    echo "  Plotting DNase-seq & ATAC-seq heatmap for $dataset using matrix file: ${MATRIX_FILE_DA}"
    plotHeatmap -m "$MATRIX_FILE_DA" \
        -o "$OUTFILE_DA" \
        --dpi 300 --colorMap GnBu --alpha 0.8 \
        --yMin 0 --yMax 9 --zMin 0 --zMax 9 \
        --heatmapHeight 10 --heatmapWidth 5 \
        --legendLocation best \
        --refPointLabel Center --samplesLabel "DNase-seq" "ATAC-seq" \
        --regionsLabel "Original" "Processed" &
    
    # ----------------------------------------------------------------
    # 2. Histone modifications heatmap (using BuPu colormap).
    # ----------------------------------------------------------------
    MATRIX_FILE_H="${MATRIX_HISTONE_DIR}/${base}.gz"
    OUTFILE_H="${OUT_HEATMAP_HISTONE}/${base}.pdf"
    
    echo "  Plotting Histone modifications heatmap for $dataset using matrix file: ${MATRIX_FILE_H}"
    plotHeatmap -m "$MATRIX_FILE_H" \
        -o "$OUTFILE_H" \
        --dpi 300 --colorMap BuPu --alpha 0.8 \
        --yMin 0 --yMax 20 --zMin 0 --zMax 5 \
        --heatmapHeight 10 --heatmapWidth 5 \
        --legendLocation best \
        --refPointLabel Center --samplesLabel "H3K4me3" "H3K27ac" \
        --regionsLabel "Original" "Processed" &
done

# Wait for all background plotHeatmap jobs to complete.
wait

echo "All heatmaps have been generated."