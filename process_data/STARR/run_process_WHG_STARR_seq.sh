#!/bin/bash
#SBATCH --job-name=WHG_STARR_BamToMatrix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=48:00:00
#SBATCH --output=/fs/cbsuhy01/storage/jz855/logs/whg_starr_bam2matrix_%j.out
#SBATCH --error=/fs/cbsuhy01/storage/jz855/logs/whg_starr_bam2matrix_%j.err
#SBATCH --chdir=/fs/cbsuhy01/storage/jz855/STARR_seq_code/public_analysis_code/process_data/STARR/
#SBATCH --mail-user=jz855@cornell.edu
#SBATCH --mail-type=ALL
set -o errexit
set -o pipefail
set -o nounset

echo "Job started at: $(date)"
# before conda activate
set +u
source /home/jz855/miniconda3/etc/profile.d/conda.sh
conda activate py39
set -u

hy1_path=/fs/cbsuhy01/storage/jz855/

# --------------------------------------------------
# Parameters
# --------------------------------------------------
DATA_DIR=${hy1_path}/STARR_seq_dataset/WHG_STARR_seq_TR/bamFiles/    # contains DNA/ and RNA/ BAMs
OUTPUT_DIR=${hy1_path}/STARR_seq_dataset/WHG_STARR_seq_TR/processed_out/
CHROM_SIZES=${hy1_path}/STARR_seq_code/Final_Code_Sharing/data/reference/hg38/hg38.chrom.sizes.txt

MIN_FRAG_SIZE=200
MAX_FRAG_SIZE=800
BIN_SIZE=100
STEP_SIZE=10
COUNT_TYPE=bin
UMI_USAGE=None
N_CPU=8

# Per-sample count filename produced in bin mode:
COUNT_FILENAME="genomic_bin_count_binSize_${BIN_SIZE}_stepSize_${STEP_SIZE}_full_overlap.sorted.bed"

# --------------------------------------------------
# Stage 1: BAM -> per-sample genomic bin counts
# --------------------------------------------------
echo "[1/2] Running STARR_BamToCount.py ..."
python STARR_BamToCount.py \
  --data_dir "$DATA_DIR" \
  --output_dir "$OUTPUT_DIR" \
  --chrom_sizes "$CHROM_SIZES" \
  --min_frag_size $MIN_FRAG_SIZE \
  --max_frag_size $MAX_FRAG_SIZE \
  --bin_size $BIN_SIZE \
  --step_size $STEP_SIZE \
  --count_type $COUNT_TYPE \
  --umi_usage $UMI_USAGE \
  --n_cpu $N_CPU

echo "Stage 1 completed at: $(date)"

# --------------------------------------------------
# Stage 2: Merge per-sample counts -> final matrix
# --------------------------------------------------
DNA_OUT_DIR="${OUTPUT_DIR}/DNA"
RNA_OUT_DIR="${OUTPUT_DIR}/RNA"
MERGE_OUT="${OUTPUT_DIR}/final_matrix"
mkdir -p "$MERGE_OUT"

echo "[2/2] Running STARR_CountToMatrix.py ..."
python STARR_CountToMatrix.py \
  --dna_dir "$DNA_OUT_DIR" \
  --rna_dir "$RNA_OUT_DIR" \
  --count_filename "$COUNT_FILENAME" \
  --output_dir "$MERGE_OUT" \
  --count_type "$COUNT_TYPE"

echo "Final matrix is in: $MERGE_OUT"
echo "Job ended at: $(date)"