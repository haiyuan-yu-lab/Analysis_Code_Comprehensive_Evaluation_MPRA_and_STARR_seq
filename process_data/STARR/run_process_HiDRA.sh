#!/bin/bash
#SBATCH --job-name=STARR_BamToMatrix
#SBATCH --nodes=1                    # one physical node
#SBATCH --ntasks=1                   # one main process (your Python script)
#SBATCH --cpus-per-task=8            # 8 CPU cores used by your script via multiprocessing
#SBATCH --mem=256G
#SBATCH --time=24:00:00
#SBATCH --output=/fs/cbsuhy01/storage/jz855/logs/hidra_bam2matrix_%j.out
#SBATCH --error=/fs/cbsuhy01/storage/jz855/logs/hidra_bam2matrix_%j.err
#SBATCH --chdir=/fs/cbsuhy01/storage/jz855/STARR_seq_code/public_analysis_code/process_data/STARR/
#SBATCH --mail-user=jz855@cornell.edu
#SBATCH --mail-type=ALL
set -o errexit;

# --------------------------------------------------
# Environment setup
# --------------------------------------------------
echo "Job started at: $(date)"
echo "Setting up conda environment..."

# init conda env
source /home/jz855/miniconda3/etc/profile.d/conda.sh;
conda activate py39;

# env variables
hy1_path=/fs/cbsuhy01/storage/jz855/

# --------------------------------------------------
# Define inputs
# --------------------------------------------------
echo "Defining input variables..."
DATA_DIR=${hy1_path}/STARR_seq_dataset/GM12878/HiDRA/bamFiles  # must contain 'DNA/' and 'RNA/' subdirs
OUTPUT_DIR=${hy1_path}/STARR_seq_dataset/GM12878/HiDRA/processed_out/
CHROM_SIZES=${hy1_path}/STARR_seq_code/Final_Code_Sharing/data/reference/hg38/hg38.chrom.sizes.txt


MIN_FRAG_SIZE=100
MAX_FRAG_SIZE=600
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