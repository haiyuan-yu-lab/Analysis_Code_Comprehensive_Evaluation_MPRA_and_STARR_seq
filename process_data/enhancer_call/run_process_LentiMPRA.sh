#!/bin/bash
#SBATCH --job-name=LentiMPRA_process_limma_out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=120:00:00
#SBATCH --output=/fs/cbsuhy01/storage/jz855/logs/LentiMPRA_process_limma_out_%j.out
#SBATCH --error=/fs/cbsuhy01/storage/jz855/logs/LentiMPRA_process_limma_out_%j.err
#SBATCH --chdir=/fs/cbsuhy01/storage/jz855/STARR_seq_code/public_analysis_code/process_data/enhancer_call/
#SBATCH --mail-user=jz855@cornell.edu
#SBATCH --mail-type=ALL

set -euo pipefail

echo "Job started at: $(date)"

# Conda
set +u
source /home/jz855/miniconda3/etc/profile.d/conda.sh
conda activate py39
set -u

# Paths & params
HY1=/fs/cbsuhy01/storage/jz855
INPUT_DIR="${HY1}/STARR_seq_code/Final_Code_Sharing/data/limma_out/LentiMPRA/"
OUTPUT_DIR="${HY1}/STARR_seq_code/Final_Code_Sharing/data/uniform_processed_data/LentiMPRA/"

NUM_REP_DNA=3
NUM_REP_RNA=3
DATA_TYPE="MPRA"
SCRIPT="process_limma_out_to_make_call.py"

# Checks
[[ -d "${INPUT_DIR}" ]] || { echo "ERROR: INPUT_DIR not found: ${INPUT_DIR}"; exit 1; }
[[ -f "${NEG_CTRL_REF}" ]] || { echo "ERROR: NEG_CTRL_REF not found: ${NEG_CTRL_REF}"; exit 1; }
mkdir -p "${OUTPUT_DIR}"

run_mode () {
  local MODE="$1"
  echo "[$(date '+%F %T')] Running ${SCRIPT} (mode=${MODE})"
  python "${SCRIPT}" \
    -i "${INPUT_DIR}" \
    -o "${OUTPUT_DIR}" \
    -d "${NUM_REP_DNA}" \
    -r "${NUM_REP_RNA}" \
    -m "${MODE}" \
    --data_type "${DATA_TYPE}"
}

# Run: both first, then either
run_mode both
run_mode either

echo "All done. Output directory: ${OUTPUT_DIR}"
echo "Job ended at: $(date)"