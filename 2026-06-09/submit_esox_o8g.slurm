#!/bin/bash
###############################################################################
# SLURM array job: run esox O8G pipeline, one array task per sample.
#
# Submit from the directory that contains run_esox_o8g.sh:
#     mkdir -p logs
#     sbatch submit_esox_o8g.slurm
#
# 5 samples -> array indices 0..4. "%2" = at most 2 running at once.
###############################################################################
#SBATCH -p gpu100
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -J esox_o8g
#SBATCH --output logs/esox_o8g_%A_%a.out
#SBATCH --error  logs/esox_o8g_%A_%a.err
#SBATCH --gres=gpu:1
#SBATCH --mem=200G
#SBATCH --cpus-per-task=30
#SBATCH --array=0-4%2

set -euo pipefail

# --- activate conda env -----------------------------------------------------
# pick whichever works on your cluster:
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate esox_env
# or: source ~/.bashrc && conda activate esox_env

# --- map array index -> sample ----------------------------------------------
SAMPLES=(OE1 RDM_-6bp_SAM_DNA RDM_-8bp_SAM_DNA RDM_OE2_SAM_DNA WT_SAM_DNA)
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

# directory where run_esox_o8g.sh lives (= dir you ran sbatch from)
RUN_SCRIPT="${SLURM_SUBMIT_DIR}/run_esox_o8g.sh"

echo "==============================================================="
echo "Array task : ${SLURM_ARRAY_TASK_ID}"
echo "Sample     : ${SAMPLE}"
echo "Node       : $(hostname)"
echo "GPU(s)     : ${CUDA_VISIBLE_DEVICES:-unset}"
echo "Start      : $(date)"
echo "==============================================================="

bash "${RUN_SCRIPT}" "${SAMPLE}"

echo "Finished ${SAMPLE} at $(date)"
