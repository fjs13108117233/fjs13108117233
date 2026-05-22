#!/bin/bash
#SBATCH -p com300
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -J sRNA_diff
#SBATCH --output logs/sRNA_diff_%j.out
#SBATCH --error  logs/sRNA_diff_%j.err

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"
mkdir -p logs

source "$HOME/anaconda3/etc/profile.d/conda.sh"
conda activate RNA

python3 sRNA_diff_pipeline.py

echo "=== 完成: $(date) ==="
