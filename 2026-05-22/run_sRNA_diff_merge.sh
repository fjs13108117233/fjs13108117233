#!/bin/bash
#SBATCH -p com300
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -J sRNA_merge
#SBATCH --output logs/sRNA_merge_%j.out
#SBATCH --error  logs/sRNA_merge_%j.err

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"
mkdir -p logs

source "$HOME/anaconda3/etc/profile.d/conda.sh"
conda activate RNA

# =========================
# 用户参数 (与 run_sRNA_diff.sh 保持一致)
# =========================
WORKDIR="/public/home/h14166/fang/Heyufei/smrna_seq"
TOOLS="${WORKDIR}/tools"
OUTDIR="${WORKDIR}/sRNA_diff_analysis"
N_TASKS=10
MIN_SAMPLES=2

# =========================
# 合并 + DESeq2
# =========================
echo "=== 合并 chunk 并运行 DESeq2: $(date) ==="

python3 ${TOOLS}/sRNA_diff_pipeline.py \
    --workdir ${WORKDIR} \
    --bed ${WORKDIR}/all_sRNA.bed \
    --list ${WORKDIR}/list \
    --outdir ${OUTDIR} \
    --min-samples ${MIN_SAMPLES} \
    --rscript ${TOOLS}/sRNA_deseq2_analysis.R \
    --merge \
    --n-tasks ${N_TASKS}

echo ""
echo "=== 全部完成: $(date) ==="
echo "输出目录: ${OUTDIR}"
echo ""
echo "结果文件:"
ls -lh ${OUTDIR}/results/*.csv 2>/dev/null || echo "  (无)"
echo ""
echo "图片文件:"
ls -lh ${OUTDIR}/plots/* 2>/dev/null || echo "  (无)"
