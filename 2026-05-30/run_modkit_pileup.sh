#!/bin/bash
#SBATCH -p com300
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=30
#SBATCH -J modkit_pileup
#SBATCH --output log/modkit_pileup_%A_%a.out
#SBATCH --error log/modkit_pileup_%A_%a.err
#SBATCH --array=0-111%20

source "$HOME/anaconda3/etc/profile.d/conda.sh"
conda activate Dorado

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"
mkdir -p log

# ============ 路径配置 ============
WORKDIR="/public/home/h14166/test/All/mapping"
REF="/public/home/h14166/refer/maizev4.fa"
LIST="/public/home/h14166/test/All/list"
THREADS="${SLURM_CPUS_PER_TASK:-30}"

# ============ 读取样本名 ============
NLINES=$(wc -l < "$LIST")
TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"

if [ "$TASK_ID" -ge "$NLINES" ]; then
  echo "Task ${TASK_ID} exceeds list lines ${NLINES}, exit."
  exit 0
fi

sample="$(sed -n "$((TASK_ID+1))p" "$LIST" | tr -d '\r')"
[ -z "$sample" ] && { echo "Empty sample name"; exit 1; }

echo "=========================================="
echo "[$(date)] 样本: ${sample}"
echo "=========================================="

SAMPLE_DIR="${WORKDIR}/${sample}"

if [ ! -f "${SAMPLE_DIR}/aligned.sorted.bam" ]; then
  echo "File ${SAMPLE_DIR}/aligned.sorted.bam does not exist"
  exit 0
fi

# ============ Step 1: modkit pileup ============
pushd "$SAMPLE_DIR" >/dev/null || { echo "Cannot enter ${SAMPLE_DIR}"; exit 1; }

if [ ! -f "pileup.bed" ]; then
  echo "[$(date)] Running modkit pileup ..."
  modkit pileup --threads "$THREADS" aligned.sorted.bam pileup.bed --ref "$REF"
  echo "[$(date)] modkit pileup done."
else
  echo "[$(date)] pileup.bed already exists, skip."
fi

# ============ Step 2: 提取各修饰比例 ============
codes=(m a 17596 17802 19227 19228 19229 69426)
for code in "${codes[@]}"; do
  out="${code}_ratio.tsv"
  awk -v code="$code" -v min_cov=1 '
    $4==code && $10>=min_cov {
      p=$11+0
      if(p>1) p/=100
      if(p>=0) printf "%s\t%d\t%d\t%.4f\tcov=%d\n",$1,$2,$3,p,$10
    }' pileup.bed > "$out"
  sed -i '1i#chrom\tstart\tend\tfraction\tcov' "$out"
done

popd >/dev/null

echo "=========================================="
echo "[$(date)] 样本 ${sample} 完成"
echo "=========================================="
