#!/bin/bash
#SBATCH -p com300
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=2
#SBATCH -J extract_ratio
#SBATCH --output log/extract_ratio_%A_%a.out
#SBATCH --error log/extract_ratio_%A_%a.err
#SBATCH --array=0-111%20

# ============================================================
# 后续提取脚本：pileup.bed 已生成，跳过 modkit pileup。
# 优化点：单次扫描 pileup.bed，一次性拆分到各 code 的输出文件，
#         避免对 2.4 亿行大文件重复读取 8 次，显著节约算力/IO。
# ============================================================

source "$HOME/anaconda3/etc/profile.d/conda.sh"
conda activate Dorado

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"
mkdir -p log

# ============ 路径配置 ============
WORKDIR="/public/home/h14166/test/All/mapping"
LIST="/public/home/h14166/test/All/list"
MIN_COV=1

# ============ 读取样本名 ============
NLINES=$(wc -l < "$LIST")
TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"

if [ "$TASK_ID" -ge "$NLINES" ]; then
  echo "Task ${TASK_ID} exceeds list lines ${NLINES}, exit."
  exit 0
fi

sample="$(sed -n "$((TASK_ID+1))p" "$LIST" | tr -d '\r')"
[ -z "$sample" ] && { echo "Empty sample name"; exit 1; }

SAMPLE_DIR="${WORKDIR}/${sample}"

echo "=========================================="
echo "[$(date)] 样本: ${sample}"
echo "=========================================="

# ============ 检查 pileup.bed ============
if [ ! -s "${SAMPLE_DIR}/pileup.bed" ]; then
  echo "pileup.bed 不存在或为空: ${SAMPLE_DIR}/pileup.bed，跳过。"
  exit 0
fi

cd "$SAMPLE_DIR"

# ============ 单次扫描，一次性拆分各修饰比例 ============
echo "[$(date)] 提取修饰比例 (single-pass) ..."
awk -v min_cov="$MIN_COV" '
  BEGIN{
    n = split("m a 17596 17802 19227 19228 19229 69426", arr, " ")
    for (i = 1; i <= n; i++) {
      code = arr[i]
      want[code] = 1
      out = code "_ratio.tsv"
      # 写表头（首次写会创建/清空文件，后续数据行追加写入）
      print "#chrom\tstart\tend\tfraction\tcov" > out
    }
  }
  ($4 in want) && ($10 + 0 >= min_cov) {
    p = $11 + 0
    if (p > 1) p /= 100
    if (p >= 0) {
      out = $4 "_ratio.tsv"
      printf "%s\t%d\t%d\t%.4f\tcov=%d\n", $1, $2, $3, p, $10 > out
    }
  }
  END{
    for (c in want) close(c "_ratio.tsv")
  }
' pileup.bed

# ============ 统计各文件位点数 ============
echo "[$(date)] 提取完成，各 code 位点数（已扣除 1 行表头）:"
for code in m a 17596 17802 19227 19228 19229 69426; do
  f="${code}_ratio.tsv"
  if [ -f "$f" ]; then
    n=$(($(wc -l < "$f") - 1))
    echo "  ${f}: ${n} sites"
  fi
done

echo "=========================================="
echo "[$(date)] 样本 ${sample} 完成"
echo "=========================================="
