#!/bin/bash
#SBATCH -p com300
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -J data_process
#SBATCH --output logs/run_data_process_%a.out
#SBATCH --error logs/run_data_process_%a.err
#SBATCH --mem=512G
#SBATCH --cpus-per-task=16
#SBATCH --array=0-3%4

set -eo pipefail
source /public/home/h14166/anaconda3/etc/profile.d/conda.sh
conda activate DL
mkdir -p logs

# ============ 路径配置 ============
WORK_DIR="/public/home/h14166/fang/RNA-modeling/WuRi"
TOOL="${WORK_DIR}/tools/infer_preprocess.py"
REF_FA="/public/home/h14166/fang/RNA-modeling/WuRi/Mapping/reference.fa"

# ============ 读取样本名 ============
mapfile -t SAMPLES < "${WORK_DIR}/list"
Sample=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

echo "=========================================="
echo "[$(date)] 开始处理样本: ${Sample}"
echo "=========================================="

# ============ 输入路径 ============
POD5_DIR="${WORK_DIR}/BasecalingResult/${Sample}"
UNALIGNED_BAM="${WORK_DIR}/BasecalingResult/${Sample}/unmapped.bam"
ALIGNED_BAM="${WORK_DIR}/Mapping/${Sample}/sorted.bam"

# ============ 输出路径 ============
OUTPUT_H5="${WORK_DIR}/${Sample}/infer_tensor.h5"
mkdir -p "${WORK_DIR}/${Sample}"

# ============ 运行预处理 ============
python ${TOOL} \
    --pod5 ${POD5_DIR} \
    --unaligned-bam ${UNALIGNED_BAM} \
    --aligned-bam ${ALIGNED_BAM} \
    --ref-fa ${REF_FA} \
    --output-h5 ${OUTPUT_H5} \
    --workers 16

echo "=========================================="
echo "[$(date)] 样本 ${Sample} 处理完成"
echo "输出文件: ${OUTPUT_H5}"
echo "=========================================="
