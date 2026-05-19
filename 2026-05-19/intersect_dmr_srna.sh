#!/bin/bash
#SBATCH -p com300
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH -J dmr_srna
#SBATCH --output logs/dmr_srna_%j.out
#SBATCH --error  logs/dmr_srna_%j.err

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
mkdir -p logs

source "$HOME/anaconda3/etc/profile.d/conda.sh"
conda activate RNA

# =========================
# 路径
# =========================
SRNA_BED="/public/home/h14166/fang/Heyufei/smrna_seq/all_sRNA.bed"
DMR_DIR="/public/home/h14166/fang/Heyufei/Nanopore/filter_location"
OUTDIR="/public/home/h14166/fang/Heyufei/Nanopore/srna_overlap"

mkdir -p "$OUTDIR"

# 4 种修饰
MODS=("4mc" "5hmc" "5mc" "6ma")

# =========================
# 准备 sRNA BED (含完整信息)
# =========================
echo "=== DMR vs sRNA 交集分析: $(date) ==="
echo "sRNA BED: $SRNA_BED"
echo ""

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

# sRNA: chrom start end name type sample sample_count DicerCall PhaseScore_max
tail -n +2 "$SRNA_BED" | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$10}' \
    | sort -k1,1 -k2,2n > "${TMP}/srna_full.bed"

awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' "${TMP}/srna_full.bed" > "${TMP}/srna_coord.bed"

SRNA_TOTAL=$(wc -l < "${TMP}/srna_full.bed")
echo "sRNA loci 总数: $SRNA_TOTAL"
echo ""

# =========================
# 输出文件
# =========================
OVERLAP_ALL="${OUTDIR}/overlap_all.tsv"
OVERLAP_STATS="${OUTDIR}/overlap_stats.tsv"

# 统计表表头
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "mod_type" "region" "DMR_count" "overlap_DMR" "pct_DMR" \
    "overlap_sRNA" "pct_sRNA" \
    > "$OVERLAP_STATS"

# 大表表头标记
HEADER_WRITTEN="no"

# =========================
# 逐修饰逐文件处理
# =========================
for MOD in "${MODS[@]}"; do
    echo "============================================"
    echo ">>> 修饰: $MOD"

    MOD_DIR="${DMR_DIR}/${MOD}"
    if [ ! -d "$MOD_DIR" ]; then
        echo "  [WARNING] 目录不存在: $MOD_DIR"
        continue
    fi

    for DMR_FILE in "${MOD_DIR}"/*_afterselect_*.tsv; do
        [ -f "$DMR_FILE" ] || continue

        FNAME=$(basename "$DMR_FILE")
        if echo "$FNAME" | grep -q "intergenic"; then
            REGION="intergenic"
        else
            REGION="genic"
        fi

        echo ""
        echo "  --- $MOD / $REGION ---"
        echo "  文件: $FNAME"

        # =============================================
        # 获取 DMR 表头
        # =============================================
        DMR_HEADER=$(head -1 "$DMR_FILE")
        DMR_NCOLS=$(echo "$DMR_HEADER" | awk -F'\t' '{print NF}')
        echo "  DMR 列数: $DMR_NCOLS"

        # 写大表表头 (只写一次)
        if [ "$HEADER_WRITTEN" = "no" ]; then
            printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
                "mod_type" "region" \
                "$DMR_HEADER" \
                "sRNA_name" "sRNA_type" "sRNA_sample" "sRNA_sample_count" \
                "sRNA_DicerCall" "sRNA_PhaseScore_max" \
                "overlap_bp" "overlap_ratio_DMR" "overlap_ratio_sRNA" \
                > "$OVERLAP_ALL"
            HEADER_WRITTEN="yes"
        fi

        # =============================================
        # DMR 全列排序
        # =============================================
        tail -n +2 "$DMR_FILE" | sort -k1,1 -k2,2n > "${TMP}/dmr_full.tsv"

        DMR_COUNT=$(wc -l < "${TMP}/dmr_full.tsv")
        echo "  DMR 位点数: $DMR_COUNT"

        if [ "$DMR_COUNT" -eq 0 ]; then
            echo "  [SKIP] 无 DMR"
            printf "%s\t%s\t%d\t%d\t%.2f\t%d\t%.2f\n" \
                "$MOD" "$REGION" 0 0 0 0 0 >> "$OVERLAP_STATS"
            continue
        fi

        awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' "${TMP}/dmr_full.tsv" > "${TMP}/dmr_coord.bed"

        # --- 统计 ---
        OVERLAP_DMR=$(bedtools intersect \
            -a "${TMP}/dmr_coord.bed" -b "${TMP}/srna_coord.bed" -u | wc -l)
        PCT_DMR=$(awk "BEGIN {printf \"%.2f\", $OVERLAP_DMR/$DMR_COUNT*100}")

        OVERLAP_SRNA=$(bedtools intersect \
            -a "${TMP}/srna_coord.bed" -b "${TMP}/dmr_coord.bed" -u | wc -l)
        PCT_SRNA=$(awk "BEGIN {printf \"%.2f\", $OVERLAP_SRNA/$SRNA_TOTAL*100}")

        echo "  与 sRNA 重叠的 DMR: $OVERLAP_DMR / $DMR_COUNT ($PCT_DMR%)"
        echo "  被 DMR 命中的 sRNA: $OVERLAP_SRNA / $SRNA_TOTAL ($PCT_SRNA%)"

        printf "%s\t%s\t%d\t%d\t%s\t%d\t%s\n" \
            "$MOD" "$REGION" "$DMR_COUNT" "$OVERLAP_DMR" "$PCT_DMR" \
            "$OVERLAP_SRNA" "$PCT_SRNA" \
            >> "$OVERLAP_STATS"

        # =============================================
        # bedtools intersect 保留双侧全列
        # =============================================
        bedtools intersect \
            -a "${TMP}/dmr_full.tsv" \
            -b "${TMP}/srna_full.bed" \
            -wa -wb \
            > "${TMP}/intersect_raw.tsv"

        HIT_COUNT=$(wc -l < "${TMP}/intersect_raw.tsv")
        echo "  交集记录数: $HIT_COUNT"

        if [ "$HIT_COUNT" -eq 0 ]; then
            continue
        fi

        # =============================================
        # 计算 overlap_bp, overlap_ratio_DMR, overlap_ratio_sRNA
        # DMR: col2=start col3=end
        # sRNA: col(nc+2)=start col(nc+3)=end
        # =============================================
        awk -v mod="$MOD" -v reg="$REGION" -v nc="$DMR_NCOLS" '
        BEGIN { OFS="\t" }
        {
            # DMR 全部列
            dmr = ""
            for (i=1; i<=nc; i++) {
                dmr = (dmr == "") ? $i : dmr OFS $i
            }

            # DMR 坐标
            dmr_start = $2 + 0
            dmr_end   = $3 + 0
            dmr_len   = dmr_end - dmr_start

            # sRNA 坐标 (在 DMR 列之后)
            srna_start = $(nc+2) + 0
            srna_end   = $(nc+3) + 0
            srna_len   = srna_end - srna_start

            # sRNA 关键列
            sname   = $(nc+4)
            stype   = $(nc+5)
            ssample = $(nc+6)
            scount  = $(nc+7)
            sdicer  = $(nc+8)
            sphase  = $(nc+9)

            # 计算重叠
            ov_start = (dmr_start > srna_start) ? dmr_start : srna_start
            ov_end   = (dmr_end < srna_end) ? dmr_end : srna_end
            ov_bp    = ov_end - ov_start
            if (ov_bp < 0) ov_bp = 0

            # 重叠比例
            if (dmr_len > 0) {
                ratio_dmr = ov_bp / dmr_len
            } else {
                ratio_dmr = 0
            }
            if (srna_len > 0) {
                ratio_srna = ov_bp / srna_len
            } else {
                ratio_srna = 0
            }

            printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%.4f\t%.4f\n", \
                mod, reg, dmr, sname, stype, ssample, scount, sdicer, sphase, \
                ov_bp, ratio_dmr, ratio_srna
        }
        ' "${TMP}/intersect_raw.tsv" >> "$OVERLAP_ALL"

        # sRNA 类型分布
        echo "  sRNA 类型分布:"
        awk -v nc="$DMR_NCOLS" '{print $(nc+5)}' "${TMP}/intersect_raw.tsv" \
            | sort | uniq -c | sort -rn | head -5 \
            | while read cnt tp; do
                printf "      %-20s %s\n" "$tp" "$cnt"
            done

        echo ""
    done
done

# =========================
# 最终输出
# =========================
echo ""
echo "============================================"
echo ">>> 合并大表: $OVERLAP_ALL"
echo "    总行数: $(wc -l < "$OVERLAP_ALL")"
echo ""
echo ">>> 统计表: $OVERLAP_STATS"
column -t -s "$(printf '\t')" "$OVERLAP_STATS"
echo ""
echo ">>> 输出目录: $OUTDIR/"
ls -lh "$OUTDIR"/*.tsv
echo ""
echo "=== 完成: $(date) ==="
