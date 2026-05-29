#!/usr/bin/env Rscript
# =============================================================================
# plot_sRNA_analysis.R - sRNA × DMR 联合分析可视化
# =============================================================================
# 生成图表:
#   1. sRNA 类型分布柱状图
#   2. sRNA 染色体分布图
#   3. sRNA 相对基因位置分布 (upstream/downstream/overlap)
#   4. sRNA 到基因距离分布直方图
#   5. 目标 sRNA 表达热图 (RPM, z-score)
#   6. Venn 图: sRNA 基因 vs 4 种 DMR 基因
#   7. 候选基因趋势一致性热图 (sRNA + DMR state)
#   8. sRNA log2FC 趋势图 (候选基因)
#
# 用法:
#   Rscript plot_sRNA_analysis.R \
#       --srna       target_sRNA_annotated.csv \
#       --rpm        sRNA_rpm.csv \
#       --meta       sample_metadata.csv \
#       --intersect  sRNA_DMR_intersect.csv \
#       --consistent sRNA_DMR_intersect_consistent.csv \
#       --outdir     plots
# =============================================================================

suppressMessages({
    library(optparse)
    library(ggplot2)
    library(pheatmap)
    library(RColorBrewer)
})

# ===== 参数 =====
option_list <- list(
    make_option("--srna", type = "character",
                default = "/public/home/h14166/fang/Heyufei/smrna_seq/sRNA_diff_analysis/results/target_sRNA_annotated.csv",
                help = "目标 sRNA 注释结果"),
    make_option("--rpm", type = "character",
                default = "/public/home/h14166/fang/Heyufei/smrna_seq/sRNA_diff_analysis/counts/sRNA_rpm.csv",
                help = "RPM 矩阵"),
    make_option("--meta", type = "character",
                default = "/public/home/h14166/fang/Heyufei/smrna_seq/sRNA_diff_analysis/counts/sample_metadata.csv",
                help = "样本分组"),
    make_option("--intersect", type = "character",
                default = "/public/home/h14166/fang/Heyufei/Unite/sRNA_DMR_intersect.csv",
                help = "sRNA × DMR 全部交集表"),
    make_option("--consistent", type = "character",
                default = "/public/home/h14166/fang/Heyufei/Unite/sRNA_DMR_intersect_consistent.csv",
                help = "趋势一致子集"),
    make_option("--outdir", type = "character",
                default = "plots",
                help = "输出目录")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

MODS <- c("4mc", "5hmc", "5mc", "6ma")

# ===== 读取数据 =====
message(">>> 读取数据 ...")
srna <- read.csv(opt$srna, stringsAsFactors = FALSE)
message("  目标 sRNA: ", nrow(srna))

# =============================================================================
# 图1: sRNA 类型分布
# =============================================================================
message(">>> [1] sRNA 类型分布 ...")
p1 <- ggplot(srna, aes(x = reorder(sRNA_type, sRNA_type, length), fill = sRNA_type)) +
    geom_bar() +
    geom_text(stat = "count", aes(label = after_stat(count)), hjust = -0.2, size = 4) +
    coord_flip() +
    labs(title = "Target sRNA Type Distribution",
         x = "sRNA Type", y = "Count") +
    theme_bw(base_size = 13) +
    theme(legend.position = "none")
ggsave(file.path(opt$outdir, "1_sRNA_type_distribution.pdf"), p1, width = 8, height = 5)
ggsave(file.path(opt$outdir, "1_sRNA_type_distribution.png"), p1, width = 8, height = 5, dpi = 150)

# =============================================================================
# 图2: 染色体分布
# =============================================================================
message(">>> [2] 染色体分布 ...")
srna$chrom <- factor(srna$chrom, levels = sort(unique(as.character(srna$chrom))))
p2 <- ggplot(srna, aes(x = chrom, fill = sRNA_type)) +
    geom_bar() +
    labs(title = "Target sRNA Chromosome Distribution",
         x = "Chromosome", y = "Count", fill = "sRNA Type") +
    theme_bw(base_size = 13)
ggsave(file.path(opt$outdir, "2_chromosome_distribution.pdf"), p2, width = 9, height = 5)
ggsave(file.path(opt$outdir, "2_chromosome_distribution.png"), p2, width = 9, height = 5, dpi = 150)

# =============================================================================
# 图3: 相对基因位置分布
# =============================================================================
message(">>> [3] 相对基因位置 ...")
srna_anno <- srna[srna$nearest_gene_id != "NA" & !is.na(srna$gene_direction), ]
if (nrow(srna_anno) > 0) {
    p3 <- ggplot(srna_anno, aes(x = gene_direction, fill = gene_direction)) +
        geom_bar() +
        geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5, size = 4) +
        scale_fill_manual(values = c(upstream = "#4C72B0", downstream = "#DD8452",
                                     overlap = "#55A868")) +
        labs(title = "sRNA Position Relative to Nearest Gene",
             x = "Position", y = "Count") +
        theme_bw(base_size = 13) +
        theme(legend.position = "none")
    ggsave(file.path(opt$outdir, "3_gene_direction.pdf"), p3, width = 7, height = 5)
    ggsave(file.path(opt$outdir, "3_gene_direction.png"), p3, width = 7, height = 5, dpi = 150)
}

# =============================================================================
# 图4: 到基因距离分布
# =============================================================================
message(">>> [4] 距离分布 ...")
srna_dist <- srna_anno[srna_anno$gene_distance != "NA", ]
srna_dist$gene_distance <- as.numeric(as.character(srna_dist$gene_distance))
srna_dist <- srna_dist[!is.na(srna_dist$gene_distance), ]
if (nrow(srna_dist) > 0) {
    p4 <- ggplot(srna_dist, aes(x = gene_distance)) +
        geom_histogram(bins = 30, fill = "#4C72B0", color = "white") +
        labs(title = "Distance from sRNA to Nearest Gene",
             x = "Distance (bp)", y = "Count") +
        theme_bw(base_size = 13)
    ggsave(file.path(opt$outdir, "4_distance_distribution.pdf"), p4, width = 8, height = 5)
    ggsave(file.path(opt$outdir, "4_distance_distribution.png"), p4, width = 8, height = 5, dpi = 150)
}

# =============================================================================
# 图5: 目标 sRNA 表达热图 (RPM z-score)
# =============================================================================
message(">>> [5] 目标 sRNA 表达热图 ...")
if (file.exists(opt$rpm) && file.exists(opt$meta)) {
    rpm  <- read.csv(opt$rpm, row.names = 1, check.names = FALSE)
    meta <- read.csv(opt$meta, row.names = 1)

    # 目标 sRNA 的 ID: chrom_start_end
    srna$srna_id <- paste(srna$chrom, srna$start, srna$end, sep = "_")
    target_ids <- intersect(srna$srna_id, rownames(rpm))
    message("  匹配到 RPM 的目标 sRNA: ", length(target_ids))

    if (length(target_ids) >= 2) {
        mat <- as.matrix(rpm[target_ids, , drop = FALSE])
        # log2(RPM+1) 后 z-score
        mat_log <- log2(mat + 1)
        mat_z <- t(scale(t(mat_log)))
        mat_z[is.na(mat_z)] <- 0

        # 列注释 (分组)
        anno_col <- data.frame(Group = meta[colnames(mat), "group"],
                               row.names = colnames(mat))

        pdf(file.path(opt$outdir, "5_target_sRNA_heatmap.pdf"),
            width = 8, height = min(20, 4 + length(target_ids) * 0.12))
        pheatmap(mat_z,
                 annotation_col = anno_col,
                 show_rownames = (length(target_ids) <= 60),
                 clustering_method = "ward.D2",
                 color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
                 main = paste0("Target sRNA Expression (z-score, n=", length(target_ids), ")"),
                 fontsize_row = 6)
        dev.off()
    }
}

# =============================================================================
# 图6: Venn 图 (sRNA 基因 vs DMR 趋势一致基因)
# =============================================================================
message(">>> [6] Venn 图 ...")
if (file.exists(opt$intersect)) {
    inter <- read.csv(opt$intersect, stringsAsFactors = FALSE)

    # sRNA 基因
    srna_genes <- unique(inter$nearest_gene_id[inter$nearest_gene_id != "NA"])

    # 各修饰趋势一致基因
    mod_gene_list <- list()
    mod_gene_list[["sRNA"]] <- srna_genes
    for (mod in MODS) {
        col <- paste0(mod, "_consistent")
        if (col %in% colnames(inter)) {
            g <- unique(inter$nearest_gene_id[inter[[col]] == 1 & inter$nearest_gene_id != "NA"])
            if (length(g) > 0) mod_gene_list[[toupper(mod)]] <- g
        }
    }

    venn_ok <- requireNamespace("VennDiagram", quietly = TRUE)
    if (venn_ok && length(mod_gene_list) >= 2 && length(mod_gene_list) <= 5) {
        suppressMessages(library(VennDiagram))
        try(futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger"), silent = TRUE)
        venn_cols <- brewer.pal(max(3, length(mod_gene_list)), "Set2")[1:length(mod_gene_list)]
        venn.diagram(
            x = mod_gene_list,
            filename = file.path(opt$outdir, "6_venn_sRNA_DMR.png"),
            imagetype = "png",
            fill = venn_cols,
            alpha = 0.5,
            cex = 1.2, cat.cex = 1.2,
            main = "sRNA genes vs DMR-consistent genes"
        )
    } else {
        # 后备: 各修饰一致基因数柱状图
        df_bar <- data.frame(
            Set = names(mod_gene_list),
            N = sapply(mod_gene_list, length)
        )
        p6 <- ggplot(df_bar, aes(x = reorder(Set, -N), y = N, fill = Set)) +
            geom_col() +
            geom_text(aes(label = N), vjust = -0.5) +
            labs(title = "Gene counts: sRNA vs DMR-consistent",
                 x = "", y = "Gene count") +
            theme_bw(base_size = 13) + theme(legend.position = "none")
        ggsave(file.path(opt$outdir, "6_gene_counts_bar.pdf"), p6, width = 7, height = 5)
        ggsave(file.path(opt$outdir, "6_gene_counts_bar.png"), p6, width = 7, height = 5, dpi = 150)
        message("  [INFO] 未安装 VennDiagram, 改用柱状图。安装: install.packages('VennDiagram')")
    }
}

# =============================================================================
# 图7: 候选基因趋势一致性热图
# =============================================================================
message(">>> [7] 候选基因趋势热图 ...")
if (file.exists(opt$consistent)) {
    cons <- read.csv(opt$consistent, stringsAsFactors = FALSE)
    if (nrow(cons) > 0) {
        # sRNA state 列
        srna_state_cols <- c("OE1_state", "OE2_state", "KO6bp_state", "KO8bp_state")
        # DMR state 列 (例如 5mc_*_state)
        dmr_state_cols <- grep("_state$", colnames(cons), value = TRUE)
        dmr_state_cols <- setdiff(dmr_state_cols, srna_state_cols)

        all_state_cols <- c(srna_state_cols, dmr_state_cols)
        all_state_cols <- all_state_cols[all_state_cols %in% colnames(cons)]

        state_mat <- as.matrix(cons[, all_state_cols, drop = FALSE])
        rownames(state_mat) <- paste0(cons$nearest_gene_id, " (", cons$sRNA_type, ")")

        # 重命名列
        colnames(state_mat) <- gsub("_state", "", colnames(state_mat))
        is_srna <- colnames(state_mat) %in% c("OE1", "OE2", "KO6bp", "KO8bp")
        colnames(state_mat)[is_srna] <- paste0("sRNA_", colnames(state_mat)[is_srna])

        state_mat[is.na(state_mat)] <- 3

        pdf(file.path(opt$outdir, "7_candidate_trend_heatmap.pdf"),
            width = 10, height = max(3, nrow(state_mat) * 0.6 + 1.5))
        pheatmap(state_mat,
                 cluster_rows = FALSE, cluster_cols = FALSE,
                 display_numbers = TRUE, number_format = "%.0f",
                 color = c("#D62728", "#1F77B4", "#CCCCCC"),
                 breaks = c(-0.5, 0.5, 1.5, 3.5),
                 main = "Candidate Genes: sRNA & DMR State\n(0=up/hyper red, 1=down/hypo blue, 3=NS grey)",
                 fontsize = 10, fontsize_row = 9)
        dev.off()
        message("  候选基因行数: ", nrow(state_mat))
    }
}

# =============================================================================
# 图8: 候选基因 sRNA log2FC 趋势图
# =============================================================================
message(">>> [8] 候选基因 log2FC 趋势 ...")
if (file.exists(opt$consistent)) {
    cons <- read.csv(opt$consistent, stringsAsFactors = FALSE)
    if (nrow(cons) > 0) {
        plot_df <- data.frame()
        groups <- c("OE1", "OE2", "KO6bp", "KO8bp")
        for (g in groups) {
            lfc_col <- paste0(g, "_log2FC")
            if (lfc_col %in% colnames(cons)) {
                tmp <- data.frame(
                    gene = paste0(cons$nearest_gene_id, "_", cons$start),
                    group = g,
                    log2FC = cons[[lfc_col]]
                )
                plot_df <- rbind(plot_df, tmp)
            }
        }
        plot_df$group <- factor(plot_df$group, levels = groups)

        p8 <- ggplot(plot_df, aes(x = group, y = log2FC, group = gene, color = gene)) +
            geom_line(alpha = 0.6) +
            geom_point(size = 2.5) +
            geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
            labs(title = "Candidate sRNA log2FC across groups",
                 x = "Group (vs WT)", y = "sRNA log2FoldChange", color = "sRNA (gene_start)") +
            theme_bw(base_size = 13)
        ggsave(file.path(opt$outdir, "8_candidate_log2FC_trend.pdf"), p8, width = 9, height = 6)
        ggsave(file.path(opt$outdir, "8_candidate_log2FC_trend.png"), p8, width = 9, height = 6, dpi = 150)
    }
}

message("\n=== 绘图完成 ===")
message("输出目录: ", opt$outdir)
message("  1_sRNA_type_distribution    - sRNA 类型分布")
message("  2_chromosome_distribution   - 染色体分布")
message("  3_gene_direction            - 相对基因位置")
message("  4_distance_distribution     - 到基因距离")
message("  5_target_sRNA_heatmap       - 目标 sRNA 表达热图")
message("  6_venn_sRNA_DMR (或柱状图)  - sRNA vs DMR 基因交集")
message("  7_candidate_trend_heatmap   - 候选基因趋势热图")
message("  8_candidate_log2FC_trend    - 候选基因 log2FC 趋势")
