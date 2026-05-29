#!/usr/bin/env Rscript
# =============================================================================
# plot_sRNA_analysis.R - sRNA × DMR 联合分析可视化 (出版级)
# =============================================================================
# 风格: Nature (NPG) 配色 + Arial 字体 + 25pt 加粗 + 清晰边框
#
# 字体说明:
#   全部图表统一使用 Arial。PDF 用 cairo_pdf 设备 (标准 pdf 设备不支持
#   Arial)。若系统无 Arial 字体, cairo 会回退到默认无衬线字体, 或:
#     conda install -c conda-forge font-ttf-liberation   # Arial 开源替代
#   也可将下方 FONT 改为 "sans"。
#
# 生成图表:
#   1. sRNA 类型分布柱状图        2. sRNA 染色体分布图
#   3. sRNA 相对基因位置分布      4. sRNA 到基因距离分布直方图
#   5. 目标 sRNA 表达热图         6. Venn 图 (sRNA vs DMR 一致基因)
#   7. 候选基因趋势一致性热图     8. 候选基因 sRNA log2FC 趋势图
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
    library(grid)
})

# ===== 参数 =====
option_list <- list(
    make_option("--srna", type = "character",
                default = "/public/home/h14166/fang/Heyufei/smrna_seq/sRNA_diff_analysis/results/target_sRNA_annotated.csv"),
    make_option("--rpm", type = "character",
                default = "/public/home/h14166/fang/Heyufei/smrna_seq/sRNA_diff_analysis/counts/sRNA_rpm.csv"),
    make_option("--meta", type = "character",
                default = "/public/home/h14166/fang/Heyufei/smrna_seq/sRNA_diff_analysis/counts/sample_metadata.csv"),
    make_option("--intersect", type = "character",
                default = "/public/home/h14166/fang/Heyufei/Unite/sRNA_DMR_intersect.csv"),
    make_option("--consistent", type = "character",
                default = "/public/home/h14166/fang/Heyufei/Unite/sRNA_DMR_intersect_consistent.csv"),
    make_option("--outdir", type = "character", default = "plots")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 出版级配置: Arial 字体 + Nature (NPG) 配色 + 25pt 加粗主题
# =============================================================================

FONT <- "Arial"   # 统一字体

# NPG (Nature Publishing Group) 经典配色
NPG <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
         "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85")

# 出版级 ggplot 主题: Arial, 25pt, 加粗, 清晰边框
BASE_SIZE <- 25
theme_pub <- function(base_size = BASE_SIZE) {
    theme_bw(base_size = base_size, base_family = FONT) +
    theme(
        text             = element_text(family = FONT, face = "bold", colour = "black"),
        plot.title       = element_text(family = FONT, face = "bold", size = base_size,
                                        hjust = 0.5, margin = margin(b = 12)),
        plot.subtitle    = element_text(family = FONT, face = "bold", size = base_size - 6,
                                        hjust = 0.5),
        axis.title       = element_text(family = FONT, face = "bold", size = base_size),
        axis.text        = element_text(family = FONT, face = "bold", size = base_size - 4,
                                        colour = "black"),
        axis.ticks       = element_line(colour = "black", linewidth = 1),
        axis.ticks.length = unit(0.25, "cm"),
        panel.border     = element_rect(colour = "black", fill = NA, linewidth = 1.5),
        panel.grid       = element_blank(),
        legend.title     = element_text(family = FONT, face = "bold", size = base_size - 4),
        legend.text      = element_text(family = FONT, face = "bold", size = base_size - 6),
        legend.key.size  = unit(0.9, "cm"),
        plot.margin      = margin(15, 15, 15, 15)
    )
}
LABEL_SIZE <- 8   # geom_text 标签 (≈ 23pt 加粗)

# 统一保存: PDF 用 cairo_pdf (支持 Arial), PNG 300 dpi
save_plot <- function(base, plot, width, height) {
    ggsave(file.path(opt$outdir, paste0(base, ".pdf")), plot,
           width = width, height = height, device = cairo_pdf, limitsize = FALSE)
    ggsave(file.path(opt$outdir, paste0(base, ".png")), plot,
           width = width, height = height, dpi = 300, limitsize = FALSE)
}

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
    geom_bar(width = 0.7, colour = "black", linewidth = 0.8) +
    geom_text(stat = "count", aes(label = after_stat(count)),
              hjust = -0.25, size = LABEL_SIZE, fontface = "bold", family = FONT) +
    coord_flip(clip = "off") +
    scale_fill_manual(values = NPG) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(title = "Target sRNA Type Distribution", x = NULL, y = "Count") +
    theme_pub() +
    theme(legend.position = "none")
save_plot("1_sRNA_type_distribution", p1, 11, 7)

# =============================================================================
# 图2: 染色体分布
# =============================================================================
message(">>> [2] 染色体分布 ...")
srna$chrom <- factor(srna$chrom, levels = sort(unique(as.character(srna$chrom))))
p2 <- ggplot(srna, aes(x = chrom, fill = sRNA_type)) +
    geom_bar(width = 0.75, colour = "black", linewidth = 0.6) +
    scale_fill_manual(values = NPG) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(title = "Target sRNA Chromosome Distribution",
         x = "Chromosome", y = "Count", fill = "sRNA Type") +
    theme_pub()
save_plot("2_chromosome_distribution", p2, 13, 7)

# =============================================================================
# 图3: 相对基因位置分布
# =============================================================================
message(">>> [3] 相对基因位置 ...")
srna_anno <- srna[srna$nearest_gene_id != "NA" & !is.na(srna$gene_direction), ]
if (nrow(srna_anno) > 0) {
    p3 <- ggplot(srna_anno, aes(x = gene_direction, fill = gene_direction)) +
        geom_bar(width = 0.7, colour = "black", linewidth = 0.8) +
        geom_text(stat = "count", aes(label = after_stat(count)),
                  vjust = -0.4, size = LABEL_SIZE, fontface = "bold", family = FONT) +
        scale_fill_manual(values = c(upstream = "#3C5488", downstream = "#E64B35",
                                     overlap = "#00A087")) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
        labs(title = "sRNA Position Relative to Nearest Gene",
             x = NULL, y = "Count") +
        theme_pub() +
        theme(legend.position = "none")
    save_plot("3_gene_direction", p3, 9, 7)
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
        geom_histogram(bins = 30, fill = "#4DBBD5", colour = "black", linewidth = 0.6) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
        labs(title = "Distance from sRNA to Nearest Gene",
             x = "Distance (bp)", y = "Count") +
        theme_pub()
    save_plot("4_distance_distribution", p4, 11, 7)
}

# =============================================================================
# 图5: 目标 sRNA 表达热图 (RPM z-score)
# =============================================================================
message(">>> [5] 目标 sRNA 表达热图 ...")
if (file.exists(opt$rpm) && file.exists(opt$meta)) {
    rpm  <- read.csv(opt$rpm, row.names = 1, check.names = FALSE)
    meta <- read.csv(opt$meta, row.names = 1)

    srna$srna_id <- paste(srna$chrom, srna$start, srna$end, sep = "_")
    target_ids <- intersect(srna$srna_id, rownames(rpm))
    message("  匹配到 RPM 的目标 sRNA: ", length(target_ids))

    if (length(target_ids) >= 2) {
        mat <- as.matrix(rpm[target_ids, , drop = FALSE])
        mat_z <- t(scale(t(log2(mat + 1))))
        mat_z[is.na(mat_z)] <- 0

        anno_col <- data.frame(Group = meta[colnames(mat), "group"],
                               row.names = colnames(mat))
        grp_levels <- unique(anno_col$Group)
        anno_colors <- list(Group = setNames(NPG[seq_along(grp_levels)], grp_levels))

        hm_color <- colorRampPalette(c("#3C5488", "#4DBBD5", "white", "#F39B7F", "#E64B35"))(100)

        # cairo_pdf(family=Arial): pheatmap 文字继承该字体
        cairo_pdf(file.path(opt$outdir, "5_target_sRNA_heatmap.pdf"),
                  width = 11, height = min(24, 5 + length(target_ids) * 0.14),
                  family = FONT)
        pheatmap(mat_z,
                 annotation_col = anno_col,
                 annotation_colors = anno_colors,
                 show_rownames = (length(target_ids) <= 60),
                 clustering_method = "ward.D2",
                 color = hm_color,
                 border_color = NA,
                 main = paste0("Target sRNA Expression (z-score, n=", length(target_ids), ")"),
                 fontsize = 18, fontsize_row = 10, fontsize_col = 16,
                 treeheight_row = 40, treeheight_col = 30,
                 angle_col = 45)
        dev.off()
    }
}

# =============================================================================
# 图6: Venn 图 (sRNA 基因 vs DMR 一致基因)
# =============================================================================
message(">>> [6] Venn 图 ...")
if (file.exists(opt$intersect)) {
    inter <- read.csv(opt$intersect, stringsAsFactors = FALSE)
    srna_genes <- unique(inter$nearest_gene_id[inter$nearest_gene_id != "NA"])

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
        venn_cols <- NPG[seq_along(mod_gene_list)]
        venn.diagram(
            x = mod_gene_list,
            filename = file.path(opt$outdir, "6_venn_sRNA_DMR.png"),
            imagetype = "png", height = 2400, width = 2400, resolution = 300,
            fill = venn_cols, alpha = 0.55,
            col = "black", lwd = 3,
            cex = 2.2, fontface = "bold", fontfamily = FONT,
            cat.cex = 2.2, cat.fontface = "bold", cat.fontfamily = FONT,
            main = "sRNA genes vs DMR-consistent genes",
            main.cex = 2.2, main.fontface = "bold", main.fontfamily = FONT
        )
    } else {
        df_bar <- data.frame(Set = names(mod_gene_list), N = sapply(mod_gene_list, length))
        p6 <- ggplot(df_bar, aes(x = reorder(Set, -N), y = N, fill = Set)) +
            geom_col(width = 0.7, colour = "black", linewidth = 0.8) +
            geom_text(aes(label = N), vjust = -0.4, size = LABEL_SIZE, fontface = "bold",
                      family = FONT) +
            scale_fill_manual(values = NPG) +
            scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
            labs(title = "Gene Counts: sRNA vs DMR-consistent", x = NULL, y = "Gene Count") +
            theme_pub() + theme(legend.position = "none")
        save_plot("6_gene_counts_bar", p6, 10, 7)
        message("  [INFO] 未安装 VennDiagram, 改用柱状图。安装: conda install -c conda-forge r-venndiagram")
    }
}

# =============================================================================
# 图7: 候选基因趋势一致性热图 (ggplot geom_tile, 出版级)
# =============================================================================
message(">>> [7] 候选基因趋势热图 ...")
if (file.exists(opt$consistent)) {
    cons <- read.csv(opt$consistent, stringsAsFactors = FALSE)
    if (nrow(cons) > 0) {
        srna_state_cols <- c("OE1_state", "OE2_state", "KO6bp_state", "KO8bp_state")
        dmr_state_cols <- grep("_state$", colnames(cons), value = TRUE)
        dmr_state_cols <- setdiff(dmr_state_cols, srna_state_cols)

        all_state_cols <- c(srna_state_cols, dmr_state_cols)
        all_state_cols <- all_state_cols[all_state_cols %in% colnames(cons)]

        col_labels <- gsub("_state", "", all_state_cols)
        is_srna <- col_labels %in% c("OE1", "OE2", "KO6bp", "KO8bp")
        col_labels[is_srna] <- paste0("sRNA_", col_labels[is_srna])

        row_labels <- paste0(cons$nearest_gene_id, "\n(", cons$sRNA_type, ")")

        hm_df <- data.frame()
        for (j in seq_along(all_state_cols)) {
            vals <- cons[[all_state_cols[j]]]
            vals[is.na(vals)] <- 3
            hm_df <- rbind(hm_df, data.frame(
                gene  = row_labels,
                group = col_labels[j],
                state = as.integer(vals)
            ))
        }
        hm_df$state_label <- factor(hm_df$state, levels = c(0, 1, 3),
                                    labels = c("Up/Hyper", "Down/Hypo", "NS"))
        hm_df$group <- factor(hm_df$group, levels = col_labels)
        hm_df$gene  <- factor(hm_df$gene, levels = rev(unique(row_labels)))
        hm_df$num   <- hm_df$state

        state_cols <- c("Up/Hyper" = "#E64B35", "Down/Hypo" = "#3C5488", "NS" = "#B0B0B0")

        p7 <- ggplot(hm_df, aes(x = group, y = gene, fill = state_label)) +
            geom_tile(colour = "white", linewidth = 1.5) +
            geom_text(aes(label = num), size = LABEL_SIZE - 1, fontface = "bold",
                      family = FONT, colour = "white") +
            scale_fill_manual(values = state_cols, name = "State") +
            scale_x_discrete(position = "top", expand = c(0, 0)) +
            scale_y_discrete(expand = c(0, 0)) +
            labs(title = "Candidate Genes: sRNA & DMR State",
                 subtitle = "0 = Up/Hyper, 1 = Down/Hypo, 3 = NS",
                 x = NULL, y = NULL) +
            coord_equal() +
            theme_pub() +
            theme(
                axis.text.x  = element_text(angle = 45, hjust = 0, face = "bold", family = FONT),
                panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
                axis.ticks   = element_blank()
            )
        n_row <- length(unique(row_labels))
        n_col <- length(all_state_cols)
        save_plot("7_candidate_trend_heatmap", p7, 3 + n_col * 1.6, 3 + n_row * 1.3)
        message("  候选基因行数: ", nrow(cons))
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
                plot_df <- rbind(plot_df, data.frame(
                    gene = paste0(cons$nearest_gene_id, "_", cons$start),
                    group = g, log2FC = cons[[lfc_col]]))
            }
        }
        plot_df$group <- factor(plot_df$group, levels = groups)

        p8 <- ggplot(plot_df, aes(x = group, y = log2FC, group = gene, colour = gene)) +
            geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 1) +
            geom_line(linewidth = 1.5, alpha = 0.85) +
            geom_point(size = 5) +
            scale_colour_manual(values = NPG) +
            labs(title = "Candidate sRNA log2FC across Groups",
                 x = "Group (vs WT)", y = "sRNA log2FoldChange", colour = "sRNA") +
            theme_pub() +
            theme(legend.text = element_text(size = BASE_SIZE - 9, family = FONT))
        save_plot("8_candidate_log2FC_trend", p8, 12, 7.5)
    }
}

message("\n=== 绘图完成 (出版级: Arial + NPG配色 + 25pt加粗) ===")
message("输出目录: ", opt$outdir)
message("  1_sRNA_type_distribution    2_chromosome_distribution")
message("  3_gene_direction            4_distance_distribution")
message("  5_target_sRNA_heatmap       6_venn_sRNA_DMR (或柱状图)")
message("  7_candidate_trend_heatmap   8_candidate_log2FC_trend")
