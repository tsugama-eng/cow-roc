# 必要なパッケージ
library(ggplot2)
library(ggnewscale)
library(ggrepel)
library(colorspace)
library(dplyr)
library(ggbeeswarm)
library(uwot)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)
library(uwot)

# アミノ酸略号とコドンの対応をリストで定義
codon_table_sorted <- list(
  F = c("TTC", "TTT"), 
  L = c("CTC", "CTG", "CTA", "CTT", "TTG", "TTA"),
  I = c("ATC", "ATA", "ATT"),
  M = c("ATG"),
  V = c("GTC", "GTG", "GTA", "GTT"),
  S = c("AGC", "AGT", "TCC", "TCG", "TCA", "TCT"),
  P = c("CCC", "CCG", "CCA", "CCT"),
  T = c("ACC", "ACG", "ACA", "ACT"),
  A = c("GCC", "GCG", "GCA", "GCT"),
  Y = c("TAC", "TAT"),
  H = c("CAC", "CAT"),
  Q = c("CAG", "CAA"),
  N = c("AAC", "AAT"),
  K = c("AAG", "AAA"),
  D = c("GAC", "GAT"),
  E = c("GAG", "GAA"),
  C = c("TGC", "TGT"),
  W = c("TGG"),
  R = c("AGG", "AGA", "CGC", "CGG", "CGA", "CGT"),
  G = c("GGC", "GGG", "GGA", "GGT"),
  Stop = c("TAG", "TAA", "TGA") # 並べ替えルールに従って C→G→A→T
)

codon_order <- unlist(codon_table_sorted)

codon_to_aa <- stack(codon_table_sorted)
colnames(codon_to_aa) <- c("codon", "amino_acid")
rownames(codon_to_aa) <- codon_to_aa$codon

aa_addition0 <- function (vdata){
    for (i in codon_to_aa$codon) {
      aa_to_use <- codon_to_aa[codon_to_aa$codon == i, "amino_acid"]
      # 条件に合う要素を抽出
      idx <- grepl(i, vdata) & !grepl("\\(", vdata)
      # 括弧付きアミノ酸を追加
      vdata[idx] <- paste0(vdata[idx], " (", aa_to_use, ")")
    }
    return(vdata)
}

#####################################
#Heatmap with the coeffs with thr0

species_labels = c(
  
  # 植物（果物・その他単子葉/双子葉）
  "Arabidopsis_thaliana"      = "Arabidopsis\nthaliana",
  "Brassica_rapa"             = "Brassica\nrapa",
  "Vitis_vinifera"            = "Vitis\nvinifera",
  "Malus_domestica"           = "Malus\ndomestica",
  "Solanum_lycopersicum"      = "Solanum\nlycopersicum",
  "Glycine_max"               = "Glycine\nmax",
  "Populus_trichocarpa"       = "Populus\ntrichocarpa",
  "Citrus_clementina"         = "Citrus\nclementina",
  "Chenopodium_quinoa"        = "Chenopodium\nquinoa",
  "Allium_cepa"               = "Allium\ncepa",
  "Asparagus_officinalis"     = "Asparagus\nofficinalis",
  "Musa_acuminata"            = "Musa\nacuminata",
  "Ananas_comosus"            = "Ananas\ncomosus",

  # 植物（主要・穀物）
  "Oryza_sativa"              = "Oryza\nsativa",
  "Triticum_aestivum"         = "Triticum\naestivum",
  "Hordeum_vulgare"           = "Hordeum\nvulgare",
  "Zea_mays"                  = "Zea\nmays",
  "Cenchrus_americanus"       = "Cenchrus\namericanus",
  "Sorghum_bicolor"           = "Sorghum\nbicolor",
  
  # 植物（苔・藻類）
  "Physcomitrium_patens"      = "Physcomitrium\npatens",
  "Chlamydomonas_reinhardtii" = "Chlamydomonas\nreinhardtii",
  
  # 真菌・動物（モデル生物）
  "Homo_sapiens"              = "Homo\nsapiens",
  "Mus_musculus"              = "Mus\nmusculus",
  "Drosophila_melanogaster"   = "Drosophila\nmelanogaster",
  "Saccharomyces_cerevisiae"  = "Saccharomyces\ncerevisiae",

  # 細菌
  "Escherichia_coli"          = "Escherichia\ncoli",
  "Salmonella_enterica"       = "Salmonella\nenterica",
  "Bacillus_subtilis"         = "Bacillus\nsubtilis",
  "Mycobacterium_tuberculosis"= "Mycobacterium\ntuberculosis"
)

df <- read.table("all_species_sum.coeff.linear.txt", header=T, sep="\t")

# 1. thr0 の行だけを抽出し、生物種名を整理
plot_df <- df %>%
  dplyr::filter(grepl("_thr0$", species_thr)) %>%
  dplyr::mutate(species = gsub("_thr0$", "", species_thr)) %>%
  dplyr::select(species, everything(), -species_thr)

# 2. 行列形式に変換
mat <- as.matrix(plot_df[,-1])
rownames(mat) <- plot_df$species

# 3. 転置（縦：コドン、横：生物種）
mat_t <- t(mat)

# 1. ラベルリストの作成（改行を空白に置換）
# species_labels の名前（Oryza_sativa 等）を順序制御に使用します
target_order <- names(species_labels)
clean_labels <- gsub("\n", " ", species_labels)

# 2. データの整形
# mat_t が「縦：コドン、横：生物種」の行列であると仮定します
# 行列の列名（生物種）が target_order に含まれるものだけに絞り、順序を合わせる
common_species <- intersect(target_order, colnames(mat_t))
mat_t_ordered <- mat_t[, common_species]
rownames(mat_t_ordered) <- aa_addition0(rownames(mat_t_ordered))

# 表示用ラベルも順序を合わせる
display_labels <- clean_labels[common_species]

# 3. カラーレンジの設定
contrast_factor <- 0.6  # この値を小さくするほど色が濃くなります
limit_val <- max_val * contrast_factor

col_fun = colorRamp2(
  c(-limit_val, 0, limit_val), 
  c("royalblue", "white", "darkorange")
)

# 1. 保存用の設定
output_file <- "Codon_Weight_Heatmap_Final_v2.png"

# 2. ヒートマップの定義
ht <- Heatmap(mat_t_ordered, 
              name = "Codon Weight", 
              col = col_fun,
              cluster_columns = FALSE, 
              column_order = 1:ncol(mat_t_ordered),
              cluster_rows = TRUE, 
              column_labels = display_labels,
              show_row_names = TRUE,
              
              # 横軸タイトル（全体タイトル）を非表示
              column_title = NULL, 
              row_title = "Codons",
              
              row_names_gp = gpar(fontsize = 8, fontface = "plain"),
              column_names_gp = gpar(fontsize = 9, fontface = "italic"),
              
              # 凡例（スケールバー）の設定
              heatmap_legend_param = list(
                title = "Codon weight",
                title_position = "topcenter",
                direction = "horizontal",
                legend_width = unit(3, "cm"),
                # タイトルのフォントスタイルを plain に指定
                title_gp = gpar(fontsize = 10, fontface = "plain"),
                # ラベル（数字）のフォントスタイルも plain に指定
                labels_gp = gpar(fontsize = 9, fontface = "plain"),
                at = c(-limit_val, 0, limit_val),
                labels = c(paste0("< ", round(-limit_val, 2)), "0", paste0("> ", round(limit_val, 2)))
              ))

# 3. 描画と保存
png(output_file, width = 16, height = 27, units = "cm", res = 600)

# スケールバーを上部に配置
draw(ht, 
     heatmap_legend_side = "top", 
     padding = unit(c(10, 10, 35, 10), "mm")) # 下, 左, 上, 右

dev.off()

message(paste("Saved to:", output_file))

#############################################
#UMAP
library(uwot)
library(ggplot2)
library(ggrepel)
library(dplyr)

# 1. データの準備 (mat は行：生物種、列：コドンの行列)
# すでに mat が作成されている前提です
data <- mat
group_s <- rownames(data) # 生物種名（アンダーバー形式）
sel_species <- names(species_labels) # 系統順のリスト

# 2. 系統グループの定義（29種を系統ごとに色分け・形状分けするため）
# case_when を使って生物種をカテゴリーに分類
species_meta <- data.frame(Group_s = group_s) %>%
  mutate(Category = case_when(
    Group_s %in% c("Oryza_sativa", "Zea_mays", "Sorghum_bicolor", "Cenchrus_americanus", "Hordeum_vulgare", "Triticum_aestivum") ~ "Poaceae",
    Group_s %in% c("Ananas_comosus", "Musa_acuminata", "Asparagus_officinalis", "Allium_cepa") ~ "Other monocots",
    Group_s %in% c("Populus_trichocarpa", "Malus_domestica", "Vitis_vinifera", "Citrus_clementina", "Arabidopsis_thaliana", "Brassica_rapa", "Solanum_lycopersicum", "Glycine_max", "Chenopodium_quinoa") ~ "Dicots",
    Group_s %in% c("Physcomitrium_patens", "Chlamydomonas_reinhardtii") ~ "Non-vascular/Algae",
    Group_s %in% c("Saccharomyces_cerevisiae", "Homo_sapiens", "Mus_musculus", "Drosophila_melanogaster") ~ "Animals/Fungi",
    Group_s %in% c("Escherichia_coli", "Bacillus_subtilis", "Mycobacterium_tuberculosis", "Salmonella_enterica") ~ "Bacteria"
  ))

# 3. UMAP計算
set.seed(123)
umap_all <- uwot::umap(
    data,              # 行：生物種(29), 列：コドン(64)の行列
    n_neighbors = 15,
    min_dist = 0.15
)

# 2. プロット用データフレームの構築
# species_meta$Category が作成済みであることを前提としています
scores_all <- data.frame(
    UMAP1 = umap_all[, 1],
    UMAP2 = umap_all[, 2],
    Group_s = group_s,
    Category = factor(species_meta$Category, 
                      levels = c("Poaceae", "Other monocots", "Dicots", 
                                 "Non-vascular/Algae", "Animals/Fungi", "Bacteria"))
)

# 表示用の学名を整理（アンダーバーをスペースに）
scores_all$Label <- gsub("_", " ", scores_all$Group_s)

# 3. 可視化のパラメータ設定
strk <- 0.8
ap1  <- 0.8

# 色と形状のカスタム設定
group_colors <- c("Poaceae" = "#E41A1C", "Other monocots" = "#FF7F00", "Dicots" = "#4DAF4A", 
                  "Non-vascular/Algae" = "#984EA3", "Animals/Fungi" = "#377EB8", "Bacteria" = "#A65628")
group_shapes <- c("Poaceae" = 16, "Other monocots" = 17, "Dicots" = 18, 
                  "Non-vascular/Algae" = 15, "Animals/Fungi" = 8, "Bacteria" = 4)

# --- プロット作成 ---
up1 <- ggplot(scores_all,
               aes(x = UMAP1, y = UMAP2,
                   color = Category, shape = Category)) +

    # 種名ラベルの追加
    # color = Category を aes 内に入れることで点の色と同期
    # show.legend = FALSE で凡例の「a」を消去
    geom_text_repel(aes(label = Label, color = Category), 
                    alpha = 0.65,
                    size = 3.5, 
                    fontface = "italic",
                    max.overlaps = Inf,
                    show.legend = FALSE) + 
                    

    # プロット本体（点）
    geom_point(size = 4, alpha = ap1, stroke = strk) +
    
    # 色と形状の適用
    scale_shape_manual(values = group_shapes, name = "Group") +
    scale_color_manual(values = group_colors, name = "Group") +
    

    # デザイン調整
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      aspect.ratio = 1, # UMAPの解釈を正しくするためにアスペクト比を1に固定
      # --- 凡例をプロット内部に寄せる設定 ---
      legend.position = c(0.02, 0.02), # 左下(0,0)から右上(1,1)の座標指定。適宜調整。
      legend.justification = c("left", "bottom"), # 配置の基準点
      legend.background = element_rect(fill = alpha("white", 0), color = NA), # 半透明の背景
      legend.key.size = unit(4, "mm"), # 凡例アイコンを小さくして節約
      legend.spacing.y = unit(2.5, "mm"),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10, face = "plain")
    ) +
    guides(
      color = guide_legend(byrow = TRUE),
      shape = guide_legend(byrow = TRUE)
    ) +
    #labs(title = "Codon weight UMAP plot (2026 Analysis)", 
    #     x = "UMAP1", y = "UMAP2")
    labs(x = "UMAP1", y = "UMAP2")


# 表示と保存
print(up1)
ggsave("codon_weight_UMAP_w_species_names_2026.png", plot = up1, width = 14, height = 14, dpi = 450, units="cm", bg = "white")

# --- プロット2: 生物種名なし ---
up2 <- ggplot(scores_all,
               aes(x = UMAP1, y = UMAP2,
                   color = Category, shape = Category)) +
    geom_point(size = 4, alpha = ap1, stroke = strk) +
    scale_shape_manual(values = group_shapes, name = "Group") +
    scale_color_manual(values = group_colors, name = "Group") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      legend.text = element_text(size = 10),
      aspect.ratio = 1
    ) +
    labs(title = "Codon weight UMAP plot (without species names)", 
         x = "UMAP1", y = "UMAP2")

print(up2)
ggsave("codon_weight_UMAP_wo_species_names_2026.png", plot = up2, width = 9, height = 7, dpi = 450, bg = "white")



