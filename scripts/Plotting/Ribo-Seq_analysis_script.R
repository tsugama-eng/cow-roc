library(ggplot2)
library(ggpubr)
library(scales)
library(MASS)       # rlm
library(quantreg)   # rq
library(fields)
library(MethComp)
library(ggpointdensity)
library(pheatmap)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ComplexHeatmap)
library(circlize) # 色の範囲指定用
library(grid)
library(ggrepel) 

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
codon_order_wo_stop <- setdiff(codon_order,c("TAG","TGA","TAA"))

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



#Testing

coeffs<-read.table("Arabidopsis_thaliana_all_nuc_gene_tpm.mean.linear.coeffs.tsv", row.names=1, header=T, sep="\t")
coeffs<-coeffs[codon_order_wo_stop,,drop=F]

counts<-read.table("testpileup_nuc.txt", row.names=1, header=T, sep="\t")
counts001<-read.table("testpileup_nuc_sub001.txt", row.names=1, header=T, sep="\t")

counts<-counts[codon_order_wo_stop,,drop=F]
counts001<-counts001[codon_order_wo_stop,,drop=F]

for (i in 1:ncol(counts)){
  res<-cor.test(as.numeric(counts[,i]),as.numeric(unlist(coeffs)))
  print(res)
}
plot(unlist(coeffs),counts[,101])

for (i in 1:ncol(counts001)){
  res<-cor.test(as.numeric(counts001[,i]),as.numeric(unlist(coeffs)))
  print(res)
}
plot(unlist(coeffs),counts001[,101])

#For analyses with all
species<-c("Arabidopsis thaliana","Arabidopsis thaliana pt","Oryza sativa","Homo sapiens","Mus musculus","Drosophila melanogaster","Saccharomyces cerevisiae","Escherichia coli")
species_key<-gsub(" ","_",species)

counts<-read.table("all_ribo_pileup.txt", row.names=NULL, header=T, sep="\t")
coeffs_all<-read.table("all_species_sum.coeff.linear.txt", row.names=1, header=T, sep="\t")

spkey<-"Arabidopsis_thaliana"

# ComplexHeatmap のメッセージ（raster等）を非表示
ht_opt$message = FALSE

for (spkey in species_key) {
  message(paste0(">>> Processing: ", spkey))
  
  # --- 1. データの抽出 ---
  counts_sub <- counts[grepl(spkey, counts[,ncol(counts)]), , drop=F]
  if (spkey == "Arabidopsis_thaliana") {
      counts_sub <- counts_sub[!grepl("_pt", counts_sub[,ncol(counts_sub)]), , drop=F]
  } else if (spkey == "Arabidopsis_thaliana_pt") {
      counts_sub <- counts_sub[grepl("_pt", counts_sub[,ncol(counts_sub)]), , drop=F]
  }
  if (nrow(counts_sub) == 0) next
  
  # --- 2. データ整形と平均計算 ---
  df <- counts_sub
  n_col <- ncol(df)
  # 解析用に純粋なコドン名を保持しておく
  raw_codons <- df[, 1] 
  
  rownames(df) <- make.unique(paste0("row_", 1:nrow(df), "_", df[, n_col]))
  colnames(df) <- c("codon", 1:(n_col-1))
  
  values_only <- df[, 2:(n_col-1)]
  values_only[] <- lapply(values_only, function(x) as.numeric(as.character(x)))
  
  # 平均計算（aa_addition前の純粋なコドン名で集計）
  codon_means_all <- aggregate(values_only, list(Codon = raw_codons), mean, na.rm=TRUE)
  rownames(codon_means_all) <- codon_means_all$Codon
  plot_data_means <- as.matrix(codon_means_all[,-1])

  # 順序決定 (クラスタリング)
  if (nrow(plot_data_means) >= 2) {
      m_scaled <- scale(plot_data_means)
      m_scaled[is.nan(m_scaled)] <- 0 
      current_codon_order <- rownames(m_scaled)[hclust(dist(m_scaled))$order]
  } else {
      current_codon_order <- rownames(plot_data_means)
  }

  sp_name_clean <- gsub("_", " ", spkey)
  my_color <- colorRampPalette(c("royalblue", "white", "darkorange"))(100)

  # --- 3. pheatmap 出力 (全データ版) ---
  # ソート済みの行列を作成
  sorted_indices <- order(factor(df$codon, levels = current_codon_order))
  df_sorted <- df[sorted_indices, ]
  data_matrix_sorted <- as.matrix(df_sorted[, 2:(n_col-1)]) # ★ここで focused 用の変数を確定

  # 表示用にアミノ酸情報を付加したコドン名を作成（解析には使わない）
  df_display <- df_sorted
  df_display$codon <- aa_addition0(df_display$codon)
  
  # 注釈の作成
  anno_row <- data.frame(Codon = factor(df_display$codon))
  rownames(anno_row) <- rownames(df_display)

  graphics.off()
  pheatmap::pheatmap(data_matrix_sorted,
           cluster_rows = FALSE, cluster_cols = FALSE, scale = "row", show_rownames = FALSE,
           main = paste0(sp_name_clean, ": Gene-normalized counts (All)"),
           annotation_row = anno_row, # アミノ酸付きのラベルを注釈に使用
           color = my_color, filename = paste0(spkey, "_pheatmap_all_scaled.png"), width=15, height=18)

  # --- 4. pheatmap 出力 (平均版) ---
  plot_data_means_sorted <- plot_data_means[current_codon_order, , drop=FALSE]
  # 表示用にアミノ酸付加
  rownames(plot_data_means_sorted) <- aa_addition0(rownames(plot_data_means_sorted))
  
  graphics.off()
  pheatmap::pheatmap(plot_data_means_sorted,
           cluster_rows = FALSE, cluster_cols = FALSE, scale = "row",
           main = paste0(sp_name_clean, ": Mean scaled counts"),
           color = my_color, filename = paste0(spkey, "_pheatmap_mean_scaled.png"), width=15, height=12)

  # --- 5. ComplexHeatmap 出力 (スケーリングなし平均版) ---
  plot_data_means_none <- plot_data_means[codon_order_wo_stop, , drop=FALSE]
  # 表示用にアミノ酸付加
  rownames(plot_data_means_none) <- aa_addition0(rownames(plot_data_means_none))
  
  mat_none <- as.matrix(plot_data_means_none)
  max_val <- max(mat_none, na.rm = TRUE)
  col_fun_none = colorRamp2(c(0, max_val/2, max_val), c("royalblue", "white", "darkorange"))

  ht <- Heatmap(mat_none, name = "Density", column_title = NULL,
                cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = TRUE, 
                row_names_gp = gpar(fontsize = 10, fontface = "plain"),
                col = col_fun_none, use_raster = FALSE,
                heatmap_legend_param = list(title = "Normalized Count", title_position = "topcenter",
                  direction = "horizontal", legend_width = unit(10, "cm"),
                  title_gp = gpar(fontsize = 14, fontface = "plain"),
                  labels_gp = gpar(fontsize = 12, fontface = "plain"), at = c(0, round(max_val, 1))))

  graphics.off()
  png(paste0(spkey, "_ComplexHeatmap_mean_none.png"), width = 7200, height = 6000, res = 450)
  draw(ht, heatmap_legend_side = "top", padding = unit(c(35, 10, 10, 10), "mm"))
  grid.text(paste0(sp_name_clean, ": Mean unscaled counts per sample"), x = 0.5, y = 0.99, gp = gpar(fontsize = 22, fontface = "plain"))
  dev.off()
  
  # --- 6. 特定種の Focused Heatmap ---
  target_cols <- NULL
  if (spkey == "Arabidopsis_thaliana" || spkey == "Homo_sapiens") {
    target_cols <- 37:65 #P site is the position 51
  } else if (spkey == "Escherichia_coli") {
    target_cols <- 42:70 #P site is the position 56
  }

  if (!is.null(target_cols)) {
    message(paste0(">>> Creating focused heatmap for: ", spkey))
    
    mat_focus <- as.matrix(plot_data_means_none[, target_cols, drop = FALSE])
    mat_focus[is.na(mat_focus)] <- 0

    max_val <- max(mat_focus, na.rm = TRUE)
    col_fun_focus = colorRamp2(c(0, max_val/2, max_val), c("royalblue", "white", "darkorange"))

    P_site <- (length(target_cols)+1)/2
    highlight_start_rel <- P_site - 5
    if (spkey == "Arabidopsis_thaliana"){highlight_start_rel <- P_site - 4}
    highlight_end_rel <- P_site + 5

    # 強調色を設定するベクトルの作成
    num_cols_focus <- ncol(mat_focus)
    anno_colors_vec <- rep("transparent", num_cols_focus)
    
    # 相対位置に基づいて色を塗る
    if (highlight_start_rel >= 1 && highlight_end_rel <= num_cols_focus) {
        anno_colors_vec[highlight_start_rel:highlight_end_rel] <- "gold"
    }

    # アノテーションの作成
    anno_gold <- HeatmapAnnotation(
      focus_range = anno_colors_vec, 
      col = list(focus_range = c("gold" = "gold", "transparent" = "transparent")),
      show_legend = FALSE,
      show_annotation_name = FALSE,
      simple_anno_size = unit(2, "mm") # バーの太さ
    )

    ht <- Heatmap(mat_focus, name = "Density", column_title = NULL,
                cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = TRUE, show_column_names = FALSE, 
                row_names_gp = gpar(fontsize = 10, fontface = "plain"),
                col = col_fun_focus, use_raster = FALSE,
                bottom_annotation = anno_gold, # 金色のバーを追加
                heatmap_legend_param = list(title = "Count", title_position = "topcenter",
                  direction = "horizontal", legend_width = unit(3, "cm"),
                  title_gp = gpar(fontsize = 14, fontface = "plain"),
                  labels_gp = gpar(fontsize = 12, fontface = "plain"), at = c(0, round(max_val, 1))))

    graphics.off()
    png(paste0(spkey, "_ComplexHeatmap_mean_none_Focused.png"), width = 1350, height = 6000, res = 450)
    draw(ht, heatmap_legend_side = "top", padding = unit(c(35, 10, 10, 10), "mm"))
    grid.text(paste0(sp_name_clean, ": Mean unscaled counts per sample"), x = 0.5, y = 0.99, gp = gpar(fontsize = 22, fontface = "plain"))

    dev.off()
  }
}

# コドン分類の定義（事前に一度だけ作成）
is_gc3 <- grepl("..[GC]$", codon_order_wo_stop)
is_at3 <- grepl("..[AT]$", codon_order_wo_stop)
is_gt3 <- grepl("..[GT]$", codon_order_wo_stop)
is_ac3 <- grepl("..[AC]$", codon_order_wo_stop)
is_ga3 <- grepl("..[GA]$", codon_order_wo_stop)
is_tc3 <- grepl("..[TC]$", codon_order_wo_stop)
is_atonly <- (nchar(gsub("[AT]", "", codon_order_wo_stop)) == 0) # ATのみ（GCが0個）
is_gcone  <- (nchar(gsub("[AT]", "", codon_order_wo_stop)) == 1) # GCが2個
is_gctwo  <- (nchar(gsub("[AT]", "", codon_order_wo_stop)) == 2) # GCが2個
is_gconly <- (nchar(gsub("[AT]", "", codon_order_wo_stop)) == 3) # GCが3個


# 分類リストを作成
group_list <- list(
  All = rep(TRUE, length(codon_order_wo_stop)),
  GC3 = is_gc3,
  AT3 = is_at3
)

group_list2 <- list(
  All = rep(TRUE, length(codon_order_wo_stop)),
  GT3 = is_gt3,
  AC3 = is_ac3
)

group_list3 <- list(
  All = rep(TRUE, length(codon_order_wo_stop)),
  GA3 = is_ga3,
  TC3 = is_tc3
)

group_list4 <- list(
  All = rep(TRUE, length(codon_order_wo_stop)),
  GC_only = is_gconly,
  GC_two = is_gctwo,
  GC_one = is_gcone,
  AT_only = is_atonly
)

results_list <- list()
ratio_results_list <- list() # ここで初期化

results_list2 <- list()
ratio_results_list2 <- list() # ここで初期化

results_list3 <- list()
ratio_results_list3 <- list() # ここで初期化

stop_pos_map <- c(
  "Homo_sapiens" = 51, "Mus_musculus" = 51, "Arabidopsis_thaliana" = 51, "Oryza_sativa" = 51, 
  "Drosophila_melanogaster" = 51, "Saccharomyces_cerevisiae" = 55, "Escherichia_coli" = 56, "Arabidopsis_thaliana_pt" = 55
) #Based on the heatmaps made as the above


align_data <- function(df, species_name) {
  offset <- stop_pos_map[species_name]
  if(is.na(offset)) offset <- 50 # デフォルト
  df$aligned_pos <- df$pos - offset
  df$species <- species_name
  # 系統ラベルの付与（色のグループ分け用）
  df$lineage <- ifelse(species_name %in% c("Ecoli", "Arabidopsis_pt"), "Prokaryotic-type", "Eukaryotic-type")
  return(df)
}

combined_all_species_data <- list()
combined_all_species_data2 <- list()
combined_all_species_data3 <- list()

combined_list_ratio_gc3 <- list()
combined_list_ratio_gt3 <- list()
combined_list_ratio_ga3 <- list()

for (sp in species_key) {
  results_list <- list()      # GC3/AT3用
  results_list2 <- list()     # GT3/AC3用
  results_list3 <- list()     # GT3/AC3用
  ratio_results_list <- list()
  ratio_results_list2 <- list()
  ratio_results_list3 <- list()
  ratio_results_list4 <- list()

  coeffs <- coeffs_all[grepl(paste0(sp,"_thr0"),rownames(coeffs_all)),codon_order_wo_stop,drop=F]
  coeffs<-t(coeffs)
  if (sp == "Arabidopsis_thaliana_pt"){
      coeffs<-read.table("Arabidopsis_thaliana_all_cm_gene_tpm.mean.linear.coeffs.tsv",row.names=1,header=T,sep="\t")
      coeffs<-coeffs[codon_order_wo_stop,,drop=F]
  }
  # 1. データの読み込みとリスト化
  data_list <- list(
    Original = counts[grepl(sp, counts[, ncol(counts)]) & !grepl("_sub0", counts[, ncol(counts)]) & !grepl("flat", counts[, ncol(counts)]), , drop = F],
    Sub0     = counts[grepl(sp, counts[, ncol(counts)]) & grepl("_sub0", counts[, ncol(counts)]) & !grepl("_sub000", counts[, ncol(counts)]), , drop = F],
    Sub000   = counts[grepl(sp, counts[, ncol(counts)]) & grepl("_sub000", counts[, ncol(counts)]), , drop = F],
    Flat_regions   = counts[grepl(sp, counts[, ncol(counts)]) & grepl("_flat_regions", counts[, ncol(counts)]), , drop = F],
    Flat_genes   = counts[grepl(sp, counts[, ncol(counts)]) & grepl("_flat_genes", counts[, ncol(counts)]), , drop = F]
  )
  if (sp == "Arabidopsis_thaliana"){
      data_list <- list(
        Original = counts[grepl(sp, counts[, ncol(counts)]) & !grepl("_sub0", counts[, ncol(counts)]) & !grepl("flat", counts[, ncol(counts)]) & !grepl("_pt", counts[, ncol(counts)]), , drop = F],
        Sub0     = counts[grepl(sp, counts[, ncol(counts)]) & grepl("_sub0", counts[, ncol(counts)]) & !grepl("_sub000", counts[, ncol(counts)]) & !grepl("_pt", counts[, ncol(counts)]), , drop = F],
        Sub000   = counts[grepl(sp, counts[, ncol(counts)]) & grepl("_sub000", counts[, ncol(counts)]) & !grepl("_pt", counts[, ncol(counts)]), , drop = F],
        Flat_regions   = counts[grepl(sp, counts[, ncol(counts)]) & grepl("_flat_regions", counts[, ncol(counts)]) & !grepl("_pt", counts[, ncol(counts)]), , drop = F],
        Flat_genes   = counts[grepl(sp, counts[, ncol(counts)]) & grepl("_flat_genes", counts[, ncol(counts)]) & !grepl("_pt", counts[, ncol(counts)]), , drop = F]
      )
  } else if (sp == "Arabidopsis_thaliana_pt"){
      sp0<-gsub("_pt","",sp)
      data_list <- list(
        Original = counts[grepl(sp0, counts[, ncol(counts)]) & !grepl("_sub0", counts[, ncol(counts)]) & !grepl("flat", counts[, ncol(counts)]) & grepl("_pt", counts[, ncol(counts)]), , drop = F],
        Sub0     = counts[grepl(sp0, counts[, ncol(counts)]) & grepl("_sub0", counts[, ncol(counts)]) & !grepl("_sub000", counts[, ncol(counts)]) & grepl("_pt", counts[, ncol(counts)]), , drop = F],
        Sub000   = counts[grepl(sp0, counts[, ncol(counts)]) & grepl("_sub000", counts[, ncol(counts)]) & grepl("_pt", counts[, ncol(counts)]), , drop = F],
        Flat_regions   = counts[grepl(sp0, counts[, ncol(counts)]) & grepl("_flat_regions", counts[, ncol(counts)]) & grepl("_pt", counts[, ncol(counts)]), , drop = F],
        Flat_genes   = counts[grepl(sp0, counts[, ncol(counts)]) & grepl("_flat_genes", counts[, ncol(counts)]) & grepl("_pt", counts[, ncol(counts)]), , drop = F]
      )
  }  

  print(lapply(data_list, nrow))
  
  num_codons <- length(codon_order_wo_stop)
  num_units  <- nrow(data_list$Original) / num_codons
  num_pos    <- ncol(data_list$Original) - 2
  
  # 2. 相関計算の実行
  # 「データ種 × コドン群」の組み合わせごとに結果を格納する
  for (d_name in names(data_list)) {
    curr_data <- data_list[[d_name]]
    
    # --- 重要：1ユニットが何行か（シロイヌナズナが64行なら64にする） ---
    unit_size <- 64 
    num_units <- nrow(curr_data) / unit_size # 768 / 64 = 12 ユニットになるはず
    num_pos <- ncol(curr_data) - 2

    # --- 比率計算 (d_name ごとに1回だけ計算すればOK) ---
    ratio_mtrx <- matrix(NA, nrow = num_units, ncol = num_pos)
    for (j in 1:num_units) {
      start_idx <- (j - 1) * unit_size + 1
      end_idx   <- j * unit_size
      u_raw     <- curr_data[start_idx:end_idx, ]
      u_clean   <- u_raw[match(codon_order_wo_stop, u_raw[,1]), 2:(num_pos + 1), drop = F]
      
      u_gc3_avg <- colMeans(u_clean[group_list$GC3, , drop = F], na.rm = TRUE)
      u_at3_avg <- colMeans(u_clean[group_list$AT3, , drop = F], na.rm = TRUE)
      ratio_mtrx[j, ] <- u_gc3_avg / (u_at3_avg + 1e-10)
    }

    m_ratio <- colMeans(ratio_mtrx, na.rm = TRUE)
    s_ratio <- apply(ratio_mtrx, 2, sd, na.rm = TRUE) / sqrt(nrow(ratio_mtrx))
    
    # 比率はデータ種ごとに保存
    ratio_results_list[[d_name]] <- data.frame(
      pos = 1:length(m_ratio), 
      mean = m_ratio, 
      ymin = m_ratio - s_ratio, 
      ymax = m_ratio + s_ratio,
      data_gp = d_name,
      codon_gp = "GC3/AT3"
    )

    ratio_mtrx2 <- matrix(NA, nrow = num_units, ncol = num_pos)
    for (j in 1:num_units) {
      start_idx <- (j - 1) * unit_size + 1
      end_idx   <- j * unit_size
      u_raw     <- curr_data[start_idx:end_idx, ]
      u_clean   <- u_raw[match(codon_order_wo_stop, u_raw[,1]), 2:(num_pos + 1), drop = F]
      
      u_gt3_avg <- colMeans(u_clean[group_list2$GT3, , drop = F], na.rm = TRUE)
      u_ac3_avg <- colMeans(u_clean[group_list2$AC3, , drop = F], na.rm = TRUE)
      
      # 1ユニットごとの比率を格納
      ratio_mtrx2[j, ] <- u_gt3_avg / (u_ac3_avg + 1e-10)
    }

    # --- ユニットごとのループ（j）が終わった後に一括計算する ---
    ratio_mtrx2[is.infinite(ratio_mtrx2)] <- NA
    ratio_mtrx2[is.nan(ratio_mtrx2)]      <- NA

    # 1. 平均値の計算
    m_ratio2 <- colMeans(ratio_mtrx2, na.rm = TRUE)
    
    # 2. 有効なユニット数のカウント（ポジションごと）
    n_valid2 <- apply(ratio_mtrx2, 2, function(x) sum(!is.na(x)))
    
    # 3. 標準偏差および標準誤差の計算
    # ズレ防止：n < 2 の場合は sd を計算せず NA にする
    s_ratio2 <- rep(NA, num_pos)
    valid_pos <- which(n_valid2 >= 2)
    s_ratio2[valid_pos] <- apply(ratio_mtrx2[, valid_pos, drop=F], 2, sd, na.rm = TRUE) / sqrt(n_valid2[valid_pos])

    # 4. データ不足地点の最終クリーニング（平均値も無効化してズレを防止）
    m_ratio2[n_valid2 < 2] <- NA
    
    # 比率はデータ種ごとに保存
    ratio_results_list2[[d_name]] <- data.frame(
      pos = 1:length(m_ratio2), 
      mean = m_ratio2, 
      ymin = m_ratio2 - s_ratio2, 
      ymax = m_ratio2 + s_ratio2,
      data_gp = d_name,
      codon_gp = "GT3/AC3"
    )
    
    ratio_mtrx3 <- matrix(NA, nrow = num_units, ncol = num_pos)
    for (j in 1:num_units) {
      start_idx <- (j - 1) * unit_size + 1
      end_idx   <- j * unit_size
      u_raw     <- curr_data[start_idx:end_idx, ]
      u_clean   <- u_raw[match(codon_order_wo_stop, u_raw[,1]), 2:(num_pos + 1), drop = F]
      
      u_ga3_avg <- colMeans(u_clean[group_list3$GA3, , drop = F], na.rm = TRUE)
      u_tc3_avg <- colMeans(u_clean[group_list3$TC3, , drop = F], na.rm = TRUE)
      
      # 1ユニットごとの比率を格納
      ratio_mtrx3[j, ] <- u_ga3_avg / (u_tc3_avg + 1e-10)
    }

    # --- ユニットごとのループ（j）が終わった後に一括計算する ---
    ratio_mtrx3[is.infinite(ratio_mtrx3)] <- NA
    ratio_mtrx3[is.nan(ratio_mtrx3)]      <- NA

    # 1. 平均値の計算
    m_ratio3 <- colMeans(ratio_mtrx3, na.rm = TRUE)
    
    # 2. 有効なユニット数のカウント（ポジションごと）
    n_valid3 <- apply(ratio_mtrx3, 2, function(x) sum(!is.na(x)))
    
    # 3. 標準偏差および標準誤差の計算
    # ズレ防止：n < 2 の場合は sd を計算せず NA にする
    s_ratio3 <- rep(NA, num_pos)
    valid_pos <- which(n_valid3 >= 2)
    s_ratio3[valid_pos] <- apply(ratio_mtrx3[, valid_pos, drop=F], 2, sd, na.rm = TRUE) / sqrt(n_valid3[valid_pos])

    # 4. データ不足地点の最終クリーニング（平均値も無効化してズレを防止）
    m_ratio3[n_valid3 < 2] <- NA
    
    # 比率はデータ種ごとに保存
    ratio_results_list3[[d_name]] <- data.frame(
      pos = 1:length(m_ratio3), 
      mean = m_ratio3, 
      ymin = m_ratio3 - s_ratio3, 
      ymax = m_ratio3 + s_ratio3,
      data_gp = d_name,
      codon_gp = "GA3/TC3"
    )

    for (g_name in names(group_list)) {
      mask <- group_list[[g_name]]
      emtrx <- matrix(NA, nrow = num_units, ncol = num_pos)
      pmtrx <- matrix(NA, nrow = num_units, ncol = num_pos)

      for (j in 1:num_units) {
        # --- ここで正確に「上から64行」を切り出す ---
        start_idx <- (j - 1) * unit_size + 1
        end_idx   <- j * unit_size
        u_raw     <- curr_data[start_idx:end_idx, ]
        
        # 行名設定時の重複エラーを回避
        # 1列目にコドン名がある前提。もし重複があれば make.unique で回避
        rownames(u_raw) <- make.unique(as.character(u_raw[, 1]))
        
        # 解析に使用するコドン（codon_order_wo_stop）だけに絞り込み
        # これにより 64行の中から 61行だけを正しい順序で取り出せます
        u_clean <- u_raw[codon_order_wo_stop, , drop = F]
        
        # サブセット（GC3/AT3）の抽出
        u_sub <- u_clean[mask, , drop = F]
        c_sub <- coeffs[mask, , drop = F]
        
        if(d_name == "Original"){org_u<-u_clean}

        for (i in 1:num_pos) {
           x <- as.numeric(u_sub[, i + 1])
           y <- as.numeric(unlist(c_sub))
  
           # 標準偏差が0、またはデータが不足している場合
           sx <- sd(x, na.rm = TRUE)
           sy <- sd(y, na.rm = TRUE)

           # 1. いずれかのsdがNA（データ不足や欠損）
           # 2. いずれかのsdが0（全ての値が同じで相関計算不能）
           # のいずれかに当てはまる場合は、計算をスキップして「相関0」とする
           if (is.na(sx) || is.na(sy) || sx == 0 || sy == 0) {
             emtrx[j, i] <- 0
             pmtrx[j, i] <- 1
           } else {
             # ここまで来れば安全に cor.test が実行可能
             res <- cor.test(x, y, use = "complete.obs") 
             emtrx[j, i] <- res$estimate
             pmtrx[j, i] <- res$p.value
           }
         }
      }
            
      # 統計量の算出 (-log10Q * R)
      qmtrx <- matrix(p.adjust(pmtrx, method = "fdr"), nrow = num_units)
      pimtrx <- -log10(qmtrx + 1e-10) * emtrx
      
      # プロット用の集計（平均と標準誤差）
      m <- colMeans(pimtrx, na.rm = TRUE)
      s <- apply(pimtrx, 2, sd, na.rm = TRUE) / sqrt(nrow(pimtrx))
      
      results_list[[paste(d_name, g_name)]] <- data.frame(
        pos = 1:length(m), mean = m, ymin = m - s, ymax = m + s,
        data_gp = d_name, codon_gp = g_name
      )
      rownames(qmtrx) <- 1:num_units
      colnames(qmtrx) <- colnames(curr_data)[2:(num_pos+1)] # ポジション名
      rownames(pimtrx) <- 1:num_units
      colnames(pimtrx) <- colnames(curr_data)[2:(num_pos+1)]
      rownames(emtrx) <- 1:num_units
      colnames(emtrx) <- colnames(curr_data)[2:(num_pos+1)]
      outqname <- paste0(sp, "_", d_name, "_GC3AT3_ribbon_plot_qval.txt")
      outpiname <- paste0(sp, "_", d_name, "_GC3AT3_ribbon_plot_pival.txt")
      outename <- paste0(sp, "_", d_name, "_GC3AT3_ribbon_plot_eval.txt")
      write.table(qmtrx, file = outqname, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
      write.table(pimtrx, file = outpiname, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
      write.table(emtrx, file = outename, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
    }
    
    for (g_name in names(group_list2)) {
      mask <- group_list2[[g_name]]
      emtrx2 <- matrix(NA, nrow = num_units, ncol = num_pos)
      pmtrx2 <- matrix(NA, nrow = num_units, ncol = num_pos)

      for (j in 1:num_units) {
        # --- ここで正確に「上から64行」を切り出す ---
        start_idx <- (j - 1) * unit_size + 1
        end_idx   <- j * unit_size
        u_raw     <- curr_data[start_idx:end_idx, ]
        
        # 行名設定時の重複エラーを回避
        # 1列目にコドン名がある前提。もし重複があれば make.unique で回避
        rownames(u_raw) <- make.unique(as.character(u_raw[, 1]))
        
        # 解析に使用するコドン（codon_order_wo_stop）だけに絞り込み
        # これにより 64行の中から 61行だけを正しい順序で取り出せます
        u_clean <- u_raw[codon_order_wo_stop, , drop = F]
        
        u_sub <- u_clean[mask, , drop = F]
        c_sub <- coeffs[mask, , drop = F]
        
        if(d_name == "Original"){org_u<-u_clean}

        for (i in 1:num_pos) {
           x <- as.numeric(u_sub[, i + 1])
           y <- as.numeric(unlist(c_sub))
  
           # 標準偏差が0、またはデータが不足している場合
           sx <- sd(x, na.rm = TRUE)
           sy <- sd(y, na.rm = TRUE)

           # 1. いずれかのsdがNA（データ不足や欠損）
           # 2. いずれかのsdが0（全ての値が同じで相関計算不能）
           # のいずれかに当てはまる場合は、計算をスキップして「相関0」とする
           if (is.na(sx) || is.na(sy) || sx == 0 || sy == 0) {
             emtrx2[j, i] <- 0
             pmtrx2[j, i] <- 1
           } else {
             # ここまで来れば安全に cor.test が実行可能
             res <- cor.test(x, y, use = "complete.obs") 
             emtrx2[j, i] <- res$estimate
             pmtrx2[j, i] <- res$p.value
           }
         }
      }
            
      # 統計量の算出 (-log10Q * R)
      qmtrx2 <- matrix(p.adjust(pmtrx2, method = "fdr"), nrow = num_units)
      pimtrx2 <- -log10(qmtrx2 + 1e-10) * emtrx2
      
      # プロット用の集計（平均と標準誤差）
      m <- colMeans(pimtrx2, na.rm = TRUE)
      s <- apply(pimtrx2, 2, sd, na.rm = TRUE) / sqrt(nrow(pimtrx2))
      
      results_list2[[paste(d_name, g_name)]] <- data.frame(
        pos = 1:length(m), mean = m, ymin = m - s, ymax = m + s,
        data_gp = d_name, codon_gp = g_name
      )
      rownames(qmtrx2) <- 1:num_units
      colnames(qmtrx2) <- colnames(curr_data)[2:(num_pos+1)] # ポジション名
      rownames(pimtrx2) <- 1:num_units
      colnames(pimtrx2) <- colnames(curr_data)[2:(num_pos+1)]
      rownames(emtrx2) <- 1:num_units
      colnames(emtrx2) <- colnames(curr_data)[2:(num_pos+1)]
      outqname<-paste0(sp, d_name, "_GT3AC3_ribbon_plot_qval.txt")
      outpiname<-paste0(sp, d_name, "_GT3AC3_ribbon_plot_pival.txt")
      outename<-paste0(sp, d_name, "_GT3AC3_ribbon_plot_eval.txt")
      write.table(qmtrx2, file = outqname, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
      write.table(pimtrx2, file = outpiname, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
      write.table(emtrx2, file = outename, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
    }

    for (g_name in names(group_list3)) {
      mask <- group_list3[[g_name]]
      emtrx3 <- matrix(NA, nrow = num_units, ncol = num_pos)
      pmtrx3 <- matrix(NA, nrow = num_units, ncol = num_pos)

      for (j in 1:num_units) {
        # --- ここで正確に「上から64行」を切り出す ---
        start_idx <- (j - 1) * unit_size + 1
        end_idx   <- j * unit_size
        u_raw     <- curr_data[start_idx:end_idx, ]
        
        # 行名設定時の重複エラーを回避
        # 1列目にコドン名がある前提。もし重複があれば make.unique で回避
        rownames(u_raw) <- make.unique(as.character(u_raw[, 1]))
        
        # 解析に使用するコドン（codon_order_wo_stop）だけに絞り込み
        # これにより 64行の中から 61行だけを正しい順序で取り出せます
        u_clean <- u_raw[codon_order_wo_stop, , drop = F]
        
        u_sub <- u_clean[mask, , drop = F]
        c_sub <- coeffs[mask, , drop = F]
        
        if(d_name == "Original"){org_u<-u_clean}

        for (i in 1:num_pos) {
           x <- as.numeric(u_sub[, i + 1])
           y <- as.numeric(unlist(c_sub))
  
           # 標準偏差が0、またはデータが不足している場合
           sx <- sd(x, na.rm = TRUE)
           sy <- sd(y, na.rm = TRUE)

           # 1. いずれかのsdがNA（データ不足や欠損）
           # 2. いずれかのsdが0（全ての値が同じで相関計算不能）
           # のいずれかに当てはまる場合は、計算をスキップして「相関0」とする
           if (is.na(sx) || is.na(sy) || sx == 0 || sy == 0) {
             emtrx3[j, i] <- 0
             pmtrx3[j, i] <- 1
           } else {
             # ここまで来れば安全に cor.test が実行可能
             res <- cor.test(x, y, use = "complete.obs") 
             emtrx3[j, i] <- res$estimate
             pmtrx3[j, i] <- res$p.value
           }
         }
      }
            
      # 統計量の算出 (-log10Q * R)
      qmtrx3 <- matrix(p.adjust(pmtrx3, method = "fdr"), nrow = num_units)
      pimtrx3 <- -log10(qmtrx3 + 1e-10) * emtrx3
      
      # プロット用の集計（平均と標準誤差）
      m <- colMeans(pimtrx3, na.rm = TRUE)
      s <- apply(pimtrx3, 2, sd, na.rm = TRUE) / sqrt(nrow(pimtrx3))
      
      results_list3[[paste(d_name, g_name)]] <- data.frame(
        pos = 1:length(m), mean = m, ymin = m - s, ymax = m + s,
        data_gp = d_name, codon_gp = g_name
      )
      rownames(qmtrx3) <- 1:num_units
      colnames(qmtrx3) <- colnames(curr_data)[2:(num_pos+1)] # ポジション名
      rownames(pimtrx3) <- 1:num_units
      colnames(pimtrx3) <- colnames(curr_data)[2:(num_pos+1)]
      rownames(emtrx3) <- 1:num_units
      colnames(emtrx3) <- colnames(curr_data)[2:(num_pos+1)]
      outqname<-paste0(sp, d_name,"_GA3TC3_ribbon_plot_qval.txt")
      outpiname<-paste0(sp, d_name,"_GA3TC3_ribbon_plot_pival.txt")
      outename<-paste0(sp, d_name,"_GA3TC3_ribbon_plot_eval.txt")
      write.table(qmtrx3, file = outqname, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
      write.table(pimtrx3, file = outpiname, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
      write.table(emtrx3, file = outename, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
    }
  }
  
  # 3. データの結合と可視化
  plot_data <- do.call(rbind, results_list)
  results_list$codon_gp<-factor(results_list$codon_gp, levels=names(group_list))
  plot_data2 <- do.call(rbind, results_list2)
  results_list2$codon_gp<-factor(results_list2$codon_gp, levels=names(group_list2))
  plot_data3 <- do.call(rbind, results_list3)
  results_list3$codon_gp<-factor(results_list3$codon_gp, levels=names(group_list3))

  # 「All」を除外して GC3/AT3 の並走を見たい場合はこちらでフィルタ
  # plot_data <- plot_data[plot_data$codon_gp != "All", ]

  rp <- ggplot(plot_data, aes(x = pos, y = mean, 
                             color = codon_gp, 
                             fill = codon_gp, 
                             linetype = data_gp, 
                             group = interaction(data_gp, codon_gp))) + # 15本を独立させる
    geom_hline(yintercept = 0, color = "gray80") +
    geom_ribbon(aes(ymin = ymin, ymax = ymax, alpha = data_gp), color = NA) +
    geom_line(linewidth = 0.8, alpha=0.5) + # linewidth に修正
    scale_alpha_manual(values = c(
      "Original" = 0.4, 
      "Sub0" = 0.2, 
      "Sub000" = 0.2, 
      "Flat_regions" = 0.3, # 視認性のために少し変える
      "Flat_genes" = 0.2    # 視認性のために少し変える
    )) +
    scale_linetype_manual(values = c(
      "Original" = "solid", 
      "Sub0" = "dashed", 
      "Sub000" = "dotted", 
      "Flat_regions" = "dotdash", 
      "Flat_genes" = "longdash"
    )) +
    scale_color_manual(values = c("All" = "black", "GC3" = "royalblue", "AT3" = "darkorange")) +
    scale_fill_manual(values = c("All" = "darkcyan", "GC3" = "royalblue", "AT3" = "darkorange")) +
    theme_bw() + # theme_minimalよりガイド線がはっきりするtheme_bwを選択
    theme(
      # 軸の数値（目盛りラベル）を大きく太くする
      axis.text.x = element_text(size = 14, face = "plain", color = "black"),
      axis.text.y = element_text(size = 14, face = "plain", color = "black"),
      
      # 軸のタイトル（Position, -log10Q...）を大きくする
      axis.title.x = element_text(size = 16, face = "plain"),
      axis.title.y = element_text(size = 16, face = "plain"),
      
      # 図のタイトルを大きくする
      plot.title = element_text(size = 18, face = "bold"),
      
      # 凡例の文字も調整
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 13, face = "plain"),
      
      # ガイド線（グリッド線）の調整
      panel.grid.major.x = element_line(color = "gray90", linewidth = 0.5), # 縦の主要線
      panel.grid.minor.x = element_blank(),                                # 縦の補助線（不要なら消す）
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.5)  # 横の主要線
    ) +
    labs(title = paste("Species:", sp), x = "Position", y = expression(-log[10](Q) %*% coefficient))
    head(results_list[["Flat_regions GC3"]]$mean)
    head(results_list[["Original GC3"]]$mean)
  print(rp)
  outplotname<-paste0(sp,"_GC3AT3_ribbon_plot.png")
  ggsave(outplotname, plot=rp, width=8, height=5, dpi=450, limitsize = FALSE)

  rp2 <- ggplot(plot_data2, aes(x = pos, y = mean, 
                             color = codon_gp, 
                             fill = codon_gp, 
                             linetype = data_gp, 
                             group = interaction(data_gp, codon_gp))) + # 15本を独立させる
    geom_hline(yintercept = 0, color = "gray80") +
    geom_ribbon(aes(ymin = ymin, ymax = ymax, alpha = data_gp), color = NA) +
    geom_line(linewidth = 0.8, alpha=0.5) + # linewidth に修正
    scale_alpha_manual(values = c(
      "Original" = 0.4, 
      "Sub0" = 0.2, 
      "Sub000" = 0.2, 
      "Flat_regions" = 0.3, # 視認性のために少し変える
      "Flat_genes" = 0.2    # 視認性のために少し変える
    )) +
    scale_linetype_manual(values = c(
      "Original" = "solid", 
      "Sub0" = "dashed", 
      "Sub000" = "dotted", 
      "Flat_regions" = "dotdash", 
      "Flat_genes" = "longdash"
    )) +
    scale_color_manual(values = c("All" = "black", "GT3" = "royalblue", "AC3" = "darkorange")) +
    scale_fill_manual(values = c("All" = "darkcyan", "GT3" = "royalblue", "AC3" = "darkorange")) +
    theme_bw() + # theme_minimalよりガイド線がはっきりするtheme_bwを選択
    theme(
      # 軸の数値（目盛りラベル）を大きく太くする
      axis.text.x = element_text(size = 14, face = "plain", color = "black"),
      axis.text.y = element_text(size = 14, face = "plain", color = "black"),
      
      # 軸のタイトル（Position, -log10Q...）を大きくする
      axis.title.x = element_text(size = 16, face = "plain"),
      axis.title.y = element_text(size = 16, face = "plain"),
      
      # 図のタイトルを大きくする
      plot.title = element_text(size = 18, face = "bold"),
      
      # 凡例の文字も調整
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 13, face = "plain"),
      
      # ガイド線（グリッド線）の調整
      panel.grid.major.x = element_line(color = "gray90", linewidth = 0.5), # 縦の主要線
      panel.grid.minor.x = element_blank(),                                # 縦の補助線（不要なら消す）
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.5)  # 横の主要線
    ) +
    labs(title = paste("Species:", sp), x = "Position", y = expression(-log[10](Q) %*% coefficient))
    head(results_list[["Flat_regions GT3"]]$mean)
    head(results_list[["Original AC3"]]$mean)
  print(rp2)

  outplotname<-paste0(sp,"_GT3AC3_ribbon_plot.png")
  ggsave(outplotname, plot=rp2, width=8, height=5, dpi=450, limitsize = FALSE)

  rp3 <- ggplot(plot_data3, aes(x = pos, y = mean, 
                             color = codon_gp, 
                             fill = codon_gp, 
                             linetype = data_gp, 
                             group = interaction(data_gp, codon_gp))) + # 15本を独立させる
    geom_hline(yintercept = 0, color = "gray80") +
    geom_ribbon(aes(ymin = ymin, ymax = ymax, alpha = data_gp), color = NA) +
    geom_line(linewidth = 0.8, alpha=0.5) + # linewidth に修正
    scale_alpha_manual(values = c(
      "Original" = 0.4, 
      "Sub0" = 0.2, 
      "Sub000" = 0.2, 
      "Flat_regions" = 0.3, # 視認性のために少し変える
      "Flat_genes" = 0.2    # 視認性のために少し変える
    )) +
    scale_linetype_manual(values = c(
      "Original" = "solid", 
      "Sub0" = "dashed", 
      "Sub000" = "dotted", 
      "Flat_regions" = "dotdash", 
      "Flat_genes" = "longdash"
    )) +
    scale_color_manual(values = c("All" = "black", "GA3" = "royalblue", "TC3" = "darkorange")) +
    scale_fill_manual(values = c("All" = "darkcyan", "GA3" = "royalblue", "TC3" = "darkorange")) +
    theme_bw() + # theme_minimalよりガイド線がはっきりするtheme_bwを選択
    theme(
      # 軸の数値（目盛りラベル）を大きく太くする
      axis.text.x = element_text(size = 14, face = "plain", color = "black"),
      axis.text.y = element_text(size = 14, face = "plain", color = "black"),
      
      # 軸のタイトル（Position, -log10Q...）を大きくする
      axis.title.x = element_text(size = 16, face = "plain"),
      axis.title.y = element_text(size = 16, face = "plain"),
      
      # 図のタイトルを大きくする
      plot.title = element_text(size = 18, face = "bold"),
      
      # 凡例の文字も調整
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 13, face = "plain"),
      
      # ガイド線（グリッド線）の調整
      panel.grid.major.x = element_line(color = "gray90", linewidth = 0.5), # 縦の主要線
      panel.grid.minor.x = element_blank(),                                # 縦の補助線（不要なら消す）
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.5)  # 横の主要線
    ) +
    labs(title = paste("Species:", sp), x = "Position", y = expression(-log[10](Q) %*% coefficient))
    head(results_list[["Flat_regions GA3"]]$mean)
    head(results_list[["Original TC3"]]$mean)
  print(rp3)

  outplotname<-paste0(sp,"_GA3TC3_ribbon_plot.png")
  ggsave(outplotname, plot=rp3, width=8, height=5, dpi=450, limitsize = FALSE)

  # Making ratio
  display_names <- c(
    "Original"     = "Original",
    "Sub0"         = "Stalling+",
    "Sub000"       = "Stalling++",
    "Flat_regions" = "Stalling+\nsubtracted",
    "Flat_genes"   = "Genes with\nno stalling+"
  )

  # データの統合
  plot_ratio_data <- do.call(rbind, ratio_results_list)
  plot_ratio_data <- plot_ratio_data[order(plot_ratio_data$data_gp, plot_ratio_data$pos), ]
  
  plot_ratio_data <- plot_ratio_data %>% 
    dplyr::filter(mean < 5 & mean > 0)
  
  plot_ratio_data$data_gp <- factor(
    plot_ratio_data$data_gp, 
    levels = names(display_names) # names() で "Original", "Sub0" ... の順を取り出す
  )

  # プロット作成
  ratio_p <- ggplot(plot_ratio_data, aes(x = pos, y = mean, color = data_gp, fill = data_gp, linetype = data_gp, group = data_gp)) +
    geom_hline(yintercept = 1, color = "black", linetype = "dashed") +
    geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2, color = NA) +
    geom_line(linewidth = 1) +
    
    # 凡例のラベルを適用（color, fill, linetypeすべて一致させる）
    scale_linetype_manual(
      values = c("Original"="solid", "Sub0"="dashed", "Sub000"="dotted", "Flat_regions"="dotdash", "Flat_genes"="longdash"),
      labels = display_names
    ) +
    scale_color_discrete(labels = display_names) +
    scale_fill_discrete(labels = display_names) +
    
    theme_bw() +
    theme(
      axis.text = element_text(size = 14, face = "plain"),
      axis.title = element_text(size = 16, face = "plain"),
      
      # --- 凡例の見栄え調整 ---
      legend.title = element_text(size = 14, face = "plain"),
      legend.text = element_text(size = 14, lineheight = 0.9), # 行間を少し詰める
      legend.key.height = unit(2.5, "lines"),                  # 各項目の高さを広げる
      legend.spacing.y = unit(0.5, "cm")                       # 項目間の隙間を追加
    ) +
    # 垂直方向の凡例配置を最適化
    guides(
      color = guide_legend(byrow = TRUE),
      fill = guide_legend(byrow = TRUE),
      linetype = guide_legend(byrow = TRUE)
    ) +
    labs(
      title = paste("Species:", gsub("_", " ", sp)),
      x = "Position",
      y = "Ratio (GC3 / AT3)",
      color = "Data catogory", linetype = "Data catogory", fill = "Data catogory"
    )
  ggsave(paste0(sp, "_GC3AT3_ratio_ribbon.png"), plot = ratio_p, width = 8, height = 5, dpi = 450)
  
  plot_ratio_data2 <- do.call(rbind, ratio_results_list2)
  plot_ratio_data2 <- plot_ratio_data2[order(plot_ratio_data2$data_gp, plot_ratio_data2$pos), ]

  plot_ratio_data2 <- plot_ratio_data2 %>% 
    dplyr::filter(mean < 5 & mean > 0)

  plot_ratio_data2$data_gp <- factor(
    plot_ratio_data2$data_gp, 
    levels = names(display_names) # names() で "Original", "Sub0" ... の順を取り出す
  )

  ratio_p2 <- ggplot(plot_ratio_data2, aes(x = pos, y = mean, color = data_gp, fill = data_gp, linetype = data_gp, group = data_gp)) +
    geom_hline(yintercept = 1, color = "black", linetype = "dashed") +
    geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2, color = NA) +
    geom_line(linewidth = 1) +
    
    # 凡例のラベルを適用（color, fill, linetypeすべて一致させる）
    scale_linetype_manual(
      values = c("Original"="solid", "Sub0"="dashed", "Sub000"="dotted", "Flat_regions"="dotdash", "Flat_genes"="longdash"),
      labels = display_names
    ) +
    scale_color_discrete(labels = display_names) +
    scale_fill_discrete(labels = display_names) +
    
    theme_bw() +
    theme(
      axis.text = element_text(size = 14, face = "plain"),
      axis.title = element_text(size = 16, face = "plain"),
      
      # --- 凡例の見栄え調整 ---
      legend.title = element_text(size = 14, face = "plain"),
      legend.text = element_text(size = 14, lineheight = 0.9), # 行間を少し詰める
      legend.key.height = unit(2.5, "lines"),                  # 各項目の高さを広げる
      legend.spacing.y = unit(0.5, "cm")                       # 項目間の隙間を追加
    ) +
    # 垂直方向の凡例配置を最適化
    guides(
      color = guide_legend(byrow = TRUE),
      fill = guide_legend(byrow = TRUE),
      linetype = guide_legend(byrow = TRUE)
    ) +
    labs(
      title = paste("Species:", gsub("_", " ", sp)),
      x = "Position",
      y = "Ratio (GT3 / AC3)",
      color = "Data catogory", linetype = "Data catogory", fill = "Data catogory"
    )
  ggsave(paste0(sp, "_GT3AC3_ratio_ribbon.png"), plot = ratio_p2, width = 8, height = 5, dpi = 450)

  plot_ratio_data3 <- do.call(rbind, ratio_results_list3)
  plot_ratio_data3 <- plot_ratio_data3[order(plot_ratio_data3$data_gp, plot_ratio_data3$pos), ]

  plot_ratio_data3 <- plot_ratio_data3 %>% 
    dplyr::filter(mean < 5 & mean > 0)

  plot_ratio_data3$data_gp <- factor(
    plot_ratio_data3$data_gp, 
    levels = names(display_names) # names() で "Original", "Sub0" ... の順を取り出す
  )

  ratio_p3 <- ggplot(plot_ratio_data3, aes(x = pos, y = mean, color = data_gp, fill = data_gp, linetype = data_gp, group = data_gp)) +
    geom_hline(yintercept = 1, color = "black", linetype = "dashed") +
    geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2, color = NA) +
    geom_line(linewidth = 1) +
    
    # 凡例のラベルを適用（color, fill, linetypeすべて一致させる）
    scale_linetype_manual(
      values = c("Original"="solid", "Sub0"="dashed", "Sub000"="dotted", "Flat_regions"="dotdash", "Flat_genes"="longdash"),
      labels = display_names
    ) +
    scale_color_discrete(labels = display_names) +
    scale_fill_discrete(labels = display_names) +
    
    theme_bw() +
    theme(
      axis.text = element_text(size = 14, face = "plain"),
      axis.title = element_text(size = 16, face = "plain"),
      
      # --- 凡例の見栄え調整 ---
      legend.title = element_text(size = 14, face = "plain"),
      legend.text = element_text(size = 14, lineheight = 0.9), # 行間を少し詰める
      legend.key.height = unit(2.5, "lines"),                  # 各項目の高さを広げる
      legend.spacing.y = unit(0.5, "cm")                       # 項目間の隙間を追加
    ) +
    # 垂直方向の凡例配置を最適化
    guides(
      color = guide_legend(byrow = TRUE),
      fill = guide_legend(byrow = TRUE),
      linetype = guide_legend(byrow = TRUE)
    ) +
    labs(
      title = paste("Species:", gsub("_", " ", sp)),
      x = "Position",
      y = "Ratio (GA3 / TC3)",
      color = "Data catogory", linetype = "Data catogory", fill = "Data catogory"
    )
  ggsave(paste0(sp, "_GA3TC3_ratio_ribbon.png"), plot = ratio_p3, width = 8, height = 4, dpi = 450)

  rownames(ratio_mtrx) <- 1:num_units
  colnames(ratio_mtrx) <- colnames(curr_data)[2:(num_pos+1)]
  rownames(ratio_mtrx2) <- 1:num_units
  colnames(ratio_mtrx2) <- colnames(curr_data)[2:(num_pos+1)]
  rownames(ratio_mtrx3) <- 1:num_units
  colnames(ratio_mtrx3) <- colnames(curr_data)[2:(num_pos+1)]
  outrationame<-paste0(sp,"_GC3AT3_ratio.txt")
  outratio2name<-paste0(sp,"_GT3AC3_ratio.txt")
  outratio3name<-paste0(sp,"_GA3TC3_ratio.txt")

  write.table(ratio_mtrx, file = outrationame, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
  write.table(ratio_mtrx2, file = outratio2name, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
  write.table(ratio_mtrx3, file = outratio3name, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

  combined_list_ratio_gc3[[sp]] <- align_data(do.call(rbind, ratio_results_list), sp)
  combined_list_ratio_gt3[[sp]] <- align_data(do.call(rbind, ratio_results_list2), sp)
  combined_list_ratio_ga3[[sp]] <- align_data(do.call(rbind, ratio_results_list3), sp)

  combined_all_species_data[[sp]] <- plot_data %>% mutate(species = sp)
  combined_all_species_data2[[sp]] <- plot_data2 %>% mutate(species = sp)
  combined_all_species_data3[[sp]] <- plot_data3 %>% mutate(species = sp)
}

#To make plots with multiple species
plot_all_gc3 <- do.call(rbind, combined_list_ratio_gc3)
plot_all_gt3 <- do.call(rbind, combined_list_ratio_gt3)
plot_all_ga3 <- do.call(rbind, combined_list_ratio_ga3)

plot_all_gc3 <- plot_all_gc3 %>% 
  dplyr::filter(mean < 5 & mean > 0)
summary(plot_all_gc3$mean)

plot_all_gt3 <- plot_all_gt3 %>% 
  dplyr::filter(mean < 5 & mean > 0)
summary(plot_all_gt3$mean)

plot_all_ga3 <- plot_all_ga3 %>% 
  dplyr::filter(mean < 5 & mean > 0)
summary(plot_all_ga3$mean)

plot_sub_gc3 <- plot_all_gc3[grepl("Sub0",plot_all_gc3$data_gp),]
plot_sub_gt3 <- plot_all_gt3[grepl("Sub0",plot_all_gt3$data_gp),]
plot_sub_ga3 <- plot_all_ga3[grepl("Sub0",plot_all_ga3$data_gp),]

stop_pos_map <- c(
  "Homo_sapiens" = 51, "Mus_musculus" = 51, "Arabidopsis_thaliana" = 51, "Oryza_sativa" = 51, 
  "Drosophila_melanogaster" = 51, "Saccharomyces_cerevisiae" = 55, "Escherichia_coli" = 56, "Arabidopsis_thaliana_pt" = 55
) #This will set the P site as the position 0

species_labels <- c(
  "Arabidopsis_thaliana" = "Arabidopsis\nthaliana",
  "Oryza_sativa" = "Oryza\nsativa",
  "Homo_sapiens" = "Homo\nsapiens",
  "Mus_musculus" = "Mus\nmusculus",
  "Drosophila_melanogaster" = "Drosophila\nmelanogaster",
  "Saccharomyces_cerevisiae" = "Saccharomyces\ncerevisiae",
  "Escherichia_coli" = "Escherichia\ncoli",
  "Arabidopsis_thaliana_pt" = "Arabidopsis\nthaliana\nplastid"
)

# 系統別カラー設定
lineage_colors <- c(
  # --- 哺乳類 (哺乳綱) ---
  "Homo_sapiens"            = "#E41A1C", # 赤
  "Mus_musculus"            = "#FF7F00", # オレンジ
  
  # --- 被子植物 (真正双子葉類・単子葉類) ---
  "Arabidopsis_thaliana"    = "#4DAF4A", # 緑
  "Oryza_sativa"            = "#006400", # 明るい黄緑
  "Arabidopsis_thaliana_pt" = "#94D050", # 深緑 (細胞内共生由来を表現)
  
  # --- その他の真核生物 (後生動物・真菌) ---
  "Drosophila_melanogaster" = "#984EA3", # 紫
  "Saccharomyces_cerevisiae" = "#F781BF", # ピンク
  
  # --- 原核生物 (真正細菌) ---
  "Escherichia_coli"        = "#377EB8"  # 青
)


plot_total_combined <- bind_rows(
  plot_sub_gc3 %>% mutate(comparison_type = "GC3 / AT3"),
  plot_sub_gt3 %>% mutate(comparison_type = "GT3 / AC3"),
  plot_sub_ga3 %>% mutate(comparison_type = "GA3 / TC3")
)

plot_total_combined$species <- factor(plot_total_combined$species, levels = names(species_labels))

# 表示順を固定したい場合
plot_total_combined$comparison_type <- factor(
  plot_total_combined$comparison_type, 
  levels = c("GC3 / AT3", "GT3 / AC3", "GA3 / TC3")
)

# --- 2. 統合プロット用関数の作成 ---
plot_total_with_all <- plot_total_combined %>%
  # aligned_posを計算済みであることを確認
  mutate(
    species_eds = species %>% 
      stringr::str_replace("_pt$", " plastid") %>% 
      str_replace_all("_", " ")
  )

# 全種をまとめたデータを作成
# --- 1. データの準備（色分け用IDを保持しつつ統合） ---
all_overlay_data <- plot_total_combined %>%
  mutate(
    color_id = species,        # 元の種名を色用に保存
    species = "All_species"    # パネル分け用を書き換え
  )

plot_final_for_facet <- bind_rows(
  plot_total_combined %>% mutate(color_id = species), 
  all_overlay_data
)

# 順序の固定
ordered_with_all <- c(names(species_labels), "All_species")
plot_final_for_facet$species <- factor(plot_final_for_facet$species, levels = ordered_with_all)
species_labels_with_all <- c(species_labels, "All_species" = "All species")

lineage_colors_with_all <- lineage_colors

# --- 2. 関数の修正（列名を一致させる） ---
plot_universal_combined_facet_v2 <- function(plot_df) {
  
  # データの集約
  plot_summary <- plot_df %>%
    # facetに使う名前を 'species' のまま固定し、色用の 'color_id' もグループに含める
    dplyr::group_by(species, color_id, comparison_type, aligned_pos) %>%
    dplyr::summarise(
      mean = mean(mean, na.rm = TRUE),
      ymin = mean(ymin, na.rm = TRUE),
      ymax = mean(ymax, na.rm = TRUE),
      .groups = "drop"
    )
    
  ggplot(plot_summary, aes(x = aligned_pos, y = mean, 
                           color = color_id,  # 個別の種名で色分け
                           fill = color_id,   # 個別の種名で塗りつぶし
                           group = color_id)) + # 重ね書きパネル内で線を分離
    
    # 金色の影
    annotate("rect", xmin = -5, xmax = 5, ymin = -Inf, ymax = Inf, alpha = 0.3, fill = "gold") + 
    geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
    
    # リボンと線
    geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.1, color = NA) +
    geom_line(linewidth = 0.6, alpha = 0.7) + 
    
    # カラー適用（All_species用のダミーを追加したリストを使用）
    scale_color_manual(values = lineage_colors_with_all) +
    scale_fill_manual(values = lineage_colors_with_all) +
    
    # facet_grid の変数名をデータ列名の 'species' と一致させる
    facet_grid(species ~ comparison_type, scales = "free_y", 
               labeller = labeller(species = species_labels_with_all)) +
    
    coord_cartesian(xlim = c(-15, 15)) + 
    theme_bw() +
    theme(
      strip.background = element_blank(),
      strip.text.y = element_text(face = "italic", hjust = 0, angle = 0, size = 14, lineheight = 0.8),
      strip.text.x = element_text(size = 14, face = "plain"),
      axis.text = element_text(size = 14, color = "black"),
      axis.title = element_text(size = 16),
      legend.position = "none",
      panel.spacing = unit(0.3, "lines")
    ) +
    labs(x = "Codon position (0 = P site)", y = "Ratio")
}

# --- 3. 実行 ---
combined_rp_v2 <- plot_universal_combined_facet_v2(plot_final_for_facet)
ggsave("all_ratios_with_overlay.png", plot = combined_rp_v2, width = 8, height = 12, dpi = 450)



plot_universal_comparison <- function(plot_df, target_groups, title_suffix) {
  
  # 1. データの整理（ここで同じ位置のデータを1つにまとめる）
  plot_sub <- plot_df %>%
    dplyr::filter(codon_gp %in% target_groups) %>%
    dplyr::group_by(species, codon_gp, aligned_pos) %>%
    dplyr::summarise(
      mean = mean(mean, na.rm = TRUE),
      ymin = mean(ymin, na.rm = TRUE),
      ymax = mean(ymax, na.rm = TRUE),
      .groups = "drop"
    )
    
  # 2. 可視化
  ggplot(plot_sub, aes(x = aligned_pos, y = mean, 
                       color = species, 
                       fill = species,
                       linetype = species,
                       group = species)) + # speciesごとに線を引く
    
    # 基準線 (Ratio=1)
    geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +

    annotate("rect", xmin = -5, xmax = 5, ymin = -Inf, ymax = Inf, 
             alpha = 0.1, fill = "orange") + 
    
    # メインの描画：リボンと線
    geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.8, alpha=0.6) + 
    
    # スケール設定（色が合わない場合はここを scale_color_discrete() 等に変更）
    # scale_color_manual(values = lineage_colors) + 
    # scale_fill_manual(values = lineage_colors) +
    
    # x軸の範囲とy軸のズーム（平坦に見えるならylimを狭くする）
    coord_cartesian(xlim = c(-25, 25)) + 
    
    theme_bw() +
    facet_wrap(~codon_gp, ncol = 1) +
    labs(title = title_suffix, x = "Position", y = "Ratio")
}

# --- 出力実行 ---
# 1. GC3 / AT3
ggsave("GC3_AT3_species_comparison.png", 
       plot_universal_comparison(plot_sub_gc3, "GC3/AT3", "GC3/AT3 ratio"), 
       width = 10, height = 7)

# 2. GT3 / AC3
ggsave("GT3_AC3_species_comparison.png", 
       plot_universal_comparison(plot_sub_gt3, "GT3/AC3", "GT3/AC3 ratio"), 
       width = 10, height = 7)

# 3. GA3 / TC3
ggsave("GA3_TC3_species_comparison.png", 
       plot_universal_comparison(plot_sub_ga3,  "GA3/TC3", "GA3/TC3 ratio"), 
       width = 10, height = 7)

#############################
#For ribbon plot with multiple species vertically aligned

full_plot_df <- bind_rows(combined_all_species_data) #After running the above big loop
full_plot_df2 <- bind_rows(combined_all_species_data2) #After running the above big loop
full_plot_df3 <- bind_rows(combined_all_species_data3) #After running the above big loop

full_plot_df <- full_plot_df %>%
  mutate(
    # _pt を (plastid) に置換し、アンダースコアを半角スペースに置換
    species_eds = species %>% 
      str_replace("_pt$", " plastid") %>%  # 末尾の _pt を置換
      str_replace_all("_", " ")            # 残りのアンダースコアをスペースに
  )

full_plot_df2 <- full_plot_df2 %>%
  mutate(
    # _pt を (plastid) に置換し、アンダースコアを半角スペースに置換
    species_eds = species %>% 
      str_replace("_pt$", " plastid") %>%  # 末尾の _pt を置換
      str_replace_all("_", " ")            # 残りのアンダースコアをスペースに
  )

full_plot_df3 <- full_plot_df3 %>%
  mutate(
    # _pt を (plastid) に置換し、アンダースコアを半角スペースに置換
    species_eds = species %>% 
      str_replace("_pt$", " plastid") %>%  # 末尾の _pt を置換
      str_replace_all("_", " ")            # 残りのアンダースコアをスペースに
  )

ordered_species_eds <- names(species_labels) %>%
  str_replace("_pt$", " plastid") %>%
  str_replace_all("_", " ")

full_plot_df$species <- factor(full_plot_df$species_eds, levels = ordered_species_eds)
full_plot_df2$species <- factor(full_plot_df2$species_eds, levels = ordered_species_eds)
full_plot_df3$species <- factor(full_plot_df3$species_eds, levels = ordered_species_eds)

data_category_map <- c(
  "Original"     = "Original",
  "Sub0"         = "Stalling+",
  "Sub000"       = "Stalling++",
  "Flat_regions" = "Stalling+\nsubtracted",
  "Flat_genes"   = "Genes with\nno stalling+"
)

# 2. 新しい列の追加と順番の固定
full_plot_df <- full_plot_df %>%
  mutate(
    # 名前を置換
    data_gp_eds = data_category_map[data_gp],
    # 順番を固定（ここに書いた順に凡例が並びます）
    data_gp_eds = factor(data_gp_eds, levels = c(
      "Original",
      "Stalling+",
      "Stalling++",
      "Stalling+\nsubtracted",
      "Genes with\nno stalling+"
    ))
  )

full_plot_df2 <- full_plot_df2 %>%
  mutate(
    # 名前を置換
    data_gp_eds = data_category_map[data_gp],
    # 順番を固定（ここに書いた順に凡例が並びます）
    data_gp_eds = factor(data_gp_eds, levels = c(
      "Original",
      "Stalling+",
      "Stalling++",
      "Stalling+\nsubtracted",
      "Genes with\nno stalling+"
    ))
  )

full_plot_df3 <- full_plot_df3 %>%
  mutate(
    # 名前を置換
    data_gp_eds = data_category_map[data_gp],
    # 順番を固定（ここに書いた順に凡例が並びます）
    data_gp_eds = factor(data_gp_eds, levels = c(
      "Original",
      "Stalling+",
      "Stalling++",
      "Stalling+\nsubtracted",
      "Genes with\nno stalling+"
    ))
  )

full_plot_df <- full_plot_df %>%
  mutate(
    # species 列をキーにして stop_pos_map から値を取得
    offset = stop_pos_map[species],
    # 元の位置 pos からオフセットを引いて aligned_pos を作成
    aligned_pos = pos - offset
  )

full_plot_df2 <- full_plot_df2 %>%
  mutate(
    # species 列をキーにして stop_pos_map から値を取得
    offset = stop_pos_map[species],
    # 元の位置 pos からオフセットを引いて aligned_pos を作成
    aligned_pos = pos - offset
  )

full_plot_df3 <- full_plot_df3 %>%
  mutate(
    # species 列をキーにして stop_pos_map から値を取得
    offset = stop_pos_map[species],
    # 元の位置 pos からオフセットを引いて aligned_pos を作成
    aligned_pos = pos - offset
  )


final_facet_plot <- ggplot(full_plot_df, aes(x = aligned_pos, y = mean, 
                                             color = codon_gp, 
                                             fill = codon_gp,
                                             linetype = data_gp_eds, 
                                             alpha = data_gp_eds,
                                             group = interaction(data_gp_eds, codon_gp))) +
  # 【影の追加】 リボソーム被覆部位などの範囲指定（例：-15から0まで）
  # geom_lineより前に書くことで最背面に配置されます
  annotate("rect", xmin = -5, xmax = 5, ymin = -Inf, ymax = Inf, 
           alpha = 0.3, fill = "gold") + 
  
  # 基準線
  geom_hline(yintercept = 0, color = "black", alpha = 0.3) +
  
  # リボンと線
  geom_ribbon(aes(ymin = ymin, ymax = ymax, alpha = data_gp_eds), color = NA) +
  geom_line(linewidth = 0.7, alpha=0.5) +
  
  # スケール設定（既存の設定を維持）
  scale_alpha_manual(values = c("Original" = 0.4, "Stalling+" = 0.2, "Stalling++" = 0.2, "Stalling+\nsubtracted" = 0.3, "Genes with\nno stalling+" = 0.2)) +
  scale_linetype_manual(values = c("Original" = "solid", "Stalling+" = "dashed", "Stalling++" = "dotted", "Stalling+\nsubtracted" = "dotdash", "Genes with\nno stalling+" = "longdash")) +
  scale_color_manual(values = c("All" = "black", "GC3" = "royalblue", "AT3" = "darkorange")) +
  scale_fill_manual(values = c("All" = "darkcyan", "GC3" = "royalblue", "AT3" = "darkorange")) +
  
  # ファセット設定
  facet_wrap(~species, ncol = 1, scales = "free_y") +
  coord_cartesian(xlim = c(-40, 40)) +
  
  # --- デザインのカスタマイズ ---
  theme_bw() +
  theme(
    # 種名をプロット領域内に寄せる設定
    strip.background = element_blank(),
    strip.text = element_text(
      face = "italic",      # 学名なのでイタリックに
      size = 16, 
      hjust = 0,         # 左端に寄せる
      vjust = 0.5           # 少し上に浮かせる
    ),
    axis.text.x = element_text(size = 14, face = "plain", color = "black"),
    axis.text.y = element_text(size = 14, face = "plain", color = "black"),
      
    # 軸のタイトル（Position, -log10Q...）を大きくする
    axis.title.x = element_text(size = 16, face = "plain"),
    axis.title.y = element_text(size = 16, face = "plain"),

    legend.text = element_text(size = 14, lineheight = 0.8),
    legend.title = element_text(size = 16, face = "plain"),

    legend.spacing.y = unit(0.5, "cm"),           # 項目間の垂直方向の隙間を広げる
    legend.key.height = unit(1.2, "cm"),          # 改行に合わせて各項目の高さを確保

    # 軸やパネルの微調整
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.3, "lines"), # パネル間の隙間を少し詰める
    legend.position = "right"
  ) +
  labs(x = "Codon position (0 = P site)",
    y = expression(-log[10](Q) %*% r),
    # 凡例のタイトルを指定
    color = "Codon group",      # codon_gp に対応
    fill = "Codon group",       # fill に対応
    linetype = "Data category", # data_gp に対応
    alpha = "Data category"     # alpha に対応)
  )
print(final_facet_plot)
ggsave("universal_weight_pival_comparison_facet_GC3AT3.png", plot = final_facet_plot, width = 6, height = 12, dpi = 450)

final_facet_plot2 <- ggplot(full_plot_df2, aes(x = aligned_pos, y = mean, 
                                             color = codon_gp, 
                                             fill = codon_gp,
                                             linetype = data_gp_eds, 
                                             alpha = data_gp_eds,
                                             group = interaction(data_gp_eds, codon_gp))) +
  # 【影の追加】 リボソーム被覆部位などの範囲指定（例：-15から0まで）
  # geom_lineより前に書くことで最背面に配置されます
  annotate("rect", xmin = -5, xmax = 5, ymin = -Inf, ymax = Inf, 
           alpha = 0.3, fill = "gold") + 
  
  # 基準線
  geom_hline(yintercept = 0, color = "black", alpha = 0.3) +
  
  # リボンと線
  geom_ribbon(aes(ymin = ymin, ymax = ymax, alpha = data_gp_eds), color = NA) +
  geom_line(linewidth = 0.7, alpha=0.5) +
  
  # スケール設定（既存の設定を維持）
  scale_alpha_manual(values = c("Original" = 0.4, "Stalling+" = 0.2, "Stalling++" = 0.2, "Stalling+\nsubtracted" = 0.3, "Genes with\nno stalling+" = 0.2)) +
  scale_linetype_manual(values = c("Original" = "solid", "Stalling+" = "dashed", "Stalling++" = "dotted", "Stalling+\nsubtracted" = "dotdash", "Genes with\nno stalling+" = "longdash")) +
  scale_color_manual(values = c("All" = "black", "GT3" = "royalblue", "AC3" = "darkorange")) +
  scale_fill_manual(values = c("All" = "darkcyan", "GT3" = "royalblue", "AC3" = "darkorange")) +
  
  # ファセット設定
  facet_wrap(~species, ncol = 1, scales = "free_y") +
  coord_cartesian(xlim = c(-40, 40)) +
  
  # --- デザインのカスタマイズ ---
  theme_bw() +
  theme(
    # 種名をプロット領域内に寄せる設定
    strip.background = element_blank(),
    strip.text = element_text(
      face = "italic",      # 学名なのでイタリックに
      size = 16, 
      hjust = 0,         # 左端に寄せる
      vjust = 0.5           # 少し上に浮かせる
    ),
    axis.text.x = element_text(size = 14, face = "plain", color = "black"),
    axis.text.y = element_text(size = 14, face = "plain", color = "black"),
      
    # 軸のタイトル（Position, -log10Q...）を大きくする
    axis.title.x = element_text(size = 16, face = "plain"),
    axis.title.y = element_text(size = 16, face = "plain"),

    legend.text = element_text(size = 14, lineheight = 0.8),
    legend.title = element_text(size = 16, face = "plain"),

    legend.spacing.y = unit(0.5, "cm"),           # 項目間の垂直方向の隙間を広げる
    legend.key.height = unit(1.2, "cm"),          # 改行に合わせて各項目の高さを確保

    # 軸やパネルの微調整
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.3, "lines"), # パネル間の隙間を少し詰める
    legend.position = "right"
  ) +
  labs(x = "Codon position (0 = P site)",
    y = expression(-log[10](Q) %*% r),
    # 凡例のタイトルを指定
    color = "Codon group",      # codon_gp に対応
    fill = "Codon group",       # fill に対応
    linetype = "Data category", # data_gp に対応
    alpha = "Data category"     # alpha に対応)
  )
print(final_facet_plot2)
ggsave("universal_weight_pival_comparison_facet_GT3AC3.png", plot = final_facet_plot2, width = 6, height = 12, dpi = 450)

final_facet_plot3 <- ggplot(full_plot_df3, aes(x = aligned_pos, y = mean, 
                                             color = codon_gp, 
                                             fill = codon_gp,
                                             linetype = data_gp_eds, 
                                             alpha = data_gp_eds,
                                             group = interaction(data_gp_eds, codon_gp))) +
  # 【影の追加】 リボソーム被覆部位などの範囲指定（例：-15から0まで）
  # geom_lineより前に書くことで最背面に配置されます
  annotate("rect", xmin = -5, xmax = 5, ymin = -Inf, ymax = Inf, 
           alpha = 0.3, fill = "gold") + 
  
  # 基準線
  geom_hline(yintercept = 0, color = "black", alpha = 0.3) +
  
  # リボンと線
  geom_ribbon(aes(ymin = ymin, ymax = ymax, alpha = data_gp_eds), color = NA) +
  geom_line(linewidth = 0.7, alpha=0.5) +
  
  # スケール設定（既存の設定を維持）
  scale_alpha_manual(values = c("Original" = 0.4, "Stalling+" = 0.2, "Stalling++" = 0.2, "Stalling+\nsubtracted" = 0.3, "Genes with\nno stalling+" = 0.2)) +
  scale_linetype_manual(values = c("Original" = "solid", "Stalling+" = "dashed", "Stalling++" = "dotted", "Stalling+\nsubtracted" = "dotdash", "Genes with\nno stalling+" = "longdash")) +
  scale_color_manual(values = c("All" = "black", "GA3" = "royalblue", "TC3" = "darkorange")) +
  scale_fill_manual(values = c("All" = "darkcyan", "GA3" = "royalblue", "TC3" = "darkorange")) +
  
  # ファセット設定
  facet_wrap(~species, ncol = 1, scales = "free_y") +
  coord_cartesian(xlim = c(-40, 40)) +
  
  # --- デザインのカスタマイズ ---
  theme_bw() +
  theme(
    # 種名をプロット領域内に寄せる設定
    strip.background = element_blank(),
    strip.text = element_text(
      face = "italic",      # 学名なのでイタリックに
      size = 16, 
      hjust = 0,         # 左端に寄せる
      vjust = 0.5           # 少し上に浮かせる
    ),
    axis.text.x = element_text(size = 14, face = "plain", color = "black"),
    axis.text.y = element_text(size = 14, face = "plain", color = "black"),
      
    # 軸のタイトル（Position, -log10Q...）を大きくする
    axis.title.x = element_text(size = 16, face = "plain"),
    axis.title.y = element_text(size = 16, face = "plain"),

    legend.text = element_text(size = 14, lineheight = 0.8),
    legend.title = element_text(size = 16, face = "plain"),

    legend.spacing.y = unit(0.5, "cm"),           # 項目間の垂直方向の隙間を広げる
    legend.key.height = unit(1.2, "cm"),          # 改行に合わせて各項目の高さを確保

    # 軸やパネルの微調整
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.3, "lines"), # パネル間の隙間を少し詰める
    legend.position = "right"
  ) +
  labs(x = "Codon position (0 = P site)",
    y = expression(-log[10](Q) %*% r),
    # 凡例のタイトルを指定
    color = "Codon group",      # codon_gp に対応
    fill = "Codon group",       # fill に対応
    linetype = "Data category", # data_gp に対応
    alpha = "Data category"     # alpha に対応)
  )
print(final_facet_plot3)
ggsave("universal_weight_pival_comparison_facet_GA3TC3.png", plot = final_facet_plot3, width = 6, height = 12, dpi = 450)




###################################

#For codons of different GC contents
get_gc_stats <- function(u_clean, group_list4, d_name) {
  # 4つのグループの平均密度を一括取得
  res <- lapply(names(group_list4), function(gn) {
    m <- colMeans(u_clean[group_list4[[gn]], , drop = F], na.rm = TRUE)
    # 簡易的にSD/SEを計算（必要に応じてユニット数nで割る処理を外側で追加）
    data.frame(pos = 1:length(m), mean = m, gc_gp = gn, data_gp = d_name)
  })
  return(do.call(rbind, res))
}

results_list_gc_grad <- list()

for (sp in species_key) {
  
  results_list_gc_grad <- list()
  
  coeffs <- coeffs_all[grepl(paste0(sp,"_thr0"),rownames(coeffs_all)),codon_order_wo_stop,drop=F]
  coeffs<-t(coeffs)
  if (sp == "Arabidopsis_thaliana_pt"){
      coeffs<-read.table("Arabidopsis_thaliana_all_cm_gene_tpm.mean.linear.coeffs.tsv",row.names=1,header=T,sep="\t")
      coeffs<-coeffs[codon_order_wo_stop,,drop=F]
  }
  # 1. データの読み込みとリスト化
  data_list <- list(
    Original = counts[grepl(sp, counts[, ncol(counts)]) & !grepl("_sub0", counts[, ncol(counts)]) & !grepl("flat", counts[, ncol(counts)]), , drop = F],
    Sub0     = counts[grepl(sp, counts[, ncol(counts)]) & grepl("_sub0", counts[, ncol(counts)]) & !grepl("_sub000", counts[, ncol(counts)]), , drop = F],
    Sub000   = counts[grepl(sp, counts[, ncol(counts)]) & grepl("_sub000", counts[, ncol(counts)]), , drop = F],
    Flat_regions   = counts[grepl(sp, counts[, ncol(counts)]) & grepl("_flat_regions", counts[, ncol(counts)]), , drop = F],
    Flat_genes   = counts[grepl(sp, counts[, ncol(counts)]) & grepl("_flat_genes", counts[, ncol(counts)]), , drop = F]
  )
  if (sp == "Arabidopsis_thaliana"){
      data_list <- list(
        Original = counts[grepl(sp, counts[, ncol(counts)]) & !grepl("_sub0", counts[, ncol(counts)]) & !grepl("flat", counts[, ncol(counts)]) & !grepl("_pt_", counts[, ncol(counts)]), , drop = F],
        Sub0     = counts[grepl(sp, counts[, ncol(counts)]) & grepl("_sub0", counts[, ncol(counts)]) & !grepl("_sub000", counts[, ncol(counts)]) & !grepl("_pt_", counts[, ncol(counts)]), , drop = F],
        Sub000   = counts[grepl(sp, counts[, ncol(counts)]) & grepl("_sub000", counts[, ncol(counts)]) & !grepl("_pt_", counts[, ncol(counts)]), , drop = F],
        Flat_regions   = counts[grepl(sp, counts[, ncol(counts)]) & grepl("_flat_regions", counts[, ncol(counts)]) & !grepl("_pt_", counts[, ncol(counts)]), , drop = F],
        Flat_genes   = counts[grepl(sp, counts[, ncol(counts)]) & grepl("_flat_genes", counts[, ncol(counts)]) & !grepl("_pt_", counts[, ncol(counts)]), , drop = F]
      )
  } else if (sp == "Arabidopsis_thaliana_pt"){
      sp0<-gsub("_pt","",sp)
      data_list <- list(
        Original = counts[grepl(sp0, counts[, ncol(counts)]) & !grepl("_sub0", counts[, ncol(counts)]) & !grepl("flat", counts[, ncol(counts)]) & grepl("_pt", counts[, ncol(counts)]), , drop = F],
        Sub0     = counts[grepl(sp0, counts[, ncol(counts)]) & grepl("_sub0", counts[, ncol(counts)]) & !grepl("_sub000", counts[, ncol(counts)]) & grepl("_pt_", counts[, ncol(counts)]), , drop = F],
        Sub000   = counts[grepl(sp0, counts[, ncol(counts)]) & grepl("_sub000", counts[, ncol(counts)]) & grepl("_pt_", counts[, ncol(counts)]), , drop = F],
        Flat_regions   = counts[grepl(sp0, counts[, ncol(counts)]) & grepl("_flat_regions", counts[, ncol(counts)]) & grepl("_pt_", counts[, ncol(counts)]), , drop = F],
        Flat_genes   = counts[grepl(sp0, counts[, ncol(counts)]) & grepl("_flat_genes", counts[, ncol(counts)]) & grepl("_pt_", counts[, ncol(counts)]), , drop = F]
      )
  }  

  for (d_name in names(data_list)) {
    curr_data <- data_list[[d_name]]
    if (nrow(curr_data) == 0) next  # データが空の場合はスキップ

    unit_size <- 64 
    num_units <- nrow(curr_data) / unit_size
    num_pos   <- ncol(curr_data) - 2

    # --- 行列の初期化はココ (d_nameごとに作成) ---
    mtrx_gc0 <- matrix(NA, nrow = num_units, ncol = num_pos)
    mtrx_gc1 <- matrix(NA, nrow = num_units, ncol = num_pos)
    mtrx_gc2 <- matrix(NA, nrow = num_units, ncol = num_pos)
    mtrx_gc3 <- matrix(NA, nrow = num_units, ncol = num_pos)

    # ユニットごとの集計
    for (j in 1:num_units) {
      start_idx <- (j - 1) * unit_size + 1
      end_idx   <- j * unit_size
      u_raw     <- curr_data[start_idx:end_idx, ]
      u_clean   <- u_raw[match(codon_order_wo_stop, u_raw[,1]), 2:(num_pos + 1), drop = F]
      
      mtrx_gc0[j, ] <- colMeans(u_clean[group_list4$AT_only, , drop = F], na.rm = TRUE)
      mtrx_gc1[j, ] <- colMeans(u_clean[group_list4$GC_one,  , drop = F], na.rm = TRUE)
      mtrx_gc2[j, ] <- colMeans(u_clean[group_list4$GC_two,  , drop = F], na.rm = TRUE)
      mtrx_gc3[j, ] <- colMeans(u_clean[group_list4$GC_only, , drop = F], na.rm = TRUE)
    }
    
    # 統計計算関数
    calc_stats <- function(mtrx, group_name, d_name) {
      m <- colMeans(mtrx, na.rm = TRUE)
      n_vec <- rowSums(!is.na(mtrx))
      s <- apply(mtrx, 2, sd, na.rm = TRUE) / sqrt(pmax(n_vec, 1)) # 0除算防止
      return(data.frame(pos = 1:length(m), mean = m, ymin = m - s, ymax = m + s, 
                        gc_gp = group_name, data_gp = d_name))
    }

    # 各グループをリストに格納
    results_list_gc_grad[[paste0(d_name, "_GC0")]] <- calc_stats(mtrx_gc0, "GC0", d_name)
    results_list_gc_grad[[paste0(d_name, "_GC1")]] <- calc_stats(mtrx_gc1, "GC1", d_name)
    results_list_gc_grad[[paste0(d_name, "_GC2")]] <- calc_stats(mtrx_gc2, "GC2", d_name)
    results_list_gc_grad[[paste0(d_name, "_GC3")]] <- calc_stats(mtrx_gc3, "GC3", d_name)
  }

# --- プロットは d_name ループの外、sp ループの内側で行う ---
  if (length(results_list_gc_grad) > 0) {
    
    # 【重要】リストをデータフレームに変換する工程を追加
    plot_gc_grad <- do.call(rbind, results_list_gc_grad)

    # 1. 表示名と順番の定義
    display_names <- c(
      "Original"     = "Original",
      "Sub0"         = "Stalling+",
      "Sub000"       = "Stalling++",
      "Flat_regions" = "Stalling+\nsubtracted",
      "Flat_genes"   = "Genes with\nno stalling+"
    )

    # 2. データの因子型（factor）化（これをしないと順番がバラバラになります）
    plot_gc_grad$data_gp <- factor(
      plot_gc_grad$data_gp, 
      levels = names(display_names)
    )
    plot_gc_grad$gc_gp <- factor(plot_gc_grad$gc_gp, levels = c("GC0", "GC1", "GC2", "GC3"))

    # 3. プロット
    g_plot <- ggplot(plot_gc_grad, aes(x = pos, y = mean, color = gc_gp, fill = gc_gp, 
                                     group = interaction(data_gp, gc_gp))) +
      geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.1, color = NA) +
      geom_line(aes(linetype = data_gp), linewidth = 0.8) +
      
      scale_color_viridis_d(option = "plasma", end = 0.8) +
      scale_fill_viridis_d(option = "plasma", end = 0.8) +
      
      scale_linetype_manual(
        values = c("Original"="solid", "Sub0"="dashed", "Sub000"="dotted", 
                   "Flat_regions"="dotdash", "Flat_genes"="longdash"),
        labels = display_names,
        breaks = names(display_names)
      ) +
      
      theme_bw() +
      theme(
        axis.text = element_text(size = 14, face = "plain"),
        axis.title = element_text(size = 16, face = "plain"),
        legend.title = element_text(size = 14, face = "plain"),
        legend.text = element_text(size = 14, lineheight = 0.9),
        legend.key.height = unit(2.5, "lines"),
        legend.spacing.y = unit(0.4, "cm")
      ) +
      guides(
        linetype = guide_legend(byrow = TRUE, order = 1),
        color = guide_legend(order = 2),
        fill = guide_legend(order = 2)
      ) +
      labs(
        title = paste("Species:", gsub("_"," ",sp)),
        x = "Position", 
        y = "Normalized Ribosome Density", # 文脈に合わせて修正
        linetype = "Data category", 
        color = "GC count in codon",        
        fill = "GC count in codon"
      )
    
    # 保存
    ggsave(paste0(sp, "_GC_gradient_plot.png"), plot = g_plot, width = 10, height = 7, dpi = 450)
    
    # メモリ節約のため、1種終わるごとに表示してクリア
    print(g_plot)
  }
}

##################################
#For data browsing
sum(is.na(emtrx))
plot(unlist(coeffs),org_u[,56]) #This gave the best separation in rice data
ccounts_unit[ccounts_unit[,56] > 0.03,1] #GT3 was relevant to spliting the group

scplot <- function(datamtrx, coeff_data, pattern = "..[GT]$", label = "GT3", colnum) {
  # 1. プロット用データの作成
  plot_df <- data.frame(
    Codon = rownames(datamtrx),
    Coeff = as.numeric(unlist(coeff_data)), # x軸：回帰係数
    Value = as.numeric(datamtrx[, colnum])  # y軸：指定された列の値
  )
  
  # 2. 判定ラベルの作成（引数 pattern で判定）
  # label に指定した名前と "Other" でグループ分けします
  plot_df$Group <- ifelse(grepl(pattern, plot_df$Codon), label, "Other")
  
  # 3. ggplot2による散布図の描画
  dp <- ggplot(plot_df, aes(x = Coeff, y = Value, color = Group)) +
    # 背景に回帰線を引くと、グループごとの「切片の差」がより明確になります
    geom_smooth(method = "lm", se = FALSE, size = 0.5, alpha = 0.5) +
    geom_point(size = 3, alpha = 0.8) +
    geom_text(aes(label = Codon), vjust = -1, size = 3, check_overlap = TRUE) +
    # 色の設定（label の名前に合わせて動的に設定）
    scale_color_manual(values = setNames(c("royalblue", "darkorange"), c(label, "Other"))) +
    theme_classic() +
    labs(
      title = paste("Scatter Plot: Column", colnum, "(Focus:", label, ")"),
      x = "Regression Coefficients",
      y = paste("Counts (Col", colnum, ")"),
      color = "Category"
    )
  
  print(dp)
}

scplot(org_u, coeffs, pattern = "..[GT]$", label = "GC3", 56)

#To make scatter plots for interesting sites
counts<-read.table("all_ribo_pileup.txt", row.names=NULL, header=T, sep="\t")
coeffs_all<-read.table("all_species_sum.coeff.linear.txt", row.names=1, header=T, sep="\t")

scplot2 <- function(datamtrx, coeff_data, pattern = "..[GC]$", label1 = "GC3", label2 = "AT3", colnum=1, color1="royalblue", color2="darkorange", wd=3, ht=3, prefix="") {
  # 1. データの作成
  plot_df <- data.frame(
    Codon = rownames(datamtrx),
    Coeff = as.numeric(as.matrix(coeff_data)),
    Value = as.numeric(as.matrix(datamtrx[, colnum]))
  )
  
  # 2. グループ分け（label2を"Other"の代わりに直接使用）
  plot_df$Group <- ifelse(grepl(pattern, plot_df$Codon), label1, label2)
  
  # 統計量
  res <- cor.test(plot_df$Coeff, plot_df$Value, method = "pearson")
  r_val <- round(res$estimate, 3)
  p_val_actual <- res$p.value
  
  # P値の表記：0.001以上なら小数点3位、それ未満なら指数表記で「そのまま」出す
  p_label <- if(p_val_actual < 0.001) {
               formatC(p_val_actual, format = "e", digits = 2) # 例: 1.50e-11
             } else {
               round(p_val_actual, 3)
             }
  
  # 論文で使いやすいように「r = ..., p = ...」の形式に
  stats_text <- paste0("r = ", r_val, "\n", "p = ", p_label)


  # 3. 描画
  dp <- ggplot(plot_df, aes(x = Coeff, y = Value)) +
    # 全体回帰線
    geom_smooth(method = "lm", se = FALSE, color = "orchid", linewidth = 0.5, linetype = "dashed") +
    # 点プロット
    geom_point(aes(color = Group), size = 2.5, alpha = 0.6) +
    # 【修正】コドン文字の色を点と合わせる（重なり回避付き）
    #geom_text_repel(aes(label = Codon, color = Group), 
    #                size = 2.5,             # 少し小さくすると重なりにくい
    #                alpha = 0.4,
    #                max.overlaps = Inf,     # 重なり制限を解除（必須）
    #                box.padding = 0.5,      # 文字の周りの余白（広げると重なりにくい）
    #                point.padding = 0.2,    # 点から少し離す
    #                segment.size = 0.2,     # 引き出し線の太さ
    #                segment.alpha = 0.4,    # 引き出し線の透明度
    #                force = 2,              # 反発力を強める
    #                show.legend = FALSE) + 
    # 色の設定
    scale_color_manual(values = setNames(c(color1, color2), c(label1, label2))) +
    # 統計指標
    annotate("text", x = min(plot_df$Coeff), y = Inf, label = stats_text, 
             hjust = 0, vjust = 1.2, size = 4, fontface = "italic") +
    # デザイン
    theme_classic() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      axis.line = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 10),
      legend.position = "none"
    ) +
    labs(
      title = paste("Site", (colnum-51)),
      x = "Codon weight",
      y = "Sum of normalized counts",
      color = "Codon Type"
    )
  
  print(dp)
  # dpi=450に修正
  ggsave(paste0(prefix, "scplot2_out_", label1, ".png"), plot = dp, width = wd, height = ht, dpi = 450)  
  return(res)
}

# 補助関数：P値の整形用（任意）
format.plain <- function(x) {
  # epsより小さければ "< 1e-06" 等の形式、それ以上なら小数点3位で整形
  format.pval(x, eps = 0.00001, digits = 3, scientific = FALSE)
}

#For Arabidopsis
sp<-"Arabidopsis_thaliana"
coeffs <- coeffs_all[grepl(paste0(sp,"_thr0"),rownames(coeffs_all)),codon_order_wo_stop,drop=F]
coeffs<-t(coeffs)

data_list <- list(
    Original = counts[grepl(sp, counts[, ncol(counts)]) & !grepl("_sub0", counts[, ncol(counts)]) & !grepl("flat", counts[, ncol(counts)]), , drop = F],
    Sub0     = counts[grepl(sp, counts[, ncol(counts)]) & grepl("_sub0", counts[, ncol(counts)]) & !grepl("_sub000", counts[, ncol(counts)]), , drop = F],
    Sub000   = counts[grepl(sp, counts[, ncol(counts)]) & grepl("_sub000", counts[, ncol(counts)]), , drop = F],
    Flat_regions   = counts[grepl(sp, counts[, ncol(counts)]) & grepl("_flat_regions", counts[, ncol(counts)]), , drop = F],
    Flat_genes   = counts[grepl(sp, counts[, ncol(counts)]) & grepl("_flat_genes", counts[, ncol(counts)]), , drop = F]
)
if (sp == "Arabidopsis_thaliana"){
      data_list <- list(
        Original = counts[grepl(sp, counts[, ncol(counts)]) & !grepl("_sub0", counts[, ncol(counts)]) & !grepl("flat", counts[, ncol(counts)]) & !grepl("_pt", counts[, ncol(counts)]), , drop = F],
        Sub0     = counts[grepl(sp, counts[, ncol(counts)]) & grepl("_sub0", counts[, ncol(counts)]) & !grepl("_sub000", counts[, ncol(counts)]) & !grepl("_pt", counts[, ncol(counts)]), , drop = F],
        Sub000   = counts[grepl(sp, counts[, ncol(counts)]) & grepl("_sub000", counts[, ncol(counts)]) & !grepl("_pt", counts[, ncol(counts)]), , drop = F],
        Flat_regions   = counts[grepl(sp, counts[, ncol(counts)]) & grepl("_flat_regions", counts[, ncol(counts)]) & !grepl("_pt", counts[, ncol(counts)]), , drop = F],
        Flat_genes   = counts[grepl(sp, counts[, ncol(counts)]) & grepl("_flat_genes", counts[, ncol(counts)]) & !grepl("_pt", counts[, ncol(counts)]), , drop = F]
      )
  } else if (sp == "Arabidopsis_thaliana_pt"){
      sp0<-gsub("_pt","",sp)
      data_list <- list(
        Original = counts[grepl(sp0, counts[, ncol(counts)]) & !grepl("_sub0", counts[, ncol(counts)]) & !grepl("flat", counts[, ncol(counts)]) & grepl("_pt", counts[, ncol(counts)]), , drop = F],
        Sub0     = counts[grepl(sp0, counts[, ncol(counts)]) & grepl("_sub0", counts[, ncol(counts)]) & !grepl("_sub000", counts[, ncol(counts)]) & grepl("_pt", counts[, ncol(counts)]), , drop = F],
        Sub000   = counts[grepl(sp0, counts[, ncol(counts)]) & grepl("_sub000", counts[, ncol(counts)]) & grepl("_pt", counts[, ncol(counts)]), , drop = F],
        Flat_regions   = counts[grepl(sp0, counts[, ncol(counts)]) & grepl("_flat_regions", counts[, ncol(counts)]) & grepl("_pt", counts[, ncol(counts)]), , drop = F],
        Flat_genes   = counts[grepl(sp0, counts[, ncol(counts)]) & grepl("_flat_genes", counts[, ncol(counts)]) & grepl("_pt", counts[, ncol(counts)]), , drop = F]
      )
}  

pis_at<-read.table("Arabidopsis_thaliana_Original_GC3AT3_ribbon_plot_pival.txt", row.names=1, header=T, sep="\t")
maxindices <- which(pis_at == max(pis_at, na.rm = TRUE), arr.ind = TRUE)
max_row_idx <- maxindices[1, "row"]
max_col_idx <- maxindices[1, "col"]
start_idx <- (max_row_idx - 1) * 64 + 1
end_idx   <- max_row_idx * 64
countsoi <- data_list[["Original"]][start_idx:end_idx, ]
rownames(countsoi) <-countsoi[,1]
countsoi<-countsoi[codon_order_wo_stop,2:102]

scplot2(countsoi, coeffs, pattern = "..[GC]$", label1 = "GC3", label2="AT3", colnum=max_col_idx, ht=2.5, prefix="Arabidopsis_thaliana_pimax")
scplot2(countsoi, coeffs, pattern = "..[GT]$", label1 = "GT3", label2="AC3", color1="darkcyan", color2="darkmagenta",colnum=51, ht=2.5, prefix="Arabidopsis_thaliana_Psite")

#For human data
sp<-"Homo_sapiens"
coeffs <- coeffs_all[grepl(paste0(sp,"_thr0"),rownames(coeffs_all)),codon_order_wo_stop,drop=F]
coeffs<-t(coeffs)

data_list <- list(
    Original = counts[grepl(sp, counts[, ncol(counts)]) & !grepl("_sub0", counts[, ncol(counts)]) & !grepl("flat", counts[, ncol(counts)]), , drop = F],
    Sub0     = counts[grepl(sp, counts[, ncol(counts)]) & grepl("_sub0", counts[, ncol(counts)]) & !grepl("_sub000", counts[, ncol(counts)]), , drop = F],
    Sub000   = counts[grepl(sp, counts[, ncol(counts)]) & grepl("_sub000", counts[, ncol(counts)]), , drop = F],
    Flat_regions   = counts[grepl(sp, counts[, ncol(counts)]) & grepl("_flat_regions", counts[, ncol(counts)]), , drop = F],
    Flat_genes   = counts[grepl(sp, counts[, ncol(counts)]) & grepl("_flat_genes", counts[, ncol(counts)]), , drop = F]
)
if (sp == "Arabidopsis_thaliana"){
      data_list <- list(
        Original = counts[grepl(sp, counts[, ncol(counts)]) & !grepl("_sub0", counts[, ncol(counts)]) & !grepl("flat", counts[, ncol(counts)]) & !grepl("_pt", counts[, ncol(counts)]), , drop = F],
        Sub0     = counts[grepl(sp, counts[, ncol(counts)]) & grepl("_sub0", counts[, ncol(counts)]) & !grepl("_sub000", counts[, ncol(counts)]) & !grepl("_pt", counts[, ncol(counts)]), , drop = F],
        Sub000   = counts[grepl(sp, counts[, ncol(counts)]) & grepl("_sub000", counts[, ncol(counts)]) & !grepl("_pt", counts[, ncol(counts)]), , drop = F],
        Flat_regions   = counts[grepl(sp, counts[, ncol(counts)]) & grepl("_flat_regions", counts[, ncol(counts)]) & !grepl("_pt", counts[, ncol(counts)]), , drop = F],
        Flat_genes   = counts[grepl(sp, counts[, ncol(counts)]) & grepl("_flat_genes", counts[, ncol(counts)]) & !grepl("_pt", counts[, ncol(counts)]), , drop = F]
      )
  } else if (sp == "Arabidopsis_thaliana_pt"){
      sp0<-gsub("_pt","",sp)
      data_list <- list(
        Original = counts[grepl(sp0, counts[, ncol(counts)]) & !grepl("_sub0", counts[, ncol(counts)]) & !grepl("flat", counts[, ncol(counts)]) & grepl("_pt", counts[, ncol(counts)]), , drop = F],
        Sub0     = counts[grepl(sp0, counts[, ncol(counts)]) & grepl("_sub0", counts[, ncol(counts)]) & !grepl("_sub000", counts[, ncol(counts)]) & grepl("_pt", counts[, ncol(counts)]), , drop = F],
        Sub000   = counts[grepl(sp0, counts[, ncol(counts)]) & grepl("_sub000", counts[, ncol(counts)]) & grepl("_pt", counts[, ncol(counts)]), , drop = F],
        Flat_regions   = counts[grepl(sp0, counts[, ncol(counts)]) & grepl("_flat_regions", counts[, ncol(counts)]) & grepl("_pt", counts[, ncol(counts)]), , drop = F],
        Flat_genes   = counts[grepl(sp0, counts[, ncol(counts)]) & grepl("_flat_genes", counts[, ncol(counts)]) & grepl("_pt", counts[, ncol(counts)]), , drop = F]
      )
}  

pis_hs<-read.table("Homo_sapiens_Original_GC3AT3_ribbon_plot_pival.txt", row.names=1, header=T, sep="\t")
maxindices <- which(pis_hs == max(pis_hs, na.rm = TRUE), arr.ind = TRUE)
max_row_idx <- maxindices[1, "row"]
max_col_idx <- maxindices[1, "col"]
start_idx <- (max_row_idx - 1) * 64 + 1
end_idx   <- max_row_idx * 64
countsoi <- data_list[["Original"]][start_idx:end_idx, ]
rownames(countsoi) <-countsoi[,1]
countsoi<-countsoi[codon_order_wo_stop,2:102]

scplot2(countsoi, coeffs, pattern = "..[GC]$", label1 = "GC3", label2="AT3", colnum=max_col_idx, ht=2.5, prefix="Homo_sapiens_pimax")
scplot2(countsoi, coeffs, pattern = "..[GT]$", label1 = "GT3", label2="AC3", color1="darkcyan", color2="darkmagenta",colnum=56, ht=2.5, prefix="Homo_sapiens_entrysite")
