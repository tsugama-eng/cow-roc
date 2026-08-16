library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(dplyr)

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

#########################################
#For the PA-predicted TPM correlation with the seven species

species_labels <- c(
  "Arabidopsis_thaliana" = "Arabidopsis\nthaliana",
  "Oryza_sativa" = "Oryza\nsativa",
  "Homo_sapiens" = "Homo\nsapiens",
  "Mus_musculus" = "Mus\nmusculus",
  "Drosophila_melanogaster" = "Drosophila\nmelanogaster",
  "Saccharomyces_cerevisiae" = "Saccharomyces\ncerevisiae",
  "Escherichia_coli" = "Escherichia\ncoli"
)

plot_data_list <- list()
all_pair_data_list <- list()

for (sp in names(species_labels)) {
  exp_file  <- paste0(sp, "_protein_abundance.txt")
  ind_pred_file <- paste0(sp, "_GC3_CAIhigh_predTPM.txt")
  
  if (!file.exists(exp_file) | !file.exists(ind_pred_file)) next

  # --- 読み込み方法の変更（row.names=NULLにして1列目を安全に処理） ---
  exp_raw <- read.table(exp_file, header=T, sep="\t", check.names=F, stringsAsFactors=F)
  ind_raw <- read.table(ind_pred_file, header=T, sep="\t", check.names=F, stringsAsFactors=F)

  # IDクリーニング関数：gene: や ENSB: を消し、FBgn... だけを残す
  clean_id <- function(df) {
    ids <- as.character(df[,1])
    ids <- gsub("^gene:|^ENSB:|^AT[0-9]G|^Os[0-9]{2}g", "", ids)
    ids <- gsub("^.*:", "", ids) # コロンがあればその後ろだけ取る
    ids <- trimws(ids)
    # 重複があれば最初の行だけ残す（エラー防止）
    df <- df[!duplicated(ids), ]
    rownames(df) <- ids[!duplicated(ids)]
    return(df[,-1, drop=FALSE]) # 1列目を除去して返す
  }

  exp <- clean_id(exp_raw)
  ind_pred_table <- clean_id(ind_raw)

  # ターゲット列の特定（grepで柔軟に）
  if (sp %in% c("Homo_sapiens", "Mus_musculus", "Arabidopsis_thaliana", "Oryza_sativa")) {
    mapping <- list(
      "GC3" = grep("GC3", colnames(ind_pred_table), value=T),
      "stAI" = grep("stAI", colnames(ind_pred_table), value=T),
      "CAI" = grep("CAI", colnames(ind_pred_table), value=T),
      "mRNA_pred" = grep("pred.logTPM|mRNA_pred", colnames(ind_pred_table), value=T)
    )
  } else {
    mapping <- list(
      "GC3" = grep("GC3", colnames(ind_pred_table), value=T),
      "CAI" = grep("CAI", colnames(ind_pred_table), value=T),
      "mRNA_pred" = grep("pred.logTPM|mRNA_pred", colnames(ind_pred_table), value=T)
    )
  }

  all_results <- data.frame()
  
  for (i in 1:ncol(exp)) {
    for (label in names(mapping)) {
      idx_col <- mapping[[label]]
      if (length(idx_col) == 0) next
      
      # 共通遺伝子の抽出
      # exp側が > 0 かつ NAでない
      valid_exp <- rownames(exp)[!is.na(exp[,i]) & exp[,i] > 0]
      # ind側が NAでない
      valid_ind <- rownames(ind_pred_table)[!is.na(ind_pred_table[, idx_col[1]])]
      
      genes_used <- intersect(valid_exp, valid_ind)
      
      # デバッグ表示（これで各指標の共通数が見えます）
      if(i == 1) message(paste(sp, label, "Overlap:", length(genes_used)))

      if(length(genes_used) < 20) next
      
      v_exp <- as.numeric(exp[genes_used, i])
      v_ind <- as.numeric(ind_pred_table[genes_used, idx_col[1]])
      
      ok <- !is.na(v_exp) & !is.na(v_ind)
      if(sum(ok) < 20) next
      
      corres <- try(cor.test(log10(v_exp[ok]), v_ind[ok]), silent=TRUE)
      if(class(corres) == "try-error") next
      
      all_results <- rbind(all_results, data.frame(
        Sample = colnames(exp)[i],
        IndexGroup = label,
        Pval = corres$p.value,
        Rho = corres$estimate,
        stringsAsFactors = FALSE
      ))
    }
  }

  if(nrow(all_results) > 0) {
    all_results$Qval <- p.adjust(all_results$Pval, method="fdr")
    all_results$Qval[all_results$Qval == 0] <- 1e-300
    all_results$Score <- -log10(all_results$Qval) * all_results$Rho
    all_results$Species <- sp
    plot_data_list[[sp]] <- all_results
  }
  
  # 各生物種ごとに、特定の指標（例：mRNA_pred）の全ペアデータを保存する箱
  sp_merged_data <- data.frame()

  for (i in 1:ncol(exp)) {
    # 指標は mRNA_pred など一つに絞るか、ループで回す
    label <- "mRNA_pred" 
    idx_col <- mapping[[label]]
    if (length(idx_col) == 0) next
    
    genes_used <- intersect(rownames(exp)[!is.na(exp[,i]) & exp[,i] > 0], 
                            rownames(ind_pred_table)[!is.na(ind_pred_table[, idx_col[1]])])
    
    if(length(genes_used) < 20) next
    
    # 全データを集積
    tmp_df <- data.frame(
      pred = as.numeric(ind_pred_table[genes_used, idx_col[1]]),
      obs  = log10(as.numeric(exp[genes_used, i])),
      Species = sp
    )
    sp_merged_data <- rbind(sp_merged_data, tmp_df)
  }
  
  # 400万点になる場合は、ここで軽くダウンサンプリングしておくと後の処理が楽です
  if(nrow(sp_merged_data) > 200000) {
    sp_merged_data <- sp_merged_data[sample(nrow(sp_merged_data), 200000), ]
  }
  all_pair_data_list[[sp]] <- sp_merged_data

}

# --- プロット処理 ---
final_df <- do.call(rbind, plot_data_list) %>%
  filter(!is.na(Score)) # データがない行を確実に除外
final_df$Species <- factor(final_df$Species, levels = names(species_labels), labels = species_labels)
final_df$IndexType <- factor(final_df$IndexGroup, 
                             levels = c("GC3", "stAI", "CAI", "mRNA_pred"), 
                             labels = c("GC3 content", "stAI", "CAIhigh", "Pred. logTPM"))

p <- ggplot(final_df, aes(x = IndexType, y = Score, color = IndexType)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1.2) +
  # ここを "free_y" から "free" (または "free_x") に変更
  facet_wrap(~ Species, nrow = 1, scales = "free") + 
  
  scale_color_manual(values = c("GC3 content" = "royalblue", "stAI" = "orchid", "CAIhigh" = "darkorange", "Pred. logTPM" = "darkcyan")) +
  labs(x = "", y = expression(-log[10](Q) %*% r ~ "(vs. protein abundance)")) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_text(family = "sans", face = "plain", size = 12),
    legend.position = "none",
    strip.background = element_blank(), 
    strip.text = element_text(face = "italic"),
    panel.spacing = unit(0.8, "lines") # 軸名が重ならないよう少し広げる
  )

print(p)# 保存
ggsave("correlation_jitter_facet_v2.png", p, width = 12, height = 5, dpi = 450)


# --- 最後に全生物種を結合して散布図作成 ---
final_scatter_df <- do.call(rbind, all_pair_data_list)
final_scatter_df$Species <- factor(final_scatter_df$Species, 
                                   levels = names(species_labels), 
                                   labels = species_labels)

# 散布図の作成
p_suppl <- ggplot(final_scatter_df, aes(x = pred, y = obs)) +
  # 密度表示（bin2dが最も軽く、かつ傾向が明確に見えます）
  geom_bin2d(bins = 80) + 
  scale_fill_viridis_c(option = "viridis", trans = "log10", guide = "none") +
  
  # 回帰線
  geom_smooth(method = "lm", color = "magenta", se = FALSE, linewidth = 0.8) +
  
  # 相関検定結果の追加
  # P値をそのまま（科学的表記）で出すために p.accuracy を設定
  stat_cor(method = "pearson", 
           # 0が左端、1が右端。少しだけ右に離すなら 0.02 程度
           label.x.npc = 0.001, 
           # 1が上端。1.0に近づけるほど上に寄ります
           label.y.npc = 0.995,
           size = 3.0,
           color = "black",
           label.sep = "\n",
           output.type = "text",
           # vjust = 1 でテキストの上端を y 座標に合わせる（上に詰められる）
           vjust = 1) +
  
  # 2行4列に配置（scales="free"で各生物種の軸を最適化）
  facet_wrap(~ Species, nrow = 2, ncol = 4, scales = "free") +
  
  theme_bw(base_size = 14) +
  labs(x = expression(Predicted~log[10](TPM)), 
       y = expression(Observed~log[10](protein~abundance))) +
       
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # 種名ラベル（strip）の外枠と背景を消す
    strip.background = element_blank(), 
    strip.text = element_text(face = "italic", size = 12),
    # 全体の外枠は残しつつ、パネル間の間隔を調整
    panel.spacing = unit(0.5, "lines"),
    # 不要な枠線がある場合は以下で調整
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  )

# 保存設定
ggsave("supplement_scatter_2x4.png", p_suppl, width = 10, height = 5.5, dpi = 450)

#Predicted logTPM vs observed logTPM
species_labels = c(
  # 植物（主要・穀物）
  "Oryza_sativa"              = "Oryza sativa",
  "Zea_mays"                  = "Zea mays",
  "Sorghum_bicolor"           = "Sorghum bicolor",
  "Cenchrus_americanus"       = "Cenchrus americanus",
  "Hordeum_vulgare"           = "Hordeum vulgare",
  "Triticum_aestivum"         = "Triticum aestivum",
  
  # 植物（果物・その他単子葉/双子葉）
  "Ananas_comosus"            = "Ananas comosus",
  "Musa_acuminata"            = "Musa acuminata",
  "Asparagus_officinalis"     = "Asparagus officinalis",
  "Allium_cepa"               = "Allium cepa",
  "Populus_trichocarpa"       = "Populus trichocarpa",
  "Malus_domestica"           = "Malus domestica",
  "Vitis_vinifera"            = "Vitis vinifera",
  "Citrus_clementina"         = "Citrus clementina",
  "Arabidopsis_thaliana"      = "Arabidopsis thaliana",
  "Brassica_rapa"             = "Brassica rapa",
  "Solanum_lycopersicum"      = "Solanum lycopersicum",
  "Glycine_max"               = "Glycine max",
  "Chenopodium_quinoa"        = "Chenopodium quinoa",
  
  # 植物（苔・藻類）
  "Physcomitrium_patens"      = "Physcomitrium patens",
  "Chlamydomonas_reinhardtii" = "Chlamydomonas reinhardtii",
  
  # 真菌・動物（モデル生物）
  "Homo_sapiens"              = "Homo sapiens",
  "Mus_musculus"              = "Mus musculus",
  "Drosophila_melanogaster"   = "Drosophila melanogaster",
  "Saccharomyces_cerevisiae"  = "Saccharomyces cerevisiae",
  
  # 細菌
  "Escherichia_coli"          = "Escherichia coli",
  "Bacillus_subtilis"         = "Bacillus subtilis",
  "Mycobacterium_tuberculosis"= "Mycobacterium tuberculosis",
  "Salmonella_enterica"       = "Salmonella enterica"
)

plot_data_list <- list()
all_pair_data_list <- list()

for (sp in names(species_labels)) {
  # ファイルパスのデバッグ表示
  message(paste("Checking species:", sp))
  exp_file  <- paste0(sp, "_all_gene_tpm.txt")
  ind_pred_file <- paste0(sp, "_GC3_CAIhigh_predTPM.txt")
  
  if (!file.exists(exp_file)) { message(paste("Missing exp_file:", exp_file)); next }
  if (!file.exists(ind_pred_file)) { message(paste("Missing ind_pred_file:", ind_pred_file)); next }
  
  # --- 読み込み方法の変更（row.names=NULLにして1列目を安全に処理） ---
  exp_raw <- read.table(exp_file, header=T, sep="\t", check.names=F, stringsAsFactors=F)
  ind_raw <- read.table(ind_pred_file, header=T, sep="\t", check.names=F, stringsAsFactors=F)

  # IDクリーニング関数：gene: や ENSB: を消し、FBgn... だけを残す
  clean_id <- function(df) {
    ids <- as.character(df[,1])
    ids <- gsub("^gene:|^ENSB:|^AT[0-9]G|^Os[0-9]{2}g", "", ids)
    ids <- gsub("^.*:", "", ids) # コロンがあればその後ろだけ取る
    ids <- trimws(ids)
    # 重複があれば最初の行だけ残す（エラー防止）
    df <- df[!duplicated(ids), ]
    rownames(df) <- ids[!duplicated(ids)]
    return(df[,-1, drop=FALSE]) # 1列目を除去して返す
  }

  exp <- clean_id(exp_raw)
  ind_pred_table <- clean_id(ind_raw)

  # ターゲット列の特定（grepで柔軟に）
  if (sp %in% c("Homo_sapiens", "Mus_musculus", "Arabidopsis_thaliana", "Oryza_sativa")) {
    mapping <- list(
      "GC3" = grep("GC3", colnames(ind_pred_table), value=T),
      "stAI" = grep("stAI", colnames(ind_pred_table), value=T),
      "CAI" = grep("CAI", colnames(ind_pred_table), value=T),
      "mRNA_pred" = grep("pred.logTPM|mRNA_pred", colnames(ind_pred_table), value=T)
    )
  } else {
    mapping <- list(
      "GC3" = grep("GC3", colnames(ind_pred_table), value=T),
      "CAI" = grep("CAI", colnames(ind_pred_table), value=T),
      "mRNA_pred" = grep("pred.logTPM|mRNA_pred", colnames(ind_pred_table), value=T)
    )
  }

  #all_results <- data.frame()
  res_list_sp  <- list() # 1. 相関計算の結果（点）を溜める
  pair_list_sp <- list() # 2. 散布図用の生データ（数万点）を溜める

  for (i in 1:ncol(exp)) {
    # --- 指標ごとのループ（ジッタープロット用） ---
    for (label in names(mapping)) {
      idx_col <- mapping[[label]]
      if (length(idx_col) == 0) next
      
      genes_used <- intersect(rownames(exp)[!is.na(exp[,i]) & exp[,i] > 0], 
                              rownames(ind_pred_table)[!is.na(ind_pred_table[, idx_col[1]])])
      
      if(length(genes_used) < 20) next
      
      v_exp <- as.numeric(exp[genes_used, i])
      v_ind <- as.numeric(ind_pred_table[genes_used, idx_col[1]])
      
      corres <- try(cor.test(log10(v_exp), v_ind), silent=TRUE)
      if(class(corres) == "try-error") next
      
      # リストに相関結果を保存（iとlabelを組み合わせたユニークな名前にして保存）
      res_key <- paste0(i, "_", label)
      res_list_sp[[res_key]] <- data.frame(
        Sample = colnames(exp)[i],
        IndexGroup = label,
        Pval = corres$p.value,
        Rho = corres$estimate,
        stringsAsFactors = FALSE
      )

      # --- 散布図用のデータ集積（mRNA_predの時だけ実行） ---
      if (label == "mRNA_pred") {
        pair_list_sp[[as.character(i)]] <- data.frame(
          pred = v_ind,
          obs  = log10(v_exp),
          Species = sp
        )
      }
    }
  }

  # --- ループ終了後、リストをデータフレームに変換 ---
  
  # 1. ジッタープロット用 (all_results)
  if(length(res_list_sp) > 0) {
    all_results_sp <- do.call(rbind, res_list_sp)
    # Q値計算などの後処理
    all_results_sp$Qval <- p.adjust(all_results_sp$Pval, method="fdr")
    all_results_sp$Qval[all_results_sp$Qval == 0] <- 1e-300
    all_results_sp$Score <- -log10(all_results_sp$Qval) * all_results_sp$Rho
    all_results_sp$Species <- sp
    plot_data_list[[sp]] <- all_results_sp
  }

  # 2. 散布図用 (sp_merged_data)
  if(length(pair_list_sp) > 0) {
    sp_merged_data <- do.call(rbind, pair_list_sp)
    if(nrow(sp_merged_data) > 200000) {
      set.seed(123)
      sp_merged_data <- sp_merged_data[sample(nrow(sp_merged_data), 200000), ]
    }
    all_pair_data_list[[sp]] <- sp_merged_data
  }
}

# --- プロット処理 ---
# Scoreが計算されているか確認し、NAを排除
final_df <- do.call(rbind, plot_data_list)
final_df <- final_df[!is.na(final_df$Score), ]

# Speciesのfactor化で、levelsにない名前があるとNAになるため、存在する種名のみに限定
valid_sp <- intersect(names(species_labels), unique(final_df$Species))
final_df$Species <- factor(final_df$Species, levels = valid_sp, labels = species_labels[valid_sp])

# IndexTypeの色指定用ラベルも、実際にデータにあるものに合わせる
final_df$IndexType <- factor(final_df$IndexGroup, 
                             levels = c("GC3", "stAI", "CAI", "mRNA_pred"), 
                             labels = c("GC3 content", "stAI", "CAIhigh", "Pred. logTPM"))

p <- ggplot(final_df, aes(x = IndexType, y = Score, color = IndexType)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1.2) +
  # ここを "free_y" から "free" (または "free_x") に変更
  facet_wrap(~ Species, ncol = 6, scales = "free") + 
  
  scale_color_manual(values = c("GC3 content" = "royalblue", "stAI" = "orchid", "CAIhigh" = "darkorange", "Pred. logTPM" = "darkcyan")) +
  labs(x = "", y = expression(-log[10](Q) %*% r ~ "(vs. observed TPM)")) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_text(family = "sans", face = "plain", size = 12),
    legend.position = "none",
    strip.background = element_blank(), 
    strip.text = element_text(face = "italic"),
    panel.spacing = unit(0.8, "lines") # 軸名が重ならないよう少し広げる
  )

print(p)# 保存
ggsave("correlation_TPM_jitter_facet_v2.png", p, width = 12, height = 20, dpi = 450)


# --- 最後に全生物種を結合して散布図作成 ---
final_scatter_df <- do.call(rbind, all_pair_data_list)
final_scatter_df$Species <- factor(final_scatter_df$Species, 
                                   levels = names(species_labels), 
                                   labels = species_labels)

# 散布図の作成
p_suppl <- ggplot(final_scatter_df, aes(x = pred, y = obs)) +
  geom_bin2d(bins = 80) + 
  scale_fill_viridis_c(option = "viridis", trans = "log10", guide = "none") +
  geom_smooth(method = "lm", color = "magenta", se = FALSE, linewidth = 0.8) +
  stat_cor(method = "pearson", 
           label.x.npc = 0.05,  # 少し内側に寄せる
           label.y.npc = 0.95, 
           size = 2.8,          # 16cm幅に合わせ少し小さく
           color = "black",
           label.sep = "\n",
           output.type = "text",
           vjust = 1) +
  
  # 30種ある場合は nrow = 2 では足りないため、自動計算に任せるか nrow = 5 等を指定
  facet_wrap(~ Species, ncol = 5, scales = "free") +
  
  theme_bw(base_size = 10) + # 16cm幅に収めるためbase_sizeを調整
  labs(x = expression(Predicted~log[10](TPM)), 
       y = expression(Observed~log[10](TPM))) + # 実測TPMに修正
       
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(), 
    strip.text = element_text(face = "italic", size = 9),
    panel.spacing = unit(0.3, "lines"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.4)
  )

# 保存設定
ggsave("supplement_scatter_TPM_multi_panel.png", p_suppl, width = 17, height = 21, units = "cm", dpi = 450)


###########################################
#For the final heatmaps

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

species_list <- c("Arabidopsis_thaliana", "Oryza_sativa", "Homo_sapiens", "Mus_musculus")

calculate_correlation_matrices <- function(species, ind_file, exp_file, out_prefix) {
  message(paste("Processing:", species))
  
  # ファイル読み込み
  indices <- read.table(ind_file, header=TRUE, row.names=1, sep="\t")
  exp_data <- read.table(exp_file, header=TRUE, row.names=1, sep="\t")
  
  # 行列の初期化
  raw_p_mtrx <- matrix(NA, nrow=ncol(indices), ncol=ncol(exp_data))
  est_mtrx <- matrix(NA, nrow=ncol(indices), ncol=ncol(exp_data))
  
  rownames(raw_p_mtrx) <- colnames(indices)
  colnames(raw_p_mtrx) <- colnames(exp_data)
  rownames(est_mtrx) <- colnames(indices)
  colnames(est_mtrx) <- colnames(exp_data)

  # 相関計算ループ
  for (i in 1:ncol(exp_data)) {
    # 【修正】列名に mRNA_HL が入っているかチェック
    is_hl_column <- grepl("mRNA_HL", colnames(exp_data)[i])
    
    for (j in 1:ncol(indices)) {
      # --- 【追加】列名がコドン（A,T,G,Cからなる3文字）か判定 ---
      # 控えめに判定する場合は "^[ATGC]{3}$"、大文字小文字を許容する場合は "[a-zA-Z]{3}"
      is_codon_column <- grepl("^[ATGC]{3}$", colnames(indices)[j])
      
      # 遺伝子の共通集合を取得
      # 対数変換を行うため、値が0より大きい行のみを使用
      genes_used <- intersect(rownames(exp_data[exp_data[,i] > 0, , drop=FALSE]), 
                             rownames(indices[indices[,j] > 0, , drop=FALSE]))
      
      if (length(genes_used) < 20) next
      
      # 発現データ側の対数変換
      e_val <- exp_data[genes_used, i]
      if (!is_hl_column && !grepl("half_lives", exp_file)) { 
        e_val <- log10(e_val) 
      }
      
      # 指標データ側の処理
      i_val <- indices[genes_used, j]
      # --- 【修正】指標がコドンの時だけ対数変換 ---
      if (is_codon_column) {
        i_val <- log10(i_val)
      }
      
      res <- cor.test(e_val, i_val, method="pearson")
      raw_p_mtrx[j, i] <- res$p.value
      est_mtrx[j, i] <- res$estimate
    }
  }
  
  # 生データを保存（後で一括補正するため）
  write.table(raw_p_mtrx, paste0(out_prefix, "_raw_p.txt"), sep="\t", quote=F)
  write.table(est_mtrx, paste0(out_prefix, "_estimates.txt"), sep="\t", quote=F)
}

# 実行ループ
for (sp in species_list) {
  # ファイル名は適宜プロジェクトの命名規則に合わせてください
  #ind_f <- paste0(sp, "_predTPM_indices_codonfreq.txt")
  ind_f <- paste0(sp, "_pred_TPM_indices_sorted.txt")
  exp_f <- paste0(sp, "_TPM_HL_PA.txt") # 統合済みファイル
  if(sp == "Oryza_sativa"){exp_f <- exp_f <- paste0(sp, "_TPM_PA.txt")}
  if(!file.exists(exp_f)) next 
  
  calculate_correlation_matrices(sp, ind_f, exp_f, sp)
}

visualize_cor <- function(spkey, sp_name_clean, codon_order) {
  message(paste("Drawing:", spkey))
  
  raw_p <- as.matrix(read.table(paste0(spkey, "_raw_p.txt"), header=T, row.names=1, sep="\t"))
  est <- as.matrix(read.table(paste0(spkey, "_estimates.txt"), header=T, row.names=1, sep="\t"))
  
  # 1. 一括FDR補正
  q_vals_vec <- p.adjust(raw_p, method="fdr")
  q_mtrx <- matrix(q_vals_vec, nrow=nrow(raw_p), ncol=ncol(raw_p), 
                   dimnames=list(rownames(raw_p), colnames(raw_p)))
  
  # 2. スコア計算
  plot_data <- (-log10(q_mtrx)) * est
  plot_data[is.na(plot_data)] <- 0
  
  # 3. データの並び替え
  target_codons <- intersect(codon_order, rownames(plot_data))
  other_indices <- setdiff(rownames(plot_data), target_codons)
  plot_data_sub <- rbind(plot_data[other_indices, , drop=FALSE],
                         plot_data[target_codons, , drop=FALSE])
  
  rownames(plot_data_sub) <- aa_addition0(rownames(plot_data_sub))
  
  # --- 【追加】サンプルの種類を判別して注釈バーを作成 ---
  col_names <- colnames(plot_data_sub)
  sample_group <- rep("Protein", length(col_names)) # デフォルトをタンパク質量に
  
# 判別ロジックの適用
  sample_group <- sapply(col_names, function(x) {
    if (grepl("SRR|ERR|DRR", x)) {
      return("TPM (mRNA abundance)")
    } else if (grepl("_mRNA_HL", x)) {
      return("mRNA half life")
    } else {
      # どちらにも該当しないものをデフォルトで Protein abundance に
      # ※必要に応じて here に "TPM" 判定を追加することも可能です
      return("Protein abundance")
    }
  })

  # カテゴリーの順序を固定（凡例の並び順に影響します）
  sample_group <- factor(sample_group, 
                         levels = c("TPM (mRNA abundance)", "mRNA half life", "Protein abundance"))

  # 注釈の色設定（表示名に合わせてキーを変更）
  group_colors <- c(
    "TPM (mRNA abundance)" = "orchid1", 
    "mRNA half life"       = "tan1", 
    "Protein abundance"    = "cyan"
  )
  
  # 下部バーの作成 (bottom_annotation)
  column_ha <- HeatmapAnnotation(
    Category = sample_group,
    col = list(Category = group_colors),
    show_legend = TRUE,
    show_annotation_name = FALSE,
    simple_anno_size = unit(2, "mm") 
  )
  # -----------------------------------------------------

  # 4. カラーレンジ
  color_limit <- 100
  col_fun = colorRamp2(
    c(-color_limit, -1.3, 0, 1.3, color_limit), 
    c("royalblue", "#B0C4DE", "white", "#FFDAB9", "darkorange")
  )

  # 5. ComplexHeatmap構築
  ht <- Heatmap(plot_data_sub, 
                name = "-log10(q) * r",
                col = col_fun,
                cluster_rows = FALSE, 
                cluster_columns = FALSE,
                show_row_names = TRUE,
                show_column_names = FALSE,
                bottom_annotation = column_ha, # 【追加】ここで下部バーを合体
                row_names_gp = gpar(fontsize = 8, fontface = "plain"),
                use_raster = FALSE,
                heatmap_legend_param = list(
                  title = "-log10(q) * r",
                  title_position = "topcenter", 
                  direction = "horizontal",
                  legend_width = unit(3, "cm"),
                  title_gp = gpar(fontsize = 10, fontface = "plain"),
                  labels_gp = gpar(fontsize = 9, fontface = "plain"),
                  at = c(-color_limit, 0, color_limit),
                  labels = c(paste0("< -", color_limit), "0", paste0("> ", color_limit))
                ))

  # 6. 出力設定
  png(paste0(spkey, "_codon_weight-expression.png"), width = 2250, height = 5400, res = 450)
  draw(ht, heatmap_legend_side = "top", padding = unit(c(45, 10, 10, 10), "mm"))
  grid.text(paste0(sp_name_clean, ": expression-codon weight"), 
            x = 0.5, y = 0.985, gp = gpar(fontsize = 16, fontface = "plain"))
  dev.off()
}

species_list <- c("Arabidopsis_thaliana", "Oryza_sativa", "Homo_sapiens", "Mus_musculus")

for (sp in species_list) {
  ind_f <- paste0(sp, "_joined_indices_sorted.txt")
  exp_f <- paste0(sp, "_TPM_HL_PA.txt")
  if (sp == "Oryza_sativa"){exp_f <- paste0(sp, "_TPM_PA.txt")}
  
  if(file.exists(ind_f) && file.exists(exp_f)){
    # ステップ1: 計算 (重いので一度やればOK)
    #calculate_and_save_stats(sp, ind_f, exp_f)
    
    # ステップ2: 描画
    visualize_cor(sp, gsub("_", " ", sp), codon_order_wo_stop)
  } else {
    message(paste("Files missing for:", sp))
  }
}

visualize_cor2 <- function(spkey, sp_name_clean, codon_order) {
  message(paste("Drawing:", spkey))
  
  raw_p <- as.matrix(read.table(paste0(spkey, "_raw_p.txt"), header=T, row.names=1, sep="\t"))
  est <- as.matrix(read.table(paste0(spkey, "_estimates.txt"), header=T, row.names=1, sep="\t"))
  
  # 1. 一括FDR補正
  q_vals_vec <- p.adjust(raw_p, method="fdr")
  q_mtrx <- matrix(q_vals_vec, nrow=nrow(raw_p), ncol=ncol(raw_p), 
                   dimnames=list(rownames(raw_p), colnames(raw_p)))
  
  # 2. スコア計算
  plot_data <- (-log10(q_mtrx)) * est
  plot_data[is.na(plot_data)] <- 0
  
  # 3. 行名のラベル書き換え (Model_T0 -> Model (T>0) 形式)
  # 既存の指標名（CAI, GC等）は維持し、予測モデル名のみ置換します
  new_rownames <- rownames(plot_data)
  new_rownames <- gsub("Linear_T0", "Linear (T>0)", new_rownames)
  new_rownames <- gsub("Linear_T1", "Linear (T>1)", new_rownames)
  new_rownames <- gsub("Ridge_T0", "Ridge (T>0)", new_rownames)
  new_rownames <- gsub("Ridge_T1", "Ridge (T>1)", new_rownames)
  new_rownames <- gsub("Lasso_T0", "Lasso (T>0)", new_rownames)
  new_rownames <- gsub("Lasso_T1", "Lasso (T>1)", new_rownames)
  new_rownames <- gsub("Weighted_T0", "Wtd (T>0)", new_rownames)
  new_rownames <- gsub("Weighted_T1", "Wtd (T>1)", new_rownames)
  new_rownames <- gsub("PCR_T0", "PCR (T>0)", new_rownames)
  new_rownames <- gsub("PCR_T1", "PCR (T>1)", new_rownames)
  
  rownames(plot_data) <- new_rownames
  
  # 4. データの並び替え
  # 書き換え後の名前でソートが必要な場合はここを調整
  target_codons <- intersect(codon_order, rownames(plot_data))
  other_indices <- setdiff(rownames(plot_data), target_codons)
  
  # 予測モデルを上の方に、コドンを下の方に配置する
  plot_data_sub <- rbind(plot_data[other_indices, , drop=FALSE],
                         plot_data[target_codons, , drop=FALSE])
  
  # アミノ酸情報の付加（コドン行のみ）
  # ※aa_addition0が内部で文字列完全一致を見る場合、書き換えた後の名前に対応させる必要があります
  rownames(plot_data_sub) <- aa_addition0(rownames(plot_data_sub))
  
  # --- サンプルの種類を判別して注釈バーを作成 (変更なし) ---
  col_names <- colnames(plot_data_sub)
  sample_group <- sapply(col_names, function(x) {
    if (grepl("SRR|ERR|DRR", x)) return("TPM (mRNA abundance)")
    else if (grepl("_mRNA_HL", x)) return("mRNA half life")
    else return("Protein abundance")
  })
  sample_group <- factor(sample_group, levels = c("TPM (mRNA abundance)", "mRNA half life", "Protein abundance"))
  group_colors <- c("TPM (mRNA abundance)" = "orchid1", "mRNA half life" = "tan1", "Protein abundance" = "cyan")
  
  column_ha <- HeatmapAnnotation(
    Category = sample_group,
    col = list(Category = group_colors),
    show_legend = TRUE,
    show_annotation_name = FALSE,
    simple_anno_size = unit(2, "mm") 
  )

  # 5. カラーレンジと描画 (微調整)
  color_limit <- 100
  col_fun = colorRamp2(
    c(-color_limit, -1.3, 0, 1.3, color_limit), 
    c("royalblue", "#B0C4DE", "white", "#FFDAB9", "darkorange")
  )

  ht <- Heatmap(plot_data_sub, 
                name = "-log10(q) * r",
                col = col_fun,
                cluster_rows = FALSE, 
                cluster_columns = FALSE,
                show_row_names = TRUE,
                show_column_names = FALSE,
                bottom_annotation = column_ha,
                row_names_gp = gpar(fontsize = 8, fontface = "plain"),
                # モデル名の行を少し強調したい場合は、ここでフォントスタイルを制御可能
                use_raster = FALSE,
                heatmap_legend_param = list(
                  title = "-log10(q) * r",
                  title_position = "topcenter", 
                  direction = "horizontal",
                  legend_width = unit(3, "cm"),
                  title_gp = gpar(fontsize = 10, fontface = "plain"),
                  labels_gp = gpar(fontsize = 9, fontface = "plain"),
                  at = c(-color_limit, 0, color_limit),
                  labels = c(paste0("< -", color_limit), "0", paste0("> ", color_limit))
                ))

  # 6. 出力
  png(paste0(spkey, "_model_comparison_heatmap.png"), width = 2500, height = 2500, res = 450)
  draw(ht, heatmap_legend_side = "top", padding = unit(c(45, 10, 10, 10), "mm"))
  grid.text(paste0(sp_name_clean, ": Model Prediction & Codon Metrics"), 
            x = 0.5, y = 0.985, gp = gpar(fontsize = 16, fontface = "plain"))
  dev.off()
}

for (sp in species_list) {
  ind_f <- paste0(sp, "_pred_TPM_indices_sorted.txt")
  exp_f <- paste0(sp, "_TPM_HL_PA.txt")
  if (sp == "Oryza_sativa"){exp_f <- paste0(sp, "_TPM_PA.txt")}
  
  if(file.exists(ind_f) && file.exists(exp_f)){
    # ステップ1: 計算 (重いので一度やればOK)
    #calculate_and_save_stats(sp, ind_f, exp_f)
    
    # ステップ2: 描画
    visualize_cor2(sp, gsub("_", " ", sp), codon_order_wo_stop)
  } else {
    message(paste("Files missing for:", sp))
  }
}

