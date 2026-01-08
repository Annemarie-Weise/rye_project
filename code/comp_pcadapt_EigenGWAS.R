library(dplyr)
library(ggplot2)
library(cowplot)

setwd("/home/mie/Schreibtisch/Unistuff/Master/Semester2/Forschungsgruppenpraktikum_Steven/Test_Git")

compare_pcadapt_EigenGWAS <- function(PC, alpha = 0.01, mark = 0.5, MAF = 0.01) {
  pcadapt <- read.csv(paste0("results/pcadapt/maf_",MAF,"/componentwise/PC",PC,"_results_",MAF,"_componentwise_pcadapt.csv"), 
                      header = TRUE)
  x <- trimws(as.character(pcadapt$CHR))
  x <- sub("(?i)^chr([1-7])r$", "\\1", x, perl = TRUE)
  x <- sub("(?i)^chrun$", "8", x, perl = TRUE)
  x <- sub("(?i)^chrb$", "9", x, perl = TRUE)
  pcadapt$CHR <- as.integer(x)
  
  eigenGWAS <- read.csv(paste0("results/EigenGWAS/maf_",MAF,"/GAPIT.Association.GWAS_Results.FarmCPU.PC",PC,"(NYC).csv"), 
                        header = TRUE)
  names(eigenGWAS) <- c("SNP_eigen", "CHR", "BP", "P_eigen", "MAF", "nobs", "Effect", "H.B.P.Value")
  
  merged <- pcadapt %>% inner_join(eigenGWAS[,2:4], by = c("CHR", "BP"))
  merged$minusLog10P_eigen <- -log10(merged$P_eigen)
  merged$Significant_eigen <- merged$P_eigen < alpha/length(merged$minusLog10P_eigen)
  significant <- merged %>% filter(Significant, Significant_eigen)
  significant <- significant %>%
    mutate(across(c(minusLog10P_eigen, minusLog10P),
                  ~ ifelse(is.infinite(.), 500, .)))
  significant$CHR <- factor(
    significant$CHR,
    levels = as.character(1:9),
    labels = c(paste0("chr", 1:7, "R"), "chrUn", "chrB")
  )
  
  p <- ggplot(significant, aes(x = BP)) +
    geom_segment(aes(xend = BP, y = minusLog10P, yend = minusLog10P_eigen),
                 color = "black", linewidth = 0.2, alpha = 0.6) +
    geom_point(aes(y = minusLog10P), color = "#E41A1C", size = 0.7) +
    geom_point(aes(y = minusLog10P_eigen), color = "#4DAF4A", size = 0.7) +
    facet_wrap(~ CHR, scales = "free_x", ncol = 4) +
    theme_bw() +
    labs(
      x = "Genome position (Mbp)",
      y = expression(-log[10](p)),
      title = paste0("PC",PC)
    ) +
    scale_x_continuous(
      labels = function(x) x / 1e6
    ) +
    theme(
      strip.text = element_text(size = 12),
      axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 9),
      legend.position = "none"
    )
  
  target_snps <- read.csv("data/ID_data/only_sylvestre_SNPs_maf.01.csv")[[1]]
  sylvestre_df <- significant[significant$SNP %in% target_snps, ]
  sylvestre_df$CHR <- factor(sylvestre_df$CHR, levels = levels(sylvestre_df$CHR))
  p <- p + 
    geom_point(data = sylvestre_df, 
               aes(x = BP, y = ifelse(minusLog10P > minusLog10P_eigen,
                              minusLog10P + mark,minusLog10P_eigen + mark)),
               inherit.aes = FALSE, shape = 25, fill = "darkblue", size = 1, colour = "darkblue")
  legend_df <- data.frame(
    Method = c("pcadapt", "EigenGWAS"),
    x = 1,
    y = c(2, 1)
  )
  legend_plot <- ggplot(legend_df, aes(x = x, y = y, color = Method)) +
    geom_point(size = 0, alpha = 0) +   # <- im Panel unsichtbar
    scale_color_manual(
      values = c(pcadapt = "#E41A1C", EigenGWAS = "#4DAF4A"),
      name = "Method"
    ) +
    guides(
      color = guide_legend(override.aes = list(size = 2.5, alpha = 1))
    ) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 10, color = "black"),
      legend.title = element_text(size = 10, color = "black"),
      legend.key.height = unit(0.6, "cm")
    )
  
  genes <- read.csv("data/ID_data/important_genes.csv", header = TRUE, stringsAsFactors = FALSE)
  genes$BP_Start <- as.numeric(genes$BP_Start)
  genes$BP_End   <- as.numeric(genes$BP_End)
  flank <- 5000
  sig_chr <- as.character(significant$CHR)
  sig_bp  <- significant$BP
  
  significant$Full_match_ids <- mapply(function(ch, bp){
    ids <- genes$ID[ genes$CHR == ch & genes$BP_Start <= bp & bp <= genes$BP_End ]
    if (length(ids)) paste(ids, collapse = ";") else NA
  }, sig_chr, sig_bp)
  significant$Border_match_ids <- mapply(function(ch, bp){
    in_left  <- genes$ID[ genes$CHR == ch & (genes$BP_Start - flank) <= bp & bp <  genes$BP_Start ]
    in_right <- genes$ID[ genes$CHR == ch &  genes$BP_End   <  bp & bp <= (genes$BP_End + flank) ]
    ids <- c(in_left, in_right)
    if (length(ids)) paste(ids, collapse = ";") else NA
  }, sig_chr, sig_bp)
  
  genes$Num_Inside_SNP          <- 0L
  genes$Val_highest_SNP         <- NA
  genes$Pos_highest_SNP         <- NA
  genes$Num_Border_SNP          <- 0L
  genes$Val_highest_SNP_border  <- NA
  genes$Pos_highest_SNP_border  <- NA
  
  id2row <- setNames(seq_len(nrow(genes)), genes$ID)
  best_stat_per_snp <- pmax(significant$minusLog10P,
                            significant$minusLog10P_eigen,
                            na.rm = TRUE)
  
  split_ids <- function(x) {
    if (is.na(x)) return(character(0))
    y <- unlist(strsplit(x, ";", fixed = TRUE))
    y[nzchar(y)]
  }
  update_max <- function(cur_val, cur_pos, cand_val, cand_pos) {
    if (is.na(cand_val) || (!is.na(cur_val) && cand_val <= cur_val)) return(list(val = cur_val, pos = cur_pos))
    else {
      return(list(val = cand_val, pos = cand_pos))
    } 
  }
  
  for (k in seq_len(nrow(significant))) {
    cand_val <- best_stat_per_snp[k]
    cand_pos <- significant$BP[k]
    ids_in <- split_ids(significant$Full_match_ids[k])
    if (length(ids_in)) {
      for (gid in ids_in) {
        i <- id2row[[gid]]
        genes$Num_Inside_SNP[i] <- genes$Num_Inside_SNP[i] + 1L
        up <- update_max(genes$Val_highest_SNP[i], genes$Pos_highest_SNP[i],
                          cand_val, cand_pos)
        genes$Val_highest_SNP[i] <- up$val
        genes$Pos_highest_SNP[i] <- up$pos
      }
    }
    ids_bd <- split_ids(significant$Border_match_ids[k])
    if (length(ids_bd)) {
      for (gid in ids_bd) {
        i <- id2row[[gid]]
        genes$Num_Border_SNP[i] <- genes$Num_Border_SNP[i] + 1L
        upb <- update_max(genes$Val_highest_SNP_border[i], genes$Pos_highest_SNP_border[i],
                          cand_val, cand_pos)
        genes$Val_highest_SNP_border[i] <- upb$val
        genes$Pos_highest_SNP_border[i] <- upb$pos
      }
    }
  }
  
  tri_df <- subset(genes, !is.na(Pos_highest_SNP) & Num_Inside_SNP > 0)
  tri_df$CHR <- factor(tri_df$CHR, levels = levels(significant$CHR))
  p <- p + 
    geom_point(data = tri_df, aes(x = Pos_highest_SNP, y = Val_highest_SNP + mark),
                      inherit.aes = FALSE, shape = 25, fill = "black", size = 1) +
    geom_text(
      data = tri_df, angle = 45, hjust = 0, vjust = 0,
      aes(x = Pos_highest_SNP, y = Val_highest_SNP + mark + (0.5*mark), label = Name),
      inherit.aes = FALSE, size = 3, show.legend = FALSE
    )
  
  p <- ggdraw() +
    draw_plot(p, x = 0, y = 0, width = 1, height = 1) +   # ← volle Fläche!
    draw_plot(legend_plot, x = 0.7, y = 0.2, width = 0.25, height = 0.25)
  
  list(plot = p,
       sylvestre_df = sylvestre_df,
       df = merged,
       significant = significant,
       genes = genes)
}


# MAF 0.01
mark <- c(60,70,6)
for (PC in 1:3) {
  p <- compare_pcadapt_EigenGWAS(PC, mark = mark[PC])
  write.csv(p$sylvestre_df, row.names = FALSE,
            paste0("results/comp_pcadapt_EigenGWAS/maf_0.01/sylvestre_only_0.01_sig_PC",PC,".csv"))
  write.csv(p$df, row.names = FALSE,
            paste0("results/comp_pcadapt_EigenGWAS/maf_0.01/all_0.01_PC",PC,".csv"))
  write.csv(p$significant, row.names = FALSE,
            paste0("results/comp_pcadapt_EigenGWAS/maf_0.01/sig_0.01_PC",PC,".csv"))
  write.csv(p$genes, row.names = FALSE,
            paste0("results/comp_pcadapt_EigenGWAS/maf_0.01/gene_test_0.01_PC",PC,".csv"))
  pdf(paste0("results/comp_pcadapt_EigenGWAS/maf_0.01/comp_0.01_manhatten_PC",PC,".pdf"), 
      width = 9, height = 5)
  print(p$plot)
  dev.off()
}

# MAF 0.05
mark <- c(15,7,6)
for (PC in 1:3) {
  p <- compare_pcadapt_EigenGWAS(PC, mark = mark[PC], MAF = 0.05)
  write.csv(p$sylvestre_df, row.names = FALSE,
            paste0("results/comp_pcadapt_EigenGWAS/maf_0.05/sylvestre_only_0.05_sig_PC",PC,".csv"))
  write.csv(p$df, row.names = FALSE,
            paste0("results/comp_pcadapt_EigenGWAS/maf_0.05/all_0.05_PC",PC,".csv"))
  write.csv(p$significant, row.names = FALSE,
            paste0("results/comp_pcadapt_EigenGWAS/maf_0.05/sig_0.05_PC",PC,".csv"))
  write.csv(p$genes, row.names = FALSE,
            paste0("results/comp_pcadapt_EigenGWAS/maf_0.05/gene_test_0.05_PC",PC,".csv"))
  pdf(paste0("results/comp_pcadapt_EigenGWAS/maf_0.05/comp_0.05_manhatten_PC",PC,".pdf"), 
      width = 9, height = 5)
  print(p$plot)
  dev.off()
}

