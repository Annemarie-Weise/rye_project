library(vcfR)
library(pcadapt)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(grid)
library(ggtext)
library(dplyr)
library(tibble)

setwd("/home/mie/Schreibtisch/Unistuff/Master/Semester2/Forschungsgruppenpraktikum_Steven/Git/rye_project")

pcadapt_data <- read.pcadapt("data/bed_PLINK/maf.01_minDP20_maxDP100_minQ40_missing.10.splitted.bed", 
                             type = "bed")

############ MAF 0.01 ################

# ------------------ Mahalanobis -----------------------
MAF <- 0.01
method <- "mahalanobis"



# Screeplot - 20 PCs
K <- 20
result <- pcadapt(pcadapt_data, K = K, min.maf = MAF ,method = method)

pcadapt_screeplot <- function(result, maf, method) {
  prop_var <- (result$singular.values^2)
  df <- data.frame(
    PC = seq_along(prop_var),
    proportion = prop_var
  )
  p <- ggplot(df, aes(x = PC, y = proportion)) +
    geom_line(color = "#377EB8") +
    geom_point(color = "#377EB8") +
    labs(
      x = "Principal component",
      y = "Variance proportion",
      title = paste0("MAF = ", maf, ", Method = ", method)
    ) +
    theme_classic(base_size = 13) +
    theme(
      axis.line         = element_line(color = "black"),
      axis.ticks        = element_line(color = "black"),
      axis.ticks.length = unit(2, "pt"),
      panel.grid        = element_blank()
    )
  return(p)
}

p <- pcadapt_screeplot(result, maf = MAF, method = method)
pdf(paste0("results/pcadapt/maf_",MAF,"/",method,"/screeplot_",MAF,"_",method,"_pcadapt.pdf"), 
    width = 6, height = 4)
print(p)
dev.off()



# PCA-Plot - 3 PCs
K <- 3
free_deg <- K
result <- pcadapt(pcadapt_data, K = K, min.maf = MAF ,method = method )

pcadapt_pc_scatter <- function(result, maf, method, legend_pos = c(0.99, 0.8)) {
  prop_var <- (result$singular.values^2)
  id_map <- read.csv("data/ID_data/species_id_map.csv", header = TRUE)
  lab_pc1 <- paste0("PC1 (", round(prop_var[1], 3), ")")
  lab_pc2 <- paste0("PC2 (", round(prop_var[2], 3), ")")
  lab_pc3 <- paste0("PC3 (", round(prop_var[3], 3), ")")
  
  fam <- read.table("data/bed_PLINK/maf.01_minDP20_maxDP100_minQ40_missing.10.splitted.fam",
                    header = FALSE, stringsAsFactors = FALSE)
  colnames(fam) <- c("FID","IID","PAT","MAT","SEX","PHENO")
  ids <- fam$IID
  
  colnames(result$scores) <- paste0("PC", seq_len(ncol(result$scores)))
  scores <- as_tibble(result$scores) %>%
    mutate(Old_ID = ids, .before = 1)
  scores <- scores %>%
    left_join(id_map %>% select(Old_ID, species), by = "Old_ID") %>%
    filter(!is.na(species), !is.na(PC1), !is.na(PC2), !is.na(PC3)) %>%
    mutate(species = factor(species))
  
  coloures <- c(
    "#4DAF4A",
    "darkblue",
    "#377EB8",
    "#984EA3",
    "#F781BF",
    "#E41A1C",
    "darkred",
    "#FF7F00"
  )
  base_opts <- list(
    scale_color_manual(values = coloures, drop = TRUE, na.translate = FALSE),
    theme_classic(base_size = 12),
    theme(
      legend.text   = element_text(face = "italic"),
      legend.title  = element_text(face = "bold"),
      plot.margin   = margin(4, 2, 2, 2),
      panel.spacing = grid::unit(1, "mm"),
      axis.line         = element_line(color = "black"),
      axis.ticks        = element_line(color = "black"),
      axis.ticks.length = grid::unit(2, "pt"),
      panel.grid        = element_blank()
    )
  )
  
  p12 <- ggplot(scores, aes(PC1, PC2, color = species)) +
    geom_point(alpha = 0.9, size = 1.7) +
    labs(
      x = lab_pc1,
      y = lab_pc2,
      title = paste0("<b>(a)</b> MAF = ", maf, ", Method = ", method),
      color = "Species"
    ) +
    base_opts +
    theme(
      legend.position = "none",
      plot.title = element_markdown(size = 14)
    )
  
  p13 <- ggplot(scores, aes(PC1, PC3, color = species)) +
    geom_point(alpha = 0.9, size = 1.7) +
    labs(
      x = lab_pc1,
      y = lab_pc3,
      title = paste0("<b>(b)</b> MAF = ", maf, ", Method = ", method),
      color = "Species"
    ) +
    base_opts +
    theme(
      legend.position      = legend_pos,
      legend.justification = c(1, 0.5),
      legend.background    = element_rect(fill = scales::alpha("white", 0.7), 
                                          colour = NA),
      legend.key.size      = grid::unit(4, "mm"),
      plot.title           = element_markdown(size = 14)
    ) +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))
  
  p <- plot_grid(p12, p13, ncol = 2, align = "hv")
  return(p)
}

p <- pcadapt_pc_scatter(result, maf = MAF, method = method)
pdf(paste0("results/pcadapt/maf_",MAF,"/",method,"/pca_",MAF,"_",method,"_pcadapt.pdf"), 
    width = 9.27, height = 4.8)
print(p)
dev.off()



# QQ-Plot - 3 PCs
# Vergleich zwischen erwarteten (unter Nullhypothese) und beobachteten P-Werten
log10p <- pchisq(result$chi2.stat, df = free_deg, lower.tail = FALSE, log.p = TRUE) / log(10)

# Testen ob das mit p-Werten passt
# p_res <- result$pvalues
# idx <- which(!is.na(p_res) & p_res > 0 & is.finite(log10p))
# all.equal(log10(p_res[idx]), log10p[idx], tolerance = 1e-12)

qq_log10p <- function(log10p, title = "", limits = c(0, 7000)) {
  obs <- sort(-log10p[is.finite(log10p)])   # beobachtet: -log10(p) (aufsteigend)
  n   <- length(obs)
  exp <- -log10(1 - ppoints(n))            # theoretische Quantile von -log10(U)
  
  ggplot(data.frame(exp, obs), aes(exp, obs)) +
    geom_point(size = 0.6, colour = "#E41A1C") +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    labs(x = expression(Expected ~ -log[10](p)), 
         y = expression(Observed ~ -log[10](p)), title = title) +
    theme_classic(base_size = 13) +
    theme(
      axis.line        = element_line(color = "black"),
      axis.ticks       = element_line(color = "black"),
      axis.ticks.length= unit(2, "pt"),
      panel.grid       = element_blank()
    ) +
    coord_cartesian(ylim = limits)
}

p <- qq_log10p(log10p, title = paste0("MAF = ", MAF, ", Method = ", method))
pdf(paste0("results/pcadapt/maf_",MAF,"/",method,"/qq_",MAF,"_",method,"_pcadapt.pdf"), 
    width = 4.8, height = 4.8)
print(p)
dev.off()



# Histogramm - 3 PCs
pdf(paste0("results/pcadapt/maf_",MAF,"/",method,"/histogram_",MAF,"_",method,"_pcadapt.pdf"), 
    width = 4.8, height = 4.8)
hist(result$pvalues, xlab = "p-values", breaks = 50, col = "#FF7F00", main = "")
title(main = paste0("MAF = ", MAF, ", Method = ", method), font.main = 1, 
      cex.main  = 0.9, adj = 0)
dev.off()


# Loading - LD-Thinning? - 3 PCs
pdf(paste0("results/pcadapt/maf_", MAF, "/", method, "/LD-Thinning_", MAF, "_", method, "_pcadapt.pdf"),
    width = 8, height = 6)
op <- par(mfrow = c(2, 2), mar   = c(4, 4, 1, 1), oma   = c(1, 1, 3, 1))
on.exit(par(op), add = TRUE)
for (i in 1:3) {
  plot(result$loadings[, i],
       pch = 19, cex = .3,
       ylab = paste0("Loadings PC", i), xlab = "")
}
mtext(paste0("   MAF = ", MAF, ", Method = ", method),
      side = 3, outer = TRUE, line = 1, adj = 0, font = 1, cex  = 1.0)
dev.off()




# Loadings-Violin - 3 PCs
loadings_violin_pcadapt <- function(result, title = "") {
  loadings_df <- as.data.frame(result$loadings[, 1:3, drop = FALSE])
  colnames(loadings_df) <- paste0("PC", seq_len(ncol(loadings_df)))
  loadings_long <- pivot_longer(
    loadings_df,
    cols = starts_with("PC"),
    names_to = "PC",
    values_to = "Loading",
    values_drop_na = TRUE
  )
  p <- ggplot(loadings_long, aes(x = PC, y = Loading)) +
    geom_violin(fill = "#E41A1C", alpha = 0.7,
                         draw_quantiles = c(0.25, 0.5, 0.75)) +
    scale_x_discrete(labels = c("PC1", "PC2", "PC3")) +
    labs(x = "Principal component",
                  y = "SNP loading",
                  title = title) +
    theme_classic(base_size = 13) +
    theme(
      axis.line         = element_line(color = "black"),
      axis.ticks        = element_line(color = "black"),
      axis.ticks.length = unit(2, "pt"),
      panel.grid        = element_blank(),
      legend.position   = c(0.75, 0.75)
    )
  return(p)
}
p <- loadings_violin_pcadapt(result, title = paste0("MAF = ", MAF, ", Method = ", method))
pdf(paste0("results/pcadapt/maf_",MAF,"/",method,"/loadings_violins_",MAF,"_",method,"_pcadapt.pdf"), 
    width = 9.27, height = 5)
print(p)
dev.off()



# Loadings von sylvestre SNPs überprüfen - 3 PCs
pcadapt_loadings_with_ids <- function(result, pcs = 1:3, free_deg = free_deg) {
  target_snps <- read.csv("data/ID_data/only_sylvestre_SNPs_maf.01.csv")[[1]]
  bim <- read.table("data/bed_PLINK/maf.01_minDP20_maxDP100_minQ40_missing.10.splitted.bim",
                    header = FALSE, stringsAsFactors = FALSE)
  colnames(bim) <- c("CHR","SNP","CM","BP","A1","A2")
  ids <- paste0(bim$CHR, "_", bim$BP)
  out <- as_tibble(result$loadings[, pcs, drop = FALSE])
  colnames(out) <- paste0("PC", pcs)
  out <- out %>% mutate(SNP = ids, A1 = bim$A1, A2 = bim$A2, .before = 1)
  
  if (is.matrix(result$pvalues)) {
    pv <- as.data.frame(result$pvalues[, pcs, drop = FALSE])
    colnames(pv) <- paste0("p_PC", pcs)
    out <- bind_cols(out, pv)
    log10p_mat <- pchisq(result$chi2.stat[, pcs, drop = FALSE],df = 1, 
                         lower.tail = FALSE, log.p = TRUE) / log(10)
    log10p_mat <- as.data.frame(log10p_mat)
    colnames(log10p_mat) <- paste0("log10p_PC", pcs)
    out <- bind_cols(out, log10p_mat)
  } else {
    out$pvalue <- result$pvalues
    out$log10p <- as.numeric(pchisq(result$chi2.stat, df = free_deg, 
                                    lower.tail = FALSE, log.p = TRUE) / log(10))
  }
  out <- filter(out, SNP %in% target_snps)
  out
}
sylvestre_loadings_df <- pcadapt_loadings_with_ids(result, free_deg = free_deg)
write.csv(sylvestre_loadings_df, 
          paste0("results/pcadapt/maf_",MAF,"/",method,"/pca_loadings_only_sylvestre_",MAF,"_",method,"_pcadapt.csv"), 
          row.names = FALSE)



# Manhatten-Plot - 3 PCs
manhattan_pcadapt <- function(result, alpha = 0.01, title, free_deg, limits = c(0, 7000), PC = 0) {
  bim <- read.table("data/bed_PLINK/maf.01_minDP20_maxDP100_minQ40_missing.10.splitted.bim",
                    header = FALSE, stringsAsFactors = FALSE)
  colnames(bim) <- c("CHR","SNP","CM","BP","A1","A2")
  log10p <- pchisq(if(PC == 0) result$chi2.stat else result$chi2.stat[,PC, drop = TRUE], 
                   df = free_deg, lower.tail = FALSE, log.p = TRUE) / log(10)
  df <- tibble(
    CHR = bim$CHR,
    BP  = bim$BP,
    SNP = paste0(bim$CHR, "_", bim$BP),
    A1 = bim$A1,
    A2 = bim$A2,
    log10p = log10p
  ) %>%
    filter(is.finite(log10p)) %>%
    mutate(minusLog10P = -log10p)
  n <- nrow(df)
  df <- df %>% mutate(Significant = minusLog10P > (-log10(alpha) + log10(n)  ))
  
  p <- ggplot(df, aes(x = BP, y = minusLog10P, color = Significant)) +
    geom_point(alpha = 0.6, size = 0.8) +
    scale_color_manual(values = c("FALSE" = "darkgrey", "TRUE" = "#E41A1C")) +
    facet_wrap(~CHR, scales = "free_x", ncol = 4) +
    labs(
      x = "Genom position (bp)",
      y = expression(-log[10](p)),
      title = title,
      color = paste0("Significance\n(alpha = ", alpha, ",\n Bonferroni)")
    ) +
    coord_cartesian(ylim = limits) +
    theme_bw() +
    theme(
      strip.text = element_text(size = 9),
      axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 6),
      legend.position = c(0.9, 0.1),
      legend.justification = c(1, 0.5)
    )
  list(plot = p, df = df)
}

p <- manhattan_pcadapt(result, free_deg = free_deg,
                      title = paste0("MAF = ", MAF, ", Method = ", method))
write.csv(p$df, row.names = FALSE,
          paste0("results/pcadapt/maf_",MAF,"/",method,"/results_",MAF,"_",method,"_pcadapt.csv"))
pdf(paste0("results/pcadapt/maf_",MAF,"/",method,"/manhatten_",MAF,"_",method,"_pcadapt.pdf"), 
    width = 11.69, height = 8.27)
print(p$plot)
dev.off()



# ------------------ Componentwise -----------------------
MAF <- 0.01
method <- "componentwise"
K <- 3
free_deg <- 1

result <- pcadapt(pcadapt_data, K = K, min.maf = MAF , method = method)

# QQ-Plots
log10p_PC1 <- pchisq(result$chi2.stat[,1], df = free_deg, lower.tail = FALSE, log.p = TRUE) / log(10)
log10p_PC2 <- pchisq(result$chi2.stat[,2], df = free_deg, lower.tail = FALSE, log.p = TRUE) / log(10)
log10p_PC3 <- pchisq(result$chi2.stat[,3], df = free_deg, lower.tail = FALSE, log.p = TRUE) / log(10)

p1 <- qq_log10p(log10p_PC1, title = bquote(atop(bold("(a)") ~ " PC1, MAF = " ~ .(MAF) ~ ",",
                                                "                  Method = " ~ .(method))), 
                limits = c(0, 2000))
p2 <- qq_log10p(log10p_PC2, title = bquote(atop(bold("(b)") ~ " PC2, MAF = " ~ .(MAF) ~ ",",
                                                "                  Method = " ~ .(method))), 
                limits = c(0, 3000))
p3 <- qq_log10p(log10p_PC3, title = bquote(atop(bold("(c)") ~ " PC3, MAF = " ~ .(MAF) ~ ",",
                                                "                  Method = " ~ .(method))), 
                limits = c(0, 80))
pg <- plot_grid(
  p1, p2,
  p3, NULL,       
  ncol = 2, nrow = 2, align = "hv"
)
pdf(paste0("results/pcadapt/maf_", MAF, "/", method, "/qqs_", MAF, "_", method, "_pcadapt.pdf"),
    width = 9.6, height = 9.6)
print(pg)
dev.off()



# Histogramme
pdf(paste0("results/pcadapt/maf_",MAF,"/",method,"/histograms_",MAF,"_",method,"_pcadapt.pdf"), 
           width = 9.6, height = 9.6)
op <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
on.exit(par(op), add = TRUE)
for (k in 1:3) {
  x <- result$pvalues[, k]
  x <- x[is.finite(x)]
  hist(x, xlab = "p-values", breaks = 50, col = "#FF7F00", main = "")
  lab <- paste0("(", letters[k], ")")
  title(
    main = bquote(
      bold(.(lab)) ~ " PC" * .(k) * "," ~ " MAF = " ~ .(MAF) * ", Method = " ~ .(method)),
    adj = 0,      # linksbündig
    cex.main = 0.9
  )
}
plot.new()
dev.off()



# Manhatten-Plot - 3 PCs
limit <- c(2000,2500,80)
for (pc in seq_len(K)) {
  p <- manhattan_pcadapt(result, free_deg = free_deg, limits = c(0, limit[pc]), PC = pc,
    title = paste0("PC", pc, ", MAF = ", MAF, ", Method = ", method)
  )
  write.csv(p$df, row.names = FALSE,
            paste0("results/pcadapt/maf_", MAF, "/", method,
                   "/PC", pc, "_results_", MAF, "_", method, "_pcadapt.csv"))
  pdf(paste0("results/pcadapt/maf_", MAF, "/", method,
             "/PC", pc, "_manhatten_", MAF, "_", method, "_pcadapt.pdf"),
      width = 11.69, height = 8.27)
  print(p$plot)
  dev.off()
}






############ MAF 0.02 ################


# ------------------ Mahalanobis -----------------------
MAF <- 0.02
method <- "mahalanobis"

# Screeplot - 20 PCs
K <- 20

result <- pcadapt(pcadapt_data, K = K, min.maf = MAF ,method = method)
p <- pcadapt_screeplot(result, maf = MAF, method = method)
pdf(paste0("results/pcadapt/maf_",MAF,"/",method,"/screeplot_",MAF,"_",method,"_pcadapt.pdf"), 
    width = 6, height = 4)
print(p)
dev.off()



# PCA-Plot - 3 PCs
K <- 3
free_deg <- K
result <- pcadapt(pcadapt_data, K = K, min.maf = MAF ,method = method )
p <- pcadapt_pc_scatter(result, maf = MAF, method = method)
pdf(paste0("results/pcadapt/maf_",MAF,"/",method,"/pca_",MAF,"_",method,"_pcadapt.pdf"), 
    width = 9.27, height = 4.8)
print(p)
dev.off()



# QQ-Plot - 3 PCs
log10p <- pchisq(result$chi2.stat, df = K, lower.tail = FALSE, log.p = TRUE) / log(10)
p <- qq_log10p(log10p, title = paste0("MAF = ", MAF, ", Method = ", method),
               limits = c(0, 300))
pdf(paste0("results/pcadapt/maf_",MAF,"/",method,"/qq_",MAF,"_",method,"_pcadapt.pdf"), 
    width = 4.8, height = 4.8)
print(p)
dev.off()



# Histogramm - 3 PCs
pdf(paste0("results/pcadapt/maf_",MAF,"/",method,"/histogram_",MAF,"_",method,"_pcadapt.pdf"), 
    width = 4.8, height = 4.8)
hist(result$pvalues, xlab = "p-values", breaks = 50, col = "#FF7F00", main = "")
title(main = paste0("MAF = ", MAF, ", Method = ", method), font.main = 1, 
      cex.main  = 0.9, adj = 0)
dev.off()



# Loading - LD-Thinning? - 3 PCs
pdf(paste0("results/pcadapt/maf_", MAF, "/", method, "/LD-Thinning_", MAF, "_", method, "_pcadapt.pdf"),
    width = 8, height = 6)
op <- par(mfrow = c(2, 2), mar   = c(4, 4, 1, 1), oma   = c(1, 1, 3, 1))
on.exit(par(op), add = TRUE)
for (i in 1:3) {
  plot(result$loadings[, i],
       pch = 19, cex = .3,
       ylab = paste0("Loadings PC", i), xlab = "")
}
mtext(paste0("   MAF = ", MAF, ", Method = ", method),
      side = 3, outer = TRUE, line = 1, adj = 0, font = 1, cex  = 1.0)
dev.off()



# Loadings-Violin - 3 PCs
p <- loadings_violin_pcadapt(result, title = paste0("MAF = ", MAF, ", Method = ", method))
pdf(paste0("results/pcadapt/maf_",MAF,"/",method,"/loadings_violins_",MAF,"_",method,"_pcadapt.pdf"), 
    width = 9.27, height = 5)
print(p)
dev.off()


# Manhatten-Plot - 3 PCs
p <- manhattan_pcadapt(result, free_deg = free_deg,
                       title = paste0("MAF = ", MAF, ", Method = ", method),
                       limits = c(0, 300))
write.csv(p$df, row.names = FALSE,
          paste0("results/pcadapt/maf_",MAF,"/",method,"/results_",MAF,"_",method,"_pcadapt.csv"))
pdf(paste0("results/pcadapt/maf_",MAF,"/",method,"/manhatten_",MAF,"_",method,"_pcadapt.pdf"), 
    width = 11.69, height = 8.27)
print(p$plot)
dev.off()


# ------------------ Componentwise -----------------------
MAF <- 0.02
method <- "componentwise"
K <- 3
free_deg <- 1

result <- pcadapt(pcadapt_data, K = K, min.maf = MAF , method = method)

# QQ-Plots
log10p_PC1 <- pchisq(result$chi2.stat[,1], df = free_deg, lower.tail = FALSE, log.p = TRUE) / log(10)
log10p_PC2 <- pchisq(result$chi2.stat[,2], df = free_deg, lower.tail = FALSE, log.p = TRUE) / log(10)
log10p_PC3 <- pchisq(result$chi2.stat[,3], df = free_deg, lower.tail = FALSE, log.p = TRUE) / log(10)

p1 <- qq_log10p(log10p_PC1, title = bquote(atop(bold("(a)") ~ " PC1, MAF = " ~ .(MAF) ~ ",",
                                                "                  Method = " ~ .(method))), 
                limits = c(0, 200))
p2 <- qq_log10p(log10p_PC2, title = bquote(atop(bold("(b)") ~ " PC2, MAF = " ~ .(MAF) ~ ",",
                                                "                  Method = " ~ .(method))), 
                limits = c(0, 25))
p3 <- qq_log10p(log10p_PC3, title = bquote(atop(bold("(c)") ~ " PC3, MAF = " ~ .(MAF) ~ ",",
                                                "                  Method = " ~ .(method))), 
                limits = c(0, 45))
pg <- plot_grid(
  p1, p2,
  p3, NULL,       
  ncol = 2, nrow = 2, align = "hv"
)
pdf(paste0("results/pcadapt/maf_", MAF, "/", method, "/qqs_", MAF, "_", method, "_pcadapt.pdf"),
    width = 9.6, height = 9.6)
print(pg)
dev.off()



# Histogramme
pdf(paste0("results/pcadapt/maf_",MAF,"/",method,"/histograms_",MAF,"_",method,"_pcadapt.pdf"), 
    width = 9.6, height = 9.6)
op <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
on.exit(par(op), add = TRUE)
for (k in 1:3) {
  x <- result$pvalues[, k]
  x <- x[is.finite(x)]
  hist(x, xlab = "p-values", breaks = 50, col = "#FF7F00", main = "")
  lab <- paste0("(", letters[k], ")")
  title(
    main = bquote(
      bold(.(lab)) ~ " PC" * .(k) * "," ~ " MAF = " ~ .(MAF) * ", Method = " ~ .(method)),
    adj = 0,      # linksbündig
    cex.main = 0.9
  )
}
plot.new()
dev.off()



# Manhatten-Plot - 3 PCs
limit <- c(200,25,30)
for (pc in seq_len(K)) {
  p <- manhattan_pcadapt(result, free_deg = free_deg, limits = c(0, limit[pc]), PC = pc,
                         title = paste0("PC", pc, ", MAF = ", MAF, ", Method = ", method)
  )
  write.csv(p$df, row.names = FALSE,
            paste0("results/pcadapt/maf_", MAF, "/", method,
                   "/PC", pc, "_results_", MAF, "_", method, "_pcadapt.csv"))
  pdf(paste0("results/pcadapt/maf_", MAF, "/", method,
             "/PC", pc, "_manhatten_", MAF, "_", method, "_pcadapt.pdf"),
      width = 11.69, height = 8.27)
  print(p$plot)
  dev.off()
}



############ MAF 0.03 ################


# ------------------ Mahalanobis -----------------------
MAF <- 0.03
method <- "mahalanobis"

# Screeplot - 20 PCs
K <- 20

result <- pcadapt(pcadapt_data, K = K, min.maf = MAF ,method = method)
p <- pcadapt_screeplot(result, maf = MAF, method = method)
pdf(paste0("results/pcadapt/maf_",MAF,"/",method,"/screeplot_",MAF,"_",method,"_pcadapt.pdf"), 
    width = 6, height = 4)
print(p)
dev.off()



# PCA-Plot - 3 PCs
K <- 3
free_deg <- K
result <- pcadapt(pcadapt_data, K = K, min.maf = MAF ,method = method )
p <- pcadapt_pc_scatter(result, maf = MAF, method = method)
pdf(paste0("results/pcadapt/maf_",MAF,"/",method,"/pca_",MAF,"_",method,"_pcadapt.pdf"), 
    width = 9.27, height = 4.8)
print(p)
dev.off()



# QQ-Plot - 3 PCs
log10p <- pchisq(result$chi2.stat, df = K, lower.tail = FALSE, log.p = TRUE) / log(10)
p <- qq_log10p(log10p, title = paste0("MAF = ", MAF, ", Method = ", method),
               limits = c(0, 400))
pdf(paste0("results/pcadapt/maf_",MAF,"/",method,"/qq_",MAF,"_",method,"_pcadapt.pdf"), 
    width = 4.8, height = 4.8)
print(p)
dev.off()



# Histogramm - 3 PCs
pdf(paste0("results/pcadapt/maf_",MAF,"/",method,"/histogram_",MAF,"_",method,"_pcadapt.pdf"), 
    width = 4.8, height = 4.8)
hist(result$pvalues, xlab = "p-values", breaks = 50, col = "#FF7F00", main = "")
title(main = paste0("MAF = ", MAF, ", Method = ", method), font.main = 1, 
      cex.main  = 0.9, adj = 0)
dev.off()



# Loading - LD-Thinning? - 3 PCs
pdf(paste0("results/pcadapt/maf_", MAF, "/", method, "/LD-Thinning_", MAF, "_", method, "_pcadapt.pdf"),
    width = 8, height = 6)
op <- par(mfrow = c(2, 2), mar   = c(4, 4, 1, 1), oma   = c(1, 1, 3, 1))
on.exit(par(op), add = TRUE)
for (i in 1:3) {
  plot(result$loadings[, i],
       pch = 19, cex = .3,
       ylab = paste0("Loadings PC", i), xlab = "")
}
mtext(paste0("   MAF = ", MAF, ", Method = ", method),
      side = 3, outer = TRUE, line = 1, adj = 0, font = 1, cex  = 1.0)
dev.off()



# Loadings-Violin - 3 PCs
p <- loadings_violin_pcadapt(result, title = paste0("MAF = ", MAF, ", Method = ", method))
pdf(paste0("results/pcadapt/maf_",MAF,"/",method,"/loadings_violins_",MAF,"_",method,"_pcadapt.pdf"), 
    width = 9.27, height = 5)
print(p)
dev.off()


# Manhatten-Plot - 3 PCs
p <- manhattan_pcadapt(result, free_deg = free_deg,
                       title = paste0("MAF = ", MAF, ", Method = ", method),
                       limits = c(0, 400))
write.csv(p$df, row.names = FALSE,
          paste0("results/pcadapt/maf_",MAF,"/",method,"/results_",MAF,"_",method,"_pcadapt.csv"))
pdf(paste0("results/pcadapt/maf_",MAF,"/",method,"/manhatten_",MAF,"_",method,"_pcadapt.pdf"), 
    width = 11.69, height = 8.27)
print(p$plot)
dev.off()


# ------------------ Componentwise -----------------------
MAF <- 0.03
method <- "componentwise"
K <- 3
free_deg <- 1

result <- pcadapt(pcadapt_data, K = K, min.maf = MAF , method = method)

# QQ-Plots
log10p_PC1 <- pchisq(result$chi2.stat[,1], df = free_deg, lower.tail = FALSE, log.p = TRUE) / log(10)
log10p_PC2 <- pchisq(result$chi2.stat[,2], df = free_deg, lower.tail = FALSE, log.p = TRUE) / log(10)
log10p_PC3 <- pchisq(result$chi2.stat[,3], df = free_deg, lower.tail = FALSE, log.p = TRUE) / log(10)

p1 <- qq_log10p(log10p_PC1, title = bquote(atop(bold("(a)") ~ " PC1, MAF = " ~ .(MAF) ~ ",",
                                                "                  Method = " ~ .(method))), 
                limits = c(0, 400))
p2 <- qq_log10p(log10p_PC2, title = bquote(atop(bold("(b)") ~ " PC2, MAF = " ~ .(MAF) ~ ",",
                                                "                  Method = " ~ .(method))), 
                limits = c(0, 70))
p3 <- qq_log10p(log10p_PC3, title = bquote(atop(bold("(c)") ~ " PC3, MAF = " ~ .(MAF) ~ ",",
                                                "                  Method = " ~ .(method))), 
                limits = c(0, 150))
pg <- plot_grid(
  p1, p2,
  p3, NULL,       
  ncol = 2, nrow = 2, align = "hv"
)
pdf(paste0("results/pcadapt/maf_", MAF, "/", method, "/qqs_", MAF, "_", method, "_pcadapt.pdf"),
    width = 9.6, height = 9.6)
print(pg)
dev.off()



# Histogramme
pdf(paste0("results/pcadapt/maf_",MAF,"/",method,"/histograms_",MAF,"_",method,"_pcadapt.pdf"), 
    width = 9.6, height = 9.6)
op <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
on.exit(par(op), add = TRUE)
for (k in 1:3) {
  x <- result$pvalues[, k]
  x <- x[is.finite(x)]
  hist(x, xlab = "p-values", breaks = 50, col = "#FF7F00", main = "")
  lab <- paste0("(", letters[k], ")")
  title(
    main = bquote(
      bold(.(lab)) ~ " PC" * .(k) * "," ~ " MAF = " ~ .(MAF) * ", Method = " ~ .(method)),
    adj = 0,      # linksbündig
    cex.main = 0.9
  )
}
plot.new()
dev.off()



# Manhatten-Plot - 3 PCs
limit <- c(400,75,150)
for (pc in seq_len(K)) {
  p <- manhattan_pcadapt(result, free_deg = free_deg, limits = c(0, limit[pc]), PC = pc,
                         title = paste0("PC", pc, ", MAF = ", MAF, ", Method = ", method)
  )
  write.csv(p$df, row.names = FALSE,
            paste0("results/pcadapt/maf_", MAF, "/", method,
                   "/PC", pc, "_results_", MAF, "_", method, "_pcadapt.csv"))
  pdf(paste0("results/pcadapt/maf_", MAF, "/", method,
             "/PC", pc, "_manhatten_", MAF, "_", method, "_pcadapt.pdf"),
      width = 11.69, height = 8.27)
  print(p$plot)
  dev.off()
}







