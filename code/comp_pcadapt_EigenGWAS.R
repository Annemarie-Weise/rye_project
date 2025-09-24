library(dplyr)
library(ggplot2)

setwd("/home/mie/Schreibtisch/Unistuff/Master/Semester2/Forschungsgruppenpraktikum_Steven/Git/rye_project")

compare_pcadapt_EigenGWAS <- function(PC, alpha = 0.01) {
  pcadapt <- read.csv(paste0("results/pcadapt/maf_0.01/componentwise/PC",PC,"_results_0.01_componentwise_pcadapt.csv"), 
                      header = TRUE)
  x <- trimws(as.character(pcadapt$CHR))
  x <- sub("(?i)^chr([1-7])r$", "\\1", x, perl = TRUE)
  x <- sub("(?i)^chrun$", "8", x, perl = TRUE)
  x <- sub("(?i)^chrb$", "9", x, perl = TRUE)
  pcadapt$CHR <- as.integer(x)
  
  eigenGWAS <- read.csv(paste0("results/EigenGWAS/GAPIT.Association.GWAS_Results.FarmCPU.PC",PC,"(NYC).csv"), 
                        header = TRUE)
  names(eigenGWAS) <- c("SNP_eigen", "CHR", "BP", "P_eigen", "MAF", "nobs", "Effect", "H.B.P.Value")
  
  merged <- pcadapt %>% inner_join(eigenGWAS[,2:4], by = c("CHR", "BP"))
  merged$minusLog10P_eigen <- -log10(eigenGWAS$P_eigen)
  merged$Significant_eigen <- eigenGWAS$P_eigen < alpha/length(eigenGWAS$minusLog10P_eigen)
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
    geom_point(aes(y = minusLog10P), color = "#E41A1C", alpha = 0.8, size = 0.7) +
    geom_point(aes(y = minusLog10P_eigen), color = "#4DAF4A",  alpha = 0.8, size = 0.7) +
    facet_wrap(~ CHR, scales = "free_x", ncol = 4) +
    theme_bw() +
    labs(
      x = "Genome position (bp)",
      y = expression(-log[10](p)),
      title = paste0("PC",PC,": pcadapt (red) vs. EigenGWAS (green)")
    ) +
    theme(
      strip.text = element_text(size = 9),
      axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 6),
      legend.position = "none"
    )
  
  target_snps <- read.csv("data/ID_data/only_sylvestre_SNPs_maf.01.csv")[[1]]
  sylvestre_df <- significant[significant$SNP %in% target_snps, ]
  
  list(plot = p, sylvestre_df = sylvestre_df, df = merged, significant = significant)
}




for (PC in 1:3) {
  p <- compare_pcadapt_EigenGWAS(PC)
  write.csv(p$sylvestre_df, row.names = FALSE,
            paste0("results/comp_pcadapt_EigenGWAS/sylvestre_only_0.01_sig_PC",PC,".csv"))
  write.csv(p$df, row.names = FALSE,
            paste0("results/comp_pcadapt_EigenGWAS/all_0.01_PC",PC,".csv"))
  write.csv(p$significant, row.names = FALSE,
            paste0("results/comp_pcadapt_EigenGWAS/sig_0.01_PC",PC,".csv"))
  pdf(paste0("results/comp_pcadapt_EigenGWAS/comp_manhatten_PC",PC,".pdf"), 
      width = 11.69, height = 8.27)
  print(p$plot)
  dev.off()
}

