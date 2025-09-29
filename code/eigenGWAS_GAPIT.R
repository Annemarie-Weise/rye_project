library(GAPIT)
library(dplyr)
library(ggplot2)

setwd("/home/mie/Schreibtisch/Unistuff/Master/Semester2/Forschungsgruppenpraktikum_Steven/Git/rye_project")

myG <- read.table("data/hmp_TASSEL/mygenotypes.hmp.txt",
                  sep="\t", 
                  comment.char = "",
                  header = FALSE)

myG[,3] <- sub("(?i)^(?:chr)?([1-7])r$", "\\1", myG[,3], perl = TRUE)
myG[,3] <- sub("(?i)^un$", "8", myG[,3], perl = TRUE)
myG[,3] <- sub("(?i)^b$", "9", myG[,3], perl = TRUE)

myY_raw <- read.table("data/pca_PLINK/maf.01_minDP20_maxDP100_minQ40_missing.10.biallelic.eigenvec", head = FALSE)
keep <- c(1, 3:min(5, ncol(myY_raw)))
myY <- myY_raw[, keep, drop = FALSE]
colnames(myY) <- c("Taxa", paste0("PC", seq_len(ncol(myY)-1)))
Y  <- myY


setwd("/home/mie/Schreibtisch/Unistuff/Master/Semester2/Forschungsgruppenpraktikum_Steven/Git/rye_project/results/EigenGWAS")

con <- file(sprintf("gapit_%s.log", format(Sys.time(), "%Y%m%d")), open = "wt")   # "at" zum Anhängen
sink(con, split = TRUE)              
sink(con, type = "message")          
on.exit({ sink(type="message"); sink(); close(con) }, add = TRUE)

myGAPIT <- GAPIT(
  Y = Y,
  G = myG,
  PCA.total = 0,
  model = "FarmCPU",
  SNP.MAF = 0,
  file.output = TRUE
) # Out-Puts enthalten p = 0 werte nur in den csvs nicht mit angezeigt bei Manhatten-Plots

sink(NULL, type = "message")  
sink(NULL)                    
close(con) 


# eigener Plot mit den p = 0 Makern
setwd("/home/mie/Schreibtisch/Unistuff/Master/Semester2/Forschungsgruppenpraktikum_Steven/Git/rye_project/results/EigenGWAS")

manhattan_EigenGWAS <- function(alpha = 0.01, title, limits, PC = 0) {
  df <- read.csv(paste0("GAPIT.Association.GWAS_Results.FarmCPU.PC",PC,"(NYC).csv"),
                 header = TRUE)
  df <- df %>% mutate(Significant = P.value < alpha/nrow(df))
  df <- df %>% mutate(minusLog10P = -log10(P.value))
  df$minusLog10P[df$P.value == 0] <- 500
  df$Chr <- factor(
    df$Chr,
    levels = as.character(1:9),
    labels = c(paste0("chr", 1:7, "R"), "chrUn", "chrB")
  )
  
  ggplot(df, aes(x = Pos, y = minusLog10P , color = Significant)) +
    geom_point(alpha = 1, size = 0.8) +
    scale_color_manual(values = c("FALSE" = "darkgrey", "TRUE" = "#4DAF4A")) +
    facet_wrap(~Chr, scales = "free_x", ncol = 4) +
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
}


limit <- c(150,500,250)
for(PC in 1:3){
  p <- manhattan_EigenGWAS(PC = PC, limits = c(0,limit[PC]),
                         title = paste0("PC", PC, ", MAF = 0.01, Model = FarmCPU"))
  pdf(paste0("manhatten_PC",PC,"_GAPIT.pdf"), 
      width = 10, height = 5)
  print(p)
  dev.off()
}



# pcs <- paste0("PC", 1:3)
# res <- setNames(vector("list", length(pcs)), pcs)
# 
# for (pc in pcs) {
#   res[[pc]] <- GAPIT(
#     Y = Y[, c("Taxa", pc), drop = FALSE],
#     G = myG,
#     PCA.total = 0,
#     model = "FarmCPU",
#     file.output = FALSE
#   )$GWAS
# }
# PC1 <- res$PC1
# PC2 <- res$PC2

