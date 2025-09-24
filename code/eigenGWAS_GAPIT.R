library(GAPIT)

setwd("/home/mie/Schreibtisch/Unistuff/Master/Semester2/Forschungsgruppenpraktikum_Steven/Git/rye_project")

myG <- read.table("data/hmp_TASSEL/mygenotypes.hmp.txt",
                  sep="\t", 
                  comment.char = "",
                  header = FALSE)
#alle_str <- toupper(gsub("\\s+", "", myG[,2]))
#alle_str <- gsub("[|/]", "-", alle_str)
#myG[-1,1] <- paste(myG[-1,1], alle_str[-1], sep="_")

myG[,3] <- sub("(?i)^(?:chr)?([1-7])r$", "\\1", myG[,3], perl = TRUE)
myG[,3] <- sub("(?i)^un$", "8", myG[,3], perl = TRUE)
myG[,3] <- sub("(?i)^b$", "9", myG[,3], perl = TRUE)

myY_raw <- read.table("data/pca_PLINK/maf.01_minDP20_maxDP100_minQ40_missing.10.biallelic.eigenvec", head = FALSE)
keep <- c(1, 3:min(5, ncol(myY_raw)))
myY <- myY_raw[, keep, drop = FALSE]
colnames(myY) <- c("Taxa", paste0("PC", seq_len(ncol(myY)-1)))
Y  <- myY

# PC2+PC3 als Kovariaten?
#CV <- myY[, c("Taxa","PC2","PC3")]

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

