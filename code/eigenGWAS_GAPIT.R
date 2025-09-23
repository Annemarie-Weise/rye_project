library(GAPIT)

setwd("/home/mie/Schreibtisch/Unistuff/Master/Semester2/Forschungsgruppenpraktikum_Steven/Git/rye_project")

myG <- read.table("data/hmp_TASSEL/mygenotypes.hmp.txt",
                  sep="\t", 
                  comment.char = "",
                  header = FALSE)
alle_str <- toupper(gsub("\\s+", "", myG[,2]))
alle_str <- gsub("[|/]", "-", alle_str)
myG[-1,1] <- paste(myG[-1,1], alle_str[-1], sep="_")

myG[,3] <- sub("(?i)^(?:chr)?([1-7])r$", "\\1", myG[,3], perl = TRUE)
myG[,3] <- sub("(?i)^un$", "8", myG[,3], perl = TRUE)
myG[,3] <- sub("(?i)^b$", "9", myG[,3], perl = TRUE)

myY_raw <- read.table("data/pca_PLINK/maf.01_minDP20_maxDP100_minQ40_missing.10.splitted.eigenvec", head = FALSE)
keep <- c(1, 3:min(5, ncol(myY_raw)))
myY <- myY_raw[, keep, drop = FALSE]
colnames(myY) <- c("Taxa", paste0("PC", seq_len(ncol(myY)-1)))
Y  <- myY

# PC2+PC3 als Kovariaten?
#CV <- myY[, c("Taxa","PC2","PC3")]

setwd("/home/mie/Schreibtisch/Unistuff/Master/Semester2/Forschungsgruppenpraktikum_Steven/Git/rye_project/results/EigenGWAS")

con <- file(sprintf("gapit_%s.log", format(Sys.time(), "%Y%m%d")), open = "wt")   # "at" zum Anhängen
sink(con, split = TRUE)              # stdout -> Datei (split=TRUE zeigt es zusätzlich in der Konsole)
sink(con, type = "message")          # stderr (warnings/messages) -> Datei
on.exit({ sink(type="message"); sink(); close(con) }, add = TRUE)

geno <- GAPIT:::GAPIT.Genotype(G = myG, SNP.MAF = 0)
GD <- geno$GD

myGAPIT <- GAPIT(
  Y = Y,
  GD = geno$GD,
  GM = geno$GI,
  PCA.total = 0,
  model = "FarmCPU"
)




