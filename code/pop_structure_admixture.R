library(tidyverse)
library(RColorBrewer)
par(mfrow = c(1, 1))

setwd("/home/mie/Schreibtisch/Unistuff/Master/Semester2/Forschungsgruppenpraktikum_Steven/Git/rye_project")

prefix <- "results/admixture/output/maf.01_minDP20_maxDP100_minQ40_missing.10.splitted_chr."
fam <- read.table("data/bed_PLINK/maf.01_minDP20_maxDP100_minQ40_missing.10.splitted.fam",
                  header = FALSE, stringsAsFactors = FALSE)
samples <- read.csv("data/ID_data/species_id_map.csv", header = TRUE, sep = ",")
cols <- c(
  "#377EB8",
  "#E41A1C",
  "#4DAF4A",
  "#984EA3"
)


pdf("results/admixture/admixture_K4_K3_K2.pdf", width=11.69, height=8.27)
par(mfcol=c(3,1),
    mar=c(0.5, 3.0, 0.5, 0.2),   # bottom,left,top,right 
    oma=c(9, 1, 1, 1),           
    mgp=c(1.4, 0.35, 0),         # Achsentitel näher an Plot
    cex.axis=0.8, las=1,
    xaxs="i",
    tcl = -0.2)                   

prep_tbl <- function(K, prefix, fam, samples){
  tbl <- read.table(paste0(prefix, K, ".Q"))
  tbl <- cbind(tbl, fam[,1])
  names(tbl) <- c(paste0("V", seq_len(K)), "Old_ID")
  tbl <- inner_join(tbl, samples[, c("Old_ID","species")], by="Old_ID")
  arrange(tbl, desc(species), across(all_of(paste0("V", 1:K))))
}

plot_K <- function(K, prefix, fam, samples, draw_x_labels=FALSE, shared_vpos=NULL){
  tbl <- prep_tbl(K, prefix, fam, samples)
  mat <- t(as.matrix(tbl[, paste0("V", seq_len(K))]))
  
  bp <- barplot(mat,
                col=cols[1:K],
                names.arg=rep("", ncol(mat)),
                ylab=paste0("Anc. Proportions, K = ", K),
                border=NA, space=0, axes=TRUE)
  runs <- rle(tbl$species); len <- runs$lengths; grp <- runs$values
  centers_idx <- cumsum(len) - (len - 1)/2
  centers_x   <- bp[centers_idx]
  
  if (length(len) > 1) {
    boundaries_idx <- cumsum(len)
    vpos_local <- (bp[boundaries_idx[-length(boundaries_idx)]] + bp[boundaries_idx[-length(boundaries_idx)]+1]) / 2
  } else {
    vpos_local <- numeric(0)
  }
  vpos <- if (!is.null(shared_vpos)) shared_vpos else vpos_local
  if (length(vpos) > 0){
    usr <- par("usr")   # Achsenlimits: usr[3]=ymin, usr[4]=ymax
    segments(x0 = vpos, y0 = usr[3]+ 0.01, x1 = vpos, y1 = usr[4],
             col = "black", lwd = 1, lend = 1)
  }
  if (draw_x_labels) {
    par(xpd=NA)
    text(centers_x,
         par("usr")[3] - 0.03*diff(par("usr")[3:4]), # kleinerer Abstand
         labels = grp,
         srt = 90,   # 45 Drehung
         cex = 0.9, # kleinere Schrift
         font = 3,   # kursiv
         adj = 1)
    par(xpd=FALSE)
  }
  invisible(vpos_local)
}

vpos_ref <- plot_K(4, prefix, fam, samples, draw_x_labels=FALSE)
plot_K(3, prefix, fam, samples, draw_x_labels=FALSE, shared_vpos=vpos_ref)
plot_K(2, prefix, fam, samples, draw_x_labels=TRUE,  shared_vpos=vpos_ref)

dev.off()


