library(tidyverse)
library(RColorBrewer)
library(ggnewscale)

par(mfrow = c(1, 1))

setwd("/home/mie/Schreibtisch/Unistuff/Master/Semester2/Forschungsgruppenpraktikum_Steven/Test_Git")

MAFs <- c(0.01, 0.05) #

for (MAF in MAFs) {
  maf_str <- sprintf(".%02d", MAF * 100)
  prefix <- paste0("results/admixture/maf_",MAF,"/output/maf",maf_str,"_minDP20_maxDP100_minQ40_missing.90.splitted_chr.")
  fam <- read.table(paste0("data/bed_PLINK/maf",maf_str,"_minDP20_maxDP100_minQ40_missing.90.splitted_chr.fam"),
                    header = FALSE, stringsAsFactors = FALSE)
  samples <- read.csv("data/ID_data/species_id_map.csv", header = TRUE, sep = ",")
  samples$species[samples$species == "Triticum araraticum"] <- ""
  cols <- c(
    "#377EB8",
    "#E41A1C",
    "#4DAF4A",
    "#984EA3"
  )
  
  
  pdf(paste0("results/admixture/admixture_K4_K3_K2_",MAF,".pdf"), width=10, height=5)
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
                  ylab=paste0("Anc. Prop., K = ", K),
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
}


# CV-Error-Plot
K <- c(2, 3, 4, 5, 6)
CV_error_maf01 <- c(0.38827, 0.37434, 0.36973, 0.36473, 0.36183)
CV_error_maf05 <- c(0.52720, 0.51230, 0.50720, 0.50343, 0.49814)

df <- data.frame(
  K = rep(c(2, 3, 4, 5, 6), 2),
  CV = c(CV_error_maf01, CV_error_maf05),
  MAF = rep(c("MAF = 0.01", "MAF = 0.05"), each = 5)
)

p <- ggplot(df, aes(K, CV, colour = MAF)) +
  geom_line() +
  geom_point() +
  scale_colour_manual(values = c("#377EB8", "#E41A1C")) +
  labs(
    x = "K",
    y = "ADMIXTURE cross-validation error",
    title = ""
  ) +
  theme_bw() +
  theme(
    legend.position = c(0.95, 0.6),   # rechts innen
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = alpha("white", 0.8), colour = "black"),
    legend.title = element_blank()
  )
pdf(paste0("results/admixture/CV_error_K2_to_6.pdf"), 
    width = 3.5, height = 3.5)
print(p)
dev.off()


# Fst-Heatmap
# MAF 0.01
maf01 <- tribble(
  ~i,     ~j,     ~fst,
  "Pop2", "Pop1", 0.291,
  "Pop3", "Pop1", 0.088,
  "Pop3", "Pop2", 0.232
)

# MAF 0.05
maf05 <- tribble(
  ~i,     ~j,     ~fst,
  "Pop1", "Pop2", 0.285,
  "Pop1", "Pop3", 0.092,
  "Pop2", "Pop3", 0.242
)

pops <- c("Pop1","Pop2","Pop3")
grid <- expand_grid(x = pops, y = pops) %>%
  mutate(
    x = factor(x, levels = pops),
    y = factor(y, levels = pops)
  )
idx <- tibble(pop = pops, idx = seq_along(pops))

maf01_plot <- maf01 %>%
  rename(y = i, x = j) %>%
  left_join(idx, by = c("x" = "pop")) %>% rename(ix = idx) %>%
  left_join(idx, by = c("y" = "pop")) %>% rename(iy = idx) %>%
  mutate(
    x = factor(x, levels = pops),
    y = factor(y, levels = pops)
  )

maf05_plot <- maf05 %>%
  rename(y = i, x = j) %>%
  left_join(idx, by = c("x" = "pop")) %>% rename(ix = idx) %>%
  left_join(idx, by = c("y" = "pop")) %>% rename(iy = idx) %>%
  mutate(
    x = factor(x, levels = pops),
    y = factor(y, levels = pops)
  )


p <-ggplot() +
  geom_tile(data = grid, aes(x, y), fill = NA, color = "grey85", linewidth = 0.4) +
  
  # Untere Dreieckshälfte (MAF 0.01) -> blau
  geom_tile(
    data = maf01_plot,
    aes(x = x, y = y, fill = fst)
  ) +
  scale_fill_gradient(
    name = "MAF = 0.01",
    low = "white", high = "#377EB8",
    limits = c(0, 0.3),
    breaks = c(0, 0.1, 0.2, 0.3),
    oob = scales::squish
  ) +
  geom_text(
    data = maf01_plot,
    aes(x = x, y = y, label = sprintf("%.3f", fst)),
    size = 4,
    colour = "black"
  ) +
  ggnewscale::new_scale_fill() +
  
  # Obere Dreieckshälfte (MAF 0.05) -> rot
  geom_tile(
    data = maf05_plot,
    aes(x = x, y = y, fill = fst)
  ) +
  scale_fill_gradient(
    name = "MAF = 0.05",
    low = "white", high = "#E41A1C",
    limits = c(0, 0.3),
    breaks = c(0, 0.1, 0.2, 0.3),
    oob = scales::squish
  ) +
  geom_text(
    data = maf05_plot,
    aes(x = x, y = y, label = sprintf("%.3f", fst)),
    size = 4,
    colour = "black"
  ) +
  
  coord_equal() +
  labs(
    title = expression("ADMIXTURE "~F[ST]~" between populations")
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    legend.position = "right"
  )

pdf(paste0("results/admixture/Fst_3_Pop_Heatmap.pdf"), 
    width = 5, height = 5)
print(p)
dev.off()



