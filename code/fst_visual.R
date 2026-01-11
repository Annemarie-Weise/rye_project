library(dplyr)
library(readr)
library(stringr)
library(ggplot2)

setwd("~/rye_project")

normalize_chr <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("(?i)^chr([1-7])r$", "\\1", x, perl = TRUE)
  x <- sub("(?i)^chrun$", "8", x, perl = TRUE)
  x <- sub("(?i)^chrb$", "9", x, perl = TRUE)
  as.integer(x)
}

chr_factor <- function(chr_int) {
  factor(
    chr_int,
    levels = 1:9,
    labels = c(paste0("chr", 1:7, "R"), "chrUn", "chrB")
  )
}

read_target_snps <- function(path) {
  read_csv(path, show_col_types = FALSE) %>%
    mutate(
      CHR = normalize_chr(str_extract(SNP, "^chr[^_]+")),
      BP  = as.numeric(str_extract(SNP, "(?<=_)\\d+"))
    ) %>%
    select(CHR, BP)
}


read_fst_windowed <- function(path) {
  df <- read_tsv(path, show_col_types = FALSE)
  
  maf <- as.numeric(str_extract(path, "maf0\\.\\d{2}") |> str_remove("maf"))
  win <- as.integer(str_extract(path, "\\d+Mb") |> str_remove("Mb"))
  
  df %>%
    transmute(
      MAF = factor(maf),
      WindowMb = factor(win, levels = c(1, 25)),
      CHR = normalize_chr(CHROM),
      BIN_START = BIN_START,
      BIN_END   = BIN_END,
      midBP = (BIN_START + BIN_END) / 2,
      FST = pmax(0, as.numeric(WEIGHTED_FST))   # <- negative auf 0
    ) %>%
    filter(!is.na(CHR), is.finite(FST))
}


plot_fst <- function(fst_df, snp_df, mark = 0.02, title = "") {
  
  fst_df <- fst_df %>%
    mutate(
      CHR_f = chr_factor(CHR),
      WindowMb = factor(WindowMb, levels = c(1, 25),
                        labels = c("1 Mbp", "25 Mbp"))
    ) %>%
    filter(CHR_f != "chrB") %>%
    filter(CHR_f != "chrUn")
  
  tri_df <- fst_df %>%
    filter(MAF == "0.01", WindowMb == "1 Mbp") %>%   # <- NUR MAF 0.01
    inner_join(snp_df, by = "CHR") %>%
    filter(BP >= BIN_START, BP < BIN_END) %>%
    mutate(
      y = FST + 0.07,
      CHR_f = chr_factor(CHR)
    )
  
  ggplot(fst_df, aes(midBP, FST, colour = WindowMb)) +
    # 25 Mbp: Linie
    geom_line(
      data = subset(fst_df, WindowMb == "25 Mbp"),
      linewidth = 0.35
    ) +
    # 1 Mbp: Punkte
    geom_point(
      data = subset(fst_df, WindowMb == "1 Mbp"),
      size = 0.4
    ) +
    # 25 Mbp: zusätzlich Punkte
    geom_point(
      data = subset(fst_df, WindowMb == "25 Mbp"),
      size = 0.4
    )+
    # Target SNPs
    geom_point(
      data = tri_df,
      aes(x = BP, y = y),
      inherit.aes = FALSE,
      shape = 25,
      fill = "darkblue",
      colour = "darkblue",
      size = 1.8
    )+
    facet_grid(
      rows = vars(MAF),
      cols = vars(CHR_f),
      scales = "free_x",
      labeller = labeller(
        MAF = function(x) paste0("MAF = ", x)
      )
    ) +
    
    scale_x_continuous(labels = function(x) x / 1e6) +
    scale_colour_manual(values = c("1 Mbp" = "#E41A1C", "25 Mbp" = "black")) +
    guides(
      colour = guide_legend(
        override.aes = list(
          linewidth = c(0, 2.0),  # 1 Mbp: keine Linie, 25 Mbp: dicke Linie
          size      = c(2.0, 2.0) # große Punkte in der Legende
        )
      )
    )+
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 14),
      strip.text.x = element_text(size = 14),
      strip.text.y.right = element_text(angle = 90, size = 13),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14)
    ) +
    
    labs(
      x = "Genome position (Mbp)",
      y = expression("Weighted "~F[ST]),
      colour = "Window size",
      title = title
    )
}


targets <- read_target_snps("data/ID_data/only_sylvestre_SNPs_maf.01.csv")

# sylvestre vs everything
fst_files <- c(
  "results/vcftools/Fst_sylvestre_other_maf0.01_25Mb.windowed.weir.fst",
  "results/vcftools/Fst_sylvestre_other_maf0.01_1Mb.windowed.weir.fst",
  "results/vcftools/Fst_sylvestre_other_maf0.05_25Mb.windowed.weir.fst",
  "results/vcftools/Fst_sylvestre_other_maf0.05_1Mb.windowed.weir.fst"
)
fst_all <- bind_rows(lapply(fst_files, read_fst_windowed))
p <- plot_fst(fst_all, targets, title = expression(italic("S. sylvestre") ~ " vs. all other individuals")) 
ggsave("results/vcftools/fst_sylvestre_vs_all_other.pdf",
       p, width = 15, height = 4)

# sylvestre vs strictum
fst_files <- c(
  "results/vcftools/Fst_strictum_sylvestre_maf0.01_25Mb.windowed.weir.fst",
  "results/vcftools/Fst_strictum_sylvestre_maf0.01_1Mb.windowed.weir.fst",
  "results/vcftools/Fst_strictum_sylvestre_maf0.05_25Mb.windowed.weir.fst",
  "results/vcftools/Fst_strictum_sylvestre_maf0.05_1Mb.windowed.weir.fst"
)
fst_all <- bind_rows(lapply(fst_files, read_fst_windowed))
p <- plot_fst(fst_all, targets, title = expression(bold("(a)") ~ italic("S. sylvestre") ~ " vs. " ~ italic("S. strictum")) ) 
ggsave("results/vcftools/fst_sylvestre_vs_strictum.pdf",
       p, width = 15, height = 4)

# sylvestre vs other minus strictum
fst_files <- c(
  "results/vcftools/Fst_sylvestre_other_minus_strictum_maf0.01_25Mb.windowed.weir.fst",
  "results/vcftools/Fst_sylvestre_other_minus_strictum_maf0.01_1Mb.windowed.weir.fst",
  "results/vcftools/Fst_sylvestre_other_minus_strictum_maf0.05_25Mb.windowed.weir.fst",
  "results/vcftools/Fst_sylvestre_other_minus_strictum_maf0.05_1Mb.windowed.weir.fst"
)
fst_all <- bind_rows(lapply(fst_files, read_fst_windowed))
p <- plot_fst(fst_all, targets, title = expression(bold("(b)") ~ italic("S. sylvestre") ~ " vs. all other individuals without " ~ italic("S. strictum")))
ggsave("results/vcftools/fst_sylvestre_vs_other_minus_strictum.pdf",
       p, width = 15, height = 4)

