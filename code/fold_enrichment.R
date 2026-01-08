library(dplyr)
library(tidyverse)
library(readr)
library(stringr)
library(ggplot2)

setwd("/home/mie/Schreibtisch/Unistuff/Master/Semester2/Forschungsgruppenpraktikum_Steven/Test_Git/results/comp_pcadapt_EigenGWAS/distribution_analysis")

read_snps <- function(path) {
  read_csv(path, col_names = FALSE, show_col_types = FALSE) %>%
    transmute(
      CHR = str_extract(X1, "^chr[^_]+"),
      BP  = as.numeric(str_extract(X1, "(?<=_)\\d+"))
    )
}

add_scaled_r <- function(df) {
  n_all_max <- max(df$n_all, na.rm = TRUE)
  df %>%
    mutate(
      r_PC1_scaled = n_PC1 / n_all_max,
      r_PC2_scaled = n_PC2 / n_all_max,
      r_PC3_scaled = n_PC3 / n_all_max
    )
}

calc_window_ratios <- function(all_df, pc1_df, pc2_df, pc3_df, windows_bp, out_prefix) {
  
  count_in_windows <- function(df, win_bp) {
    df %>%
      mutate(BIN_START = floor(BP / win_bp) * win_bp) %>%
      count(CHR, BIN_START, name = "n")
  }
  
  for (win_bp in windows_bp) {
    all_w <- count_in_windows(all_df, win_bp) %>%
      rename(n_all = n)
    
    pc1_w <- count_in_windows(pc1_df, win_bp) %>% rename(n_PC1 = n)
    pc2_w <- count_in_windows(pc2_df, win_bp) %>% rename(n_PC2 = n)
    pc3_w <- count_in_windows(pc3_df, win_bp) %>% rename(n_PC3 = n)
    
    out <- all_w %>%
      left_join(pc1_w, by = c("CHR","BIN_START")) %>%
      left_join(pc2_w, by = c("CHR","BIN_START")) %>%
      left_join(pc3_w, by = c("CHR","BIN_START")) %>%
      mutate(
        n_PC1 = ifelse(is.na(n_PC1), 0L, n_PC1),
        n_PC2 = ifelse(is.na(n_PC2), 0L, n_PC2),
        n_PC3 = ifelse(is.na(n_PC3), 0L, n_PC3),
        r_PC1 = n_PC1 / n_all,
        r_PC2 = n_PC2 / n_all,
        r_PC3 = n_PC3 / n_all,
        BIN_END = BIN_START + win_bp,
        midBP = BIN_START + win_bp/2,
        WindowBp = win_bp
      ) %>%
      select(CHR, BIN_START, BIN_END, midBP, n_all, n_PC1, r_PC1, n_PC2, r_PC2, n_PC3, r_PC3)
    
    out <- add_scaled_r(out)
    win_mb <- win_bp / 1e6
    out_file <- paste0(out_prefix, "_", win_mb, "Mb.csv")
    write_csv(out, out_file)
  }
}

chr_levels <- c(paste0("chr", 1:7, "R"), "chrUn", "chrB")

read_ri <- function(path, win_label) {
  read_csv(path, show_col_types = FALSE) %>%
    mutate(
      CHR    = factor(CHR, levels = chr_levels),
      Window = factor(win_label, levels = c("1 Mbp","25 Mbp")),
      midBP  = as.numeric(midBP)
    ) %>%
    filter(!CHR %in% c("chrB","chrUn")) %>%
    select(
      CHR, BIN_START, BIN_END, midBP, Window,
      n_all,
      n_PC1, n_PC2, n_PC3
    )
}

plot_nPC <- function(df, title = "") {
  
  df_long <- df %>%
    pivot_longer(
      cols = starts_with("n_PC"),
      names_to = "PC", values_to = "n"
    ) %>%
    mutate(
      PC = factor(PC,
                  levels = c("n_PC1","n_PC2","n_PC3"),
                  labels = c("PC1","PC2","PC3")),
      Window = factor(Window, levels = c("1 Mbp","25 Mbp"))
    )
  # Hintergrundbalken für 25 Mbp
  df_bg <- df %>%
    filter(Window == "25 Mbp") %>%
    mutate(width_bp = BIN_END - BIN_START)
  ymax_pc  <- max(df_long$n, na.rm = TRUE)
  ymax_all <- max(df_bg$n_all, na.rm = TRUE)
  bar_frac <- 0.8
  scale_fac <- ifelse(ymax_all > 0, (ymax_pc * bar_frac) / ymax_all, 1)
  
  df_bg <- df_bg %>%
    mutate(n_all_scaled = n_all * scale_fac)
  
  ggplot() +
    geom_col(
      data = df_bg,
      aes(x = midBP, y = n_all_scaled, width = width_bp, fill = "n_all (25 Mbp)"),
      inherit.aes = FALSE,
      colour = NA,
      alpha = 0.6
    ) +
    # 25 Mbp: Linie + Punkte
    geom_line(
      data = df_long %>% filter(Window == "25 Mbp"),
      aes(x = midBP, y = n, group = interaction(PC, CHR), colour = "n_PC (25 Mbp)"),
      linewidth = 0.35
    ) +
    geom_point(
      data = df_long %>% filter(Window == "25 Mbp", n > 0),
      aes(x = midBP, y = n, colour = "n_PC (25 Mbp)"),
      size = 0.6
    ) +
    # 1 Mbp:
    geom_point(
      data = df_long %>% filter(Window == "1 Mbp", n > 0),
      aes(x = midBP, y = n, colour = "n_PC (1 Mbp)"),
      size = 0.45
    ) +
    facet_grid(rows = vars(PC), cols = vars(CHR), scales = "free_x") +
    scale_x_continuous(labels = function(x) x / 1e6) +
    # Farben (Punkte/Linien)
    scale_colour_manual(
      name = "",
      values = c("n_PC (25 Mbp)" = "black", "n_PC (1 Mbp)" = "#E41A1C")
    ) +
    # Füllung (Balken)
    scale_fill_manual(
      name = "",
      values = c("n_all (25 Mbp)" = "grey70")
    ) +
    guides(
      colour = guide_legend(
        override.aes = list(size = 1.5)
      )
    ) +
    scale_y_continuous(
      name = expression(n[PC]),
      sec.axis = sec_axis(~ . / scale_fac, name = expression(n[all]~"(25 Mbp)"))
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      strip.text.x = element_text(size = 9),
      strip.text.y = element_text(size = 9),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8)
    ) +
    labs(
      x = "Genome position (Mbp, window midpoint)",
      title = title
    )
}





windows <- c(1e6, 25e6)

#MAF 0.01
all_001 <- read_snps("maf_0.01/all_0.01.csv")
pc1_001 <- read_snps("maf_0.01/sig_PC1_0.01.csv")
pc2_001 <- read_snps("maf_0.01/sig_PC2_0.01.csv")
pc3_001 <- read_snps("maf_0.01/sig_PC3_0.01.csv")

calc_window_ratios(
  all_df = all_001,
  pc1_df = pc1_001,
  pc2_df = pc2_001,
  pc3_df = pc3_001,
  windows_bp = windows,
  out_prefix = "maf0.01_ri"
)

res_25 <- read_ri("maf0.01_ri_25Mb.csv", "25 Mbp")
res_1  <- read_ri("maf0.01_ri_1Mb.csv",  "1 Mbp")
df <- bind_rows(res_1, res_25)

p <- plot_nPC(df, title = "Numbers of variants per window (MAF = 0.01)")
ggsave("n_maf0.01_PCrows_1vs25Mb.pdf", p, width = 10, height = 5)
print(p)



# MAF 0.05
all_005 <- read_snps("maf_0.05/all_0.05.csv")
pc1_005 <- read_snps("maf_0.05/sig_PC1_0.05.csv")
pc2_005 <- read_snps("maf_0.05/sig_PC2_0.05.csv")
pc3_005 <- read_snps("maf_0.05/sig_PC3_0.05.csv")

calc_window_ratios(
  all_df = all_005,
  pc1_df = pc1_005,
  pc2_df = pc2_005,
  pc3_df = pc3_005,
  windows_bp = windows,
  out_prefix = "maf0.05_ri"
)

res_25 <- read_ri("maf0.05_ri_25Mb.csv", "25 Mbp")
res_1  <- read_ri("maf0.05_ri_1Mb.csv",  "1 Mbp")
df <- bind_rows(res_1, res_25)

p <- plot_nPC(df, title = "Numbers of variants per window (MAF = 0.05)")
ggsave("n_maf0.05_PCrows_1vs25Mb.pdf", p, width = 10, height = 5)
print(p)
