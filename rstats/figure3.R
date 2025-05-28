library(conflicted)
library(tidyverse)

source("rstats/common_settings.R")


# fig3a Pi ---------------------------------------------------------------------

files = list.files("out/pi")
df_pi = dplyr::tibble()

for (f in files) {
  fpath = stringr::str_c("out", "pi", f, sep = "/")
  pop_name = stringr::str_remove(f, ".windowed.pi")
  tmp = readr::read_tsv(fpath) |>
    dplyr::mutate(FID = pop_name)
  df_pi = dplyr::bind_rows(df_pi, tmp)
}

df_pi = df_pi |>
  dplyr::mutate(group = dplyr::case_when(
    stringr::str_detect(FID, "spadiceus|gallus|jabouillei|bankiva") ~ "Wild",
    stringr::str_detect(FID, "RJF") ~ "Zoo",
    stringr::str_detect(FID, "LDH|BLH") ~ "Commercial",
    .default = "Local"
  )) |>
  # dplyr::group_by(group) |>
  dplyr::mutate(FID = forcats::fct_reorder(FID, PI))


ppi = ggplot(df_pi) +
  aes(x = FID, y = PI) +
  geom_violin(aes(fill = group), color = NA, alpha = .66) +
  geom_boxplot(fill = NA, color = "#444444", width = .33, outliers = TRUE) +
  scale_x_discrete(labels = labels) +
  scale_fill_manual(values = groups) +
  labs(y = "Nucleotide diversity") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )


# fig3b Tajima's D -------------------------------------------------------------

files = list.files("out/tajimasD")
df_tD = dplyr::tibble()

for (f in files) {
  fpath = stringr::str_c("out", "tajimasD", f, sep = "/")
  pop_name = stringr::str_remove(f, ".Tajima.D")
  tmp = readr::read_tsv(fpath) |>
    dplyr::mutate(FID = pop_name)
  df_tD = dplyr::bind_rows(df_tD, tmp)
}

df_tD = df_tD |>
  tidyr::drop_na(TajimaD) |>
  dplyr::mutate(group = dplyr::case_when(
    stringr::str_detect(FID, "spadiceus|gallus|jabouillei|bankiva") ~ "Wild",
    stringr::str_detect(FID, "RJF") ~ "Zoo",
    stringr::str_detect(FID, "LDH|BLH") ~ "Commercial",
    .default = "Local"
  )) |>
  # dplyr::group_by(group) |>
  dplyr::mutate(FID = forcats::fct_reorder(FID, TajimaD))

ptD = ggplot(df_tD) +
  aes(x = FID, y = TajimaD) +
  geom_violin(aes(fill = group), color = NA, alpha = .66) +
  geom_boxplot(fill = NA, width = .33, outliers = TRUE) +
  geom_hline(yintercept = 0, color = "#666666", linetype = "dashed") +
  scale_x_discrete(labels = labels) +
  scale_fill_manual(values = groups) +
  labs(y = "Tajima's D") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )


# PSMC -------------------------------------------------------------------------

as_fp = function(x, digits = 1) {
  x = formatC(x, digits = digits, format = "e") |> 
    stringr::str_replace("^(.*)e", "e") |> 
    stringr::str_replace("e", "10^") |> 
    stringr::str_replace('\\+', '')
  return(parse(text = x))
}

df_psmc = readr::read_tsv("out/psmc_result.tsv")

psmc_iter20 = df_psmc |>
  tidyr::separate(Run, into = c("Run", "param"), sep = "\\.") |>
  dplyr::filter(iter == 20) |>
  dplyr::inner_join(sample2group, by = dplyr::join_by(Run == IID))

ppsmc = ggplot(psmc_iter20) +
  aes(x = time, y = Ne) +
  geom_step(aes(group = sample_id, color = group)) +
  scale_x_log10(breaks = c(1e+3, 1e+4, 1e+5, 1e+6, 1e+7), labels = as_fp) +
  scale_y_log10(breaks = c(1e+4, 1e+5, 1e+6), labels = as_fp) +
  scale_color_manual(values = groups) +
  facet_grid(cols = vars(group)) +
  labs(
    x = expression(paste("Years ago (g = 1, µ = 1.91 × ", 10^{-9}, ")")), 
    y = expression(paste("Effective population size"))
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

current_Ne = psmc_iter20 |> dplyr::filter(time == 0)
anova(aov(current_Ne$Ne~current_Ne$group))
pairwise.t.test(x = current_Ne$Ne, g = current_Ne$group, p.adjust.method = "BH")

pNe = ggplot(current_Ne) +
  aes(x = group, y = Ne) +
  geom_jitter(aes(color = group), height = 0, width = .2, size = 2) +
  ggplot2::annotate("segment", x = 1, xend = 2.9, y = 5e+5, yend = 5e+5) +
  ggplot2::annotate("segment", x = 3.1, xend = 4, y = 5e+5, yend = 5e+5) +
  ggplot2::annotate("text", x = 2, y = 6e+5, label = "*") +
  ggplot2::annotate("text", x = 3.5, y = 6e+5, label = "*") +
  scale_y_log10(limits = c(1e+3, 1e+6), breaks = c(1e+3, 1e+4, 1e+5, 1e+6), labels = as_fp) +
  scale_color_manual(values = groups) +
  labs(y = expression(paste("Current ", italic(N)[e]))) +
  theme_test(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

# cowplot ----------------------------------------------------------------------

row3 = cowplot::plot_grid(ppsmc, pNe, ncol = 2, rel_widths = c(3, 1), labels = c("c", "d"), label_size = 20)
f = cowplot::plot_grid(ppi, ptD, row3, nrow = 3, labels = c("a", "b", ""), label_size = 20)
ggsave("images/figure3.png", f, w = 12, h = 10, bg = "#FFFFFF")
