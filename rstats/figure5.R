library(conflicted)
library(tidyverse)
library(myrrr)

source("./rstats/common_settings.R")


acc2chr = myrrr::grcg7b
chr_levels = c(1:39, "W", "Z")

## fst -----
fst = readr::read_tsv("out/fst/zoo_vs_spadiceus.windowed.weir.fst") |>
  dplyr::full_join(acc2chr, by = "CHROM") |>
  dplyr::mutate(chr = forcats::fct_relevel(chr, chr_levels)) |>
  dplyr::filter(N_VARIANTS > 10) |>
  dplyr::filter(!chr %in% c("Un", "MT", "W", "Z"))

.fst_threshold = mean(fst$MEAN_FST) + 3 * sd(fst$MEAN_FST)

pfst = myrrr::ggman(
  fst,
  chr = "chr",
  bp = "BIN_START",
  p = "MEAN_FST",
  col = c("#6a51a3", "#cbc9e2"),
  ylab = expression(paste({italic("F")[ST]}))
) + geom_hline(yintercept = .fst_threshold, linetype = "dashed")


## Tajima's D -----

tajimaD = readr::read_tsv("out/fst/zoo.Tajima.D") |>
  dplyr::full_join(acc2chr, by = "CHROM") |>
  dplyr::mutate(
    chr = forcats::fct_relevel(chr, chr_levels),
    BIN_START= BIN_START + 1
  ) |>
  dplyr::filter(N_SNPS > 10) |>
  dplyr::filter(!chr %in% c("Un", "MT", "W", "Z"))

ptd = myrrr::ggman(
  tajimaD,
  chr = "chr",
  bp = "BIN_START",
  p = "TajimaD",
  col = c("#6a51a3", "#cbc9e2"),
  ylab = expression(paste("Tajima's ", italic("D"), sep = "")),
  ylim = c(-3, 5)
)

.tajima_up = mean(tajimaD$TajimaD) + 3 * sd(tajimaD$TajimaD)
.tajima_low = mean(tajimaD$TajimaD) - 3 * sd(tajimaD$TajimaD)


## Overlap of Fst and TajimaD, gene annotation -----

fst_tajimaD = tajimaD |> 
  dplyr::select(CHROM, BIN_START, TajimaD) |>
  dplyr::inner_join(fst, by = c("CHROM", "BIN_START")) |>
  dplyr::mutate(sig = dplyr::if_else(MEAN_FST > .fst_threshold, "1", "0"))

fst01_genes = fst_tajimaD |>
  dplyr::filter(sig == "1") |>
  myrrr::annotate(seqnames = "CHROM", start = "BIN_START", end = "BIN_END") |>
  dplyr::filter(!is.na(gene) & !stringr::str_detect(gene, "^LOC|orf")) |>
  dplyr::distinct(gene, .keep_all = TRUE) |>
  dplyr::filter(gene=="TSHR" | TajimaD < .tajima_low | TajimaD > .tajima_up)

pf2D = ggplot(fst_tajimaD) +
  aes(x = TajimaD, y = MEAN_FST) +
  geom_point(aes(color = sig)) +
  geom_hline(yintercept = .fst_threshold, linetype = "dashed") +
  geom_vline(xintercept = c(.tajima_low, .tajima_up), linetype = "dashed") +
  ggplot2::annotate("errorbarh", y = 0, xmin = .tajima_low, xmax =  .tajima_up, height = .025) +
  ggplot2::annotate("label", y = 0, x = 1.5, label = "Mean±3×SD", label.size = 0) + 
  ggrepel::geom_label_repel(data = fst01_genes, aes(label = gene)) +
  labs(x = expression(paste("Tajima's ", italic("D"))), y = expression(paste({italic("F")["ST"]}))) +
  scale_color_manual(values = c("1" = "#333333", "0" = "#DDDDDD")) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")
pf2D


## TSHR genotyping -----

tshr_vcf = readr::read_tsv("out/TSHR.annot.vcf", comment = "##") |>
  tidyr::pivot_longer(10:130, names_to = "IID", values_to = "GT") |>
  dplyr::mutate(
    IID = stringr::str_split(IID, "\\/", simplify = TRUE)[,2],
    GT = stringr::str_split(GT, ":", simplify = TRUE)[,1],
    AC = dplyr::if_else(GT == "./.", 0, 2),
    VAR = stringr::str_extract(INFO, "p\\.\\w{3}\\d+\\w{3}")
  ) |>
  dplyr::mutate(
    AC_ALT = dplyr::case_when(
      GT == "1/1" ~ 2,
      GT == "0/1" ~ 1,
      .default = 0
    )) |>
  dplyr::select(`#CHROM`, POS, REF, ALT, IID, GT, AC, AC_ALT, VAR) |>
  dplyr::left_join(sample2group, by = "IID")

.label_order = c(
  "GGb", "GGg", "GGj", "GGs (China)", "GGs (Thailand)", "GGs (Singapore)",
  "RJF (izoo)", "RJF (Tama)", "RJF (Maruyama)", "RJF (Kumamoto)",
  "Vietnam", "Thailand", "China (Yunnan)", "China (Guangxi)", "Indonesia",
  "WL", "RIR"
)

tshr_af_arg558gly = tshr_vcf |>
  dplyr::group_by(FID, POS) |>
  summarise(AC = sum(AC), ACalt = sum(AC_ALT)) |>
  dplyr::mutate(
    AF = ACalt/AC, 
    POS = forcats::as_factor(POS), 
    FID = forcats::fct_relevel(FID, names(colors))
  ) |>
  dplyr::left_join(tshr_vcf |> 
                     dplyr::distinct(POS, REF, ALT, VAR) |>
                     dplyr::mutate(POS = forcats::as_factor(POS)), 
                   by = "POS") |>
  dplyr::filter(VAR == "p.Arg558Gly") |>
  dplyr::mutate(FID = labels[FID] |> unname() |> forcats::as_factor()) |>
  dplyr::mutate(FID = forcats::fct_relevel(FID, .label_order) |> forcats::fct_rev()) |>
  print()

a558g = ggplot(tshr_af_arg558gly) +
  aes(x = POS, y = FID) +
  geom_tile(aes(fill = AF)) +
  scale_fill_viridis_c(
    option = "E",
    breaks = c(0, .5, 1),
    labels = c("0.0 (Domestic)", "0.5", "1.0 (Wild)")
  ) +
  labs(
    x = "TSHR - Gly558Arg",
    fill = "Allele Frequency"
  ) +
  theme_test(base_size = 14) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank()
  )
a558g


## cowplot -----

.pab = cowplot::plot_grid(
  pfst, ptd, nrow = 2, align = "v", axis = "rl",
  labels = c("a", "b", ""), label_size = 20
  )
.pcd = cowplot::plot_grid(
  pf2D, a558g, ncol = 2, rel_widths = c(2.4, 1), 
  labels = c("c", "d"), label_size = 20
  )
p = cowplot::plot_grid(.pab, .pcd, nrow = 2)
p
ggsave("images/figure5.png", p, h = 10, w = 12)
