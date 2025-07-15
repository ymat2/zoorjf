library(conflicted)
library(tidyverse)
library(myrrr)

source("./rstats/common_settings.R")


acc2chr = myrrr::grcg7b
chr_levels = c(1:39, "W", "Z")


## PCA -----

colors2 = c(
  "Guangxi" = "#ccebc5",
  "Indonesia" = "#ccebc5",
  "Yunnan" = "#ccebc5",
  "Thailand" = "#ccebc5",
  "Vietnam" = "#ccebc5",
  "bankiva" = "#fbb4ae",
  "gallus" = "#fbb4ae",
  "jabouillei" = "#fbb4ae",
  "spadiceus_China" = "#e41a1c",
  "spadiceus_Thailand" = "#e41a1c",
  "spadiceus_Singapore" = "#e41a1c",
  "RJFizoo" = "#decbe4",
  #"RJFabrc" = "#984ea3",
  "RJFkmt" = "#984ea3",
  "RJFtama" = "#984ea3",
  "RJFmrym" = "#984ea3",
  "LDH" = "#b3cde3",
  "BLH" = "#b3cde3"
)

eigenvec = readr::read_tsv("out/pca/RJF.snp.pca.eigenvec") |>
  dplyr::mutate(IID = stringr::str_split(IID, "/", simplify = TRUE)[,2]) |>
  dplyr::select(!FID) |>
  dplyr::left_join(sample2group, by = "IID") |>
  dplyr::mutate(FID = forcats::fct_relevel(FID, names(colors)))

eigenval = readr::read_csv("out/pca/RJF.snp.pca.eigenval", col_names = "V1")
df = data.frame(pc = 1:nrow(eigenval), eigenval/sum(eigenval)*100)

pca12 = ggplot(eigenvec) +
  aes(PC1, PC2, color = FID, shape = FID) +
  geom_point(size = 3, alpha = 0.9) +
  #ggplot2::annotate("text", label = "3 zooRJF populations", x = .15, y = .15, color = "#984ea3") +
  #ggplot2::annotate("text", label = "Wild RJF populations", x = 0, y = .15, color = "#e41a1c") +
  ggplot2::annotate("segment", x = -.03, y = .05, xend = .12, yend = .08, color = "#333333", 
                    arrow = arrow(length = unit(0.1, "inches"), ends = "both")) +
  ggplot2::annotate("text", label = expression(paste({italic("F")["ST"]})), x = .05, y = .05, color = "#333333", size = 5) +
  scale_color_manual(values = colors2, labels = labels) +
  scale_shape_manual(values = shapes, labels = labels) +
  labs(
    x = paste0("PC1 (", round(df[1,2], digits = 2), "%)"),
    y = paste0("PC2 (", round(df[2,2], digits = 2), "%)"),
    color = "Population", shape = "Population"
  ) +
  theme_test(base_size = 14)
pca12


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
  dplyr::filter(TajimaD < .tajima_low | TajimaD > .tajima_up)

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
  theme_test(base_size = 14) +
  theme(legend.position = "none")
pf2D


## cowplot -----

.pad = cowplot::plot_grid(
  pca12, pf2D, ncol = 2, align = "h", axis = "tb",
  labels = c("a", "d"), label_size = 20
)
p = cowplot::plot_grid(.pad, pfst, ptd, nrow = 3, rel_heights = c(2, 1, 1), labels = c("", "b", "c"), label_size = 20)
p
ggsave("images/figureS1.png", p, h = 10.5, w = 14)
