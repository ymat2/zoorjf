library(conflicted)
library(tidyverse)
library(detectRUNS)

source("rstats/common_settings.R")


# fig4a heterozygosity ---------------------------------------------------------

het = readr::read_table("out/roh/RJF.snp.het") |>
  dplyr::mutate(IID = stringr::str_split(IID, "/", simplify = TRUE)[,2]) |>
  dplyr::select(!FID) |>
  dplyr::left_join(sample2group, by = "IID") |>
  #dplyr::group_by(group) |>
  dplyr::arrange(group, `F`) |>
  dplyr::mutate(FID = forcats::fct_inorder(FID))

phet = ggplot(het) +
  aes(x = FID, y = `F`) +
  geom_boxplot(aes(fill = group), alpha = .66, outliers = FALSE) +
  geom_jitter(color = "#CCCCCC", size = 1, height = 0, width = .2) +
  scale_x_discrete(labels = labels) +
  scale_fill_manual(values = groups) +
  labs(y = "Inbreeding coefficient") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

# fig4b ROH distribution -------------------------------------------------------

hom_indiv = readr::read_table("out/roh/RJF.snp.roh.hom.indiv") |>
  dplyr::mutate(IID = stringr::str_split(IID, "/", simplify = TRUE)[,2]) |>
  dplyr::select(!FID) |>
  dplyr::left_join(sample2group, by = "IID") |>
  dplyr::mutate(FID = forcats::fct_relevel(FID, names(colors)))

phom = ggplot(hom_indiv) +
  aes(x = NSEG, y = KB/1000) +
  geom_point(aes(color = FID, shape = FID), size = 3, alpha = 0.9) +
  scale_color_manual(values = colors, labels = labels) +
  scale_shape_manual(values = shapes, labels = labels) +
  labs(
    x = "Number of ROH segments", 
    y = "Total length of ROHs (Mbp)",
    color = "Population", shape = "Population"
  ) +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank())


# fig.3c Per Length ROH count --------------------------------------------------

acc2chr = myrrr::grcg7b
chr_levels = c(1:39)

runs = detectRUNS::readExternalRuns("out/roh/RJF.snp.roh.hom", program = "plink") |>
  dplyr::mutate(id = stringr::str_split(id, "/", simplify = TRUE)[,2]) |>
  dplyr::select(!group) |>
  dplyr::left_join(sample2group, by = dplyr::join_by(id == IID)) |>
  dplyr::left_join(acc2chr, by = dplyr::join_by(chrom == CHROM)) |>
  #dplyr::mutate(chrom = forcats::fct_relevel(chr, chr_levels)) |>
  dplyr::mutate(chrom = chr)

length_class = runs |>
  dplyr::mutate(length_class = dplyr::case_when(
    lengthBps < 1000000 ~ "Short (0.3Mb-1Mb)",
    lengthBps < 10000000 ~ "Medium (1Mb-10Mb)",
    .default = "Long (>10Mb)"
  )) |>
  dplyr::count(id, length_class) |>
  tidyr::complete(id, length_class, fill = list(n = 0)) |>
  dplyr::full_join(runs |> dplyr::distinct(id, group), by = "id")

pl = ggplot(length_class) +
  aes(x = group, y = n) +
  geom_jitter(aes(color = length_class), height = 0, width = .2, size = 2) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  scale_color_viridis_d(option = "cividis", begin = .9, end = .1) +
  labs(x = "", y = "Number of ROH segments") +
  facet_wrap(vars(length_class), nrow = 3, scales = "free") +
  theme_test(base_size = 14) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = NA, color = NA),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )


# cowplot ----------------------------------------------------------------------

.fab = cowplot::plot_grid(phet, phom, nrow = 2, labels = c("a", "b"), rel_heights = c(2, 3), label_size = 20)
f = cowplot::plot_grid(.fab, pl, ncol = 2, labels = c("", "c"), rel_widths = c(2, 1), label_size = 20)
ggsave("images/figure4.png", f, w = 10, h = 8.5, bg = "#FFFFFF")
