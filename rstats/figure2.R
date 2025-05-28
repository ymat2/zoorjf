library(conflicted)
library(tidyverse)
library(phytools)
library(ggtree)
library(OptM)
library(myrrr)

source("rstats/common_settings.R")


# fig.2b PCA -------------------------------------------------------------------

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
  scale_color_manual(values = colors, labels = labels) +
  scale_shape_manual(values = shapes, labels = labels) +
  labs(
    x = paste0("PC1 (", round(df[1,2], digits = 2), "%)"),
    y = paste0("PC2 (", round(df[2,2], digits = 2), "%)"),
    color = "Population", shape = "Population"
  ) +
  theme_bw(base_size = 18) +
  theme(legend.position = "none")
pca12

pca23 = ggplot(eigenvec) +
  aes(PC2, PC3, color = FID, shape = FID) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = colors, labels = labels) +
  scale_shape_manual(values = shapes, labels = labels) +
  labs(
    x = paste0("PC2 (", round(df[2,2], digits = 2), "%)"),
    y = paste0("PC3 (", round(df[3,2], digits = 2), "%)"),
    color = "Population", shape = "Population"
  ) +
  theme_bw(base_size = 18) +
  theme()
pca23


# fig1c NJtree + ADMIXTURE -----------------------------------------------------

njtree = ape::read.tree("out/RJF.sub.nj.tree")
njtree_label = data.frame(label = njtree$tip.label) |>
  dplyr::mutate(label = stringr::str_remove_all(label, "'")) |>
  dplyr::mutate(label = stringr::str_split(label, "/", simplify = TRUE)[,2]) |>
  dplyr::left_join(sample2group, by = dplyr::join_by(label == IID)) |>
  dplyr::mutate(label = sample_id)
njtree$tip.label = njtree_label$label
outgroup = stringr::str_subset(njtree$tip.label, "bankiva")

.cladelabs = dplyr::tibble(
  labels = c(labels |> stringr::str_remove("GGb"), "Indonesia"),
  nodes = c(155, 170, 147, 131, 136, 1, 228, 233, 124, 159, 195,201, 210, 222, 219, 185, 178, 174)
)
.intercepts = dplyr::tibble(
  x = 15,
  xend = Inf,
  y = c(2.2, 8.5, 17.5, 25.5, 31.5, 38.5, 47.5, 53.5, 62.5, 65.5, 73.5, 78.5, 86.5, 94.5, 102.5, 107.5, 111.5)
)

t_nj = njtree |>
  ape::root(outgroup = outgroup) |>
  dplyr::full_join(njtree_label, by = "label") |>
  ggtree(linewidth = .5, branch.length = "none") + 
  geom_tiplab(aes(color = group), size = 2) +
  geom_nodelab(aes(label = label), hjust = -0.2, node = "internal", size = 2) +
  #geom_hline(yintercept = .intercepts, linetype = "dashed") +
  geom_segment(data = .intercepts, aes(x = x, xend = xend, y = y, yend = y), linetype = "dashed", linewidth = .25) +
  geom_cladelab(
    node = .cladelabs$nodes, 
    label = .cladelabs$labels, 
    geom = "text", 
    textcolour = "#333333", 
    barcolor = "transparent", 
    size = 5, 
    label.size = 0, 
    offset = 15,
    hjust = 1,
    barsize = .5, 
    align = TRUE
  ) +
  scale_x_continuous(expand = expansion(mult = c(.05, .05))) +
  scale_color_manual(values = groups) +
  theme(legend.position = "none")
t_nj

fam = readr::read_tsv("out/admixture/RJF.snp.fam", col_names = FALSE) |>
  dplyr::mutate(Run = stringr::str_split(X1, "/", simplify = TRUE)[,2]) |>
  dplyr::full_join(sample2group, by = dplyr::join_by(Run == IID)) |>
  dplyr::mutate(sample_id = dplyr::if_else(is.na(sample_id), Run, sample_id)) |>
  tidyr::drop_na(X1) |>
  dplyr::select(Run, FID, sample_id, geo_loc_name)

q2 = readr::read_delim("out/admixture/RJF.snp.2.Q", col_names = FALSE) |>
  dplyr::bind_cols(fam) |>
  dplyr::mutate(K = "K = 2") |>
  tidyr::pivot_longer(dplyr::starts_with("X"), names_to = "anc", values_to = "prop")

q3 = readr::read_delim("out/admixture/RJF.snp.3.Q", col_names = FALSE) |>
  dplyr::bind_cols(fam) |>
  dplyr::mutate(K = "K = 3") |>
  tidyr::pivot_longer(dplyr::starts_with("X"), names_to = "anc", values_to = "prop")

q4 = readr::read_delim("out/admixture/RJF.snp.4.Q", col_names = FALSE) |>
  dplyr::bind_cols(fam) |>
  dplyr::mutate(K = "K = 4") |>
  tidyr::pivot_longer(dplyr::starts_with("X"), names_to = "anc", values_to = "prop")

df = dplyr::bind_rows(q2, q3, q4) |> dplyr::mutate(label = sample_id)
adm = ggplot(df) +
  aes(x = prop, y = label, fill = anc) +
  geom_bar(stat = "identity", position = "fill", alpha = .8) +
  scale_fill_brewer(palette = "Paired") +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  facet_grid(cols = vars(K), scales = "free_x", space = "free_x") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.title = element_blank(),
    #axis.text.x = element_text(size = 5, angle = 75, hjust = 1),
    axis.text.x = element_blank(),
    axis.text.y = element_blank()
  )

padm = adm |> aplot::insert_left(t_nj, width = 1.41) |> aplot::as.patchwork()


# fig.2d,e Supplemental plots -------------------------------------------------- 

## CV Error ----

cve = readr::read_delim("out/admixture/CV-error.txt", col_names = c("c", "e", "K", "error")) |>
  dplyr::mutate(K = stringr::str_extract(K, "[0-9]+")) |>
  ggplot() +
  aes(x = as.integer(K), y = error) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  scale_x_continuous(breaks = seq(1, 10), labels = seq(1, 10)) +
  labs(x = "K", y = "CV Error") +
  theme_test(base_size = 16)


## Opt M in Treemix ----

folder = "out/treemix"
prefix = stringr::str_c(folder, "treemix", sep  ="/")

test.optM = OptM::optM(folder, method = "Evanno")

deltaM = ggplot(test.optM) +
  aes(x = m, y = Deltam) +
  geom_point(color = "#333333", size = 3) +
  geom_line(color = "#333333", linewidth = 1) +
  scale_x_continuous(labels = seq(0, 8), breaks = seq(0, 8)) +
  labs(
    x = expression("m (migration edges)"),
    y = expression(paste({Delta}, "m", sep=""))
  ) +
  theme_test(base_size = 16) 


# fig1.f Treemix ---------------------------------------------------------------

opt_M = 2
stem = paste0(prefix, ".1.", opt_M)
tmobj = myrrr::read_treemix(stem = stem)
tmobj$layout$tips$pop = labels[tmobj$layout$tips$pop]

ptm = myrrr::plot_treemix(tmobj) +
  scale_x_continuous(expand = expansion(mult = c(.05, .25))) +
  theme_treemix(base_size = 14) +
  theme(
    legend.position = "inside",
    legend.justification = c(1, 0),
    legend.background = element_blank()
  )


# cowplot ----------------------------------------------------------------------

.fab = cowplot::plot_grid(pca12, pca23, nrow = 1, rel_widths = c(1, 1.4), labels = c("a", "b"), label_size = 20)
.fde = cowplot::plot_grid(cve, deltaM, labels = c("d", "e"), label_size = 20)
.fdef = cowplot::plot_grid(.fde, ptm, nrow = 2, rel_heights = c(1, 1.73), labels = c("", "f"), label_size = 20)
.fcdef = cowplot::plot_grid(padm, .fdef, ncol = 2, labels = c("c", ""), label_size = 20)
f = cowplot::plot_grid(.fab, .fcdef, nrow = 2, rel_heights = c(1, 1.73))
ggsave("images/figure2.png", f, w = 14, h = 15, bg = "#FFFFFF")
