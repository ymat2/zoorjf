library(conflicted)
library(tidyverse)
library(phytools)
library(ggtree)


# fig. 1a: Scheme --------------------------------------------------------------

treetext = "(('Wild RJFs':0.1, 'Zoological park RJFs':0.1):0.5, ('Local breeds':0.4, 'Commercial breeds':0.4):0.2):0.4;"
t = ape::read.tree(text = treetext)
t$tip.label = stringr::str_remove_all(t$tip.label, "'")

yscale = c(3, 4, 2, 1, 3, 3, 2)
names(yscale) = c(t$tip.label, c(5, 6, 7))
colors = c(
  "Wild RJFs" = "#e41a1c",
  "Zoological park RJFs" = "#984ea3",
  "Local breeds" = "#4daf4a",
  "Commercial breeds" = "#377eb8"
)

ptree = ggtree(t, layout = "ellipse", yscale = "label", yscale_mapping = yscale, linewidth = 2, color = "#CCCCCC") +
  geom_tiplab(color = "#333333", size = 7, hjust = -.1) +
  geom_tippoint(aes(color = label), size = 7) +
  geom_rootedge(linewidth = 2, color = "#CCCCCC") +
  scale_x_continuous(expand = expansion(mult = c(.05, 2))) +
  scale_y_continuous(expand = expansion(mult = c(.1, .1))) +
  scale_color_manual(values = colors) +
  theme(legend.position = "none")

f1a = cowplot::ggdraw() +
  cowplot::draw_plot(ptree) +
  cowplot::draw_image("images/rjf.png", x = .8, y = .25, w = .15) +
  cowplot::draw_image("images/chicken.png", x = .8, y = -.3, w = .15)

# fig. 1b: Map -----------------------------------------------------------------

world_map = ggplot2::map_data("world") |> 
  dplyr::as_tibble() |> 
  dplyr::filter(region=="Japan")

df4map=dplyr::tibble(
  zoo=c(
    "Maruyama Zoo (n=4)",
    "Tama Zoological Park (n=5)",
    "iZoo (n=8)",
    "Kumamoto City Zoo (n=10)"
  ),
  long=c(141, 139.4, 138, 130.7),
  lat=c(43, 35.6, 34.7, 32.8),
  nudge = c(-1, 1, -1, -1),
  sample_size=c(4, 5, 8, 10)
)

pmap = df4map |>
  ggplot() + 
  aes(long, lat) +
  geom_map(data = world_map, aes(map_id = region), map = world_map, 
           fill = "#CCCCCC", linewidth = 0) +
  geom_point(aes(size = sample_size), color = "#984ea3") +
  geom_label(aes(y = lat+nudge, label = zoo), size = 5, fill = NA, label.size = 0) +
  scale_x_continuous(labels = scales::label_number(suffix = "°E")) +
  scale_y_continuous(limits = c(30, NA), labels = scales::label_number(suffix = "°N")) +
  scale_size_continuous(range = c(5, 8)) +
  scale_color_brewer(palette = "Set1") +
  coord_fixed() +
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    panel.grid.minor = element_blank()
  )
pmap


# cowplot ----------------------------------------------------------------------

f1 = cowplot::plot_grid(f1a, pmap, nrow = 2, rel_heights = c(1, 3), labels = c("a", "b"), label_size = 20)
ggsave("images/figure1.png", f1, w = 8, h = 8, bg = "#FFFFFF")
