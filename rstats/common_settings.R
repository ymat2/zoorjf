.sample_info = readr::read_tsv("data/sra_accession.tsv") |>
  dplyr::rename(IID = Run) |>
  dplyr::select(IID, sample_id, geo_loc_name)

sample2group = readr::read_tsv("out/pca/RJF.snp.pca.eigenvec") |>
  dplyr::mutate(IID = stringr::str_split(IID, "/", simplify = TRUE)[,2]) |>
  dplyr::left_join(.sample_info, by = "IID") |>
  dplyr::mutate(sample_id = dplyr::if_else(is.na(sample_id), IID, sample_id)) |>
  dplyr::mutate(FID = dplyr::case_when(
    stringr::str_detect(sample_id, "IZ") ~ "RJFizoo",
    stringr::str_detect(sample_id, "TMP") ~ "RJFtama",
    stringr::str_detect(sample_id, "HoU") ~ "RJFmrym",
    stringr::str_detect(sample_id, "RJF[0-9]{4}") ~ "RJFkmt",
    stringr::str_detect(sample_id, "RJF[0-9]{1,2}") ~ "RJFabrc",
    stringr::str_detect(sample_id, "spadiceus") ~ stringr::str_c("spadiceus_", geo_loc_name),
    .default = stringr::str_split(sample_id, "[0-9]+", simplify = TRUE)[,1]
  )) |>
  dplyr::mutate(group = dplyr::case_when(
    stringr::str_detect(FID, "spadiceus|gallus|jabouillei|bankiva") ~ "Wild",
    stringr::str_detect(FID, "RJF") ~ "Zoo",
    stringr::str_detect(FID, "LDH|BLH") ~ "Commercial",
    .default = "Local"
  )) |>
  dplyr::mutate(FID = forcats::fct_relevel(FID, names(colors))) |>
  dplyr::select(!dplyr::starts_with("PC"))

groups = c(
  "Wild" = "#e41a1c",
  "Zoo" = "#984ea3",
  "Local" = "#4daf4a",
  "Commercial" = "#377eb8"
)

colors = c(
  "Guangxi" = "#4daf4a",
  "Indonesia" = "#4daf4a",
  "Yunnan" = "#4daf4a",
  "Thailand" = "#4daf4a",
  "Vietnam" = "#4daf4a",
  "bankiva" = "#e41a1c",
  "gallus" = "#e41a1c",
  "jabouillei" = "#e41a1c",
  "spadiceus_China" = "#e41a1c",
  "spadiceus_Thailand" = "#e41a1c",
  "spadiceus_Singapore" = "#e41a1c",
  "RJFizoo" = "#984ea3",
  #"RJFabrc" = "#984ea3",
  "RJFkmt" = "#984ea3",
  "RJFtama" = "#984ea3",
  "RJFmrym" = "#984ea3",
  "LDH" = "#377eb8",
  "BLH" = "#377eb8"
)

shapes = c(
  "Guangxi" = 3,
  "Indonesia" = 1,
  "Yunnan" = 5,
  "Thailand" = 4,
  "Vietnam" = 2,
  "bankiva" = 1,
  "gallus" = 2,
  "jabouillei" = 5,
  "spadiceus_China" = 3,
  "spadiceus_Thailand" = 4,
  "spadiceus_Singapore" = 6,
  "RJFizoo" = 16,
  #"RJFabrc" = 22,
  "RJFkmt" = 18,
  "RJFtama" = 15,
  "RJFmrym" = 17,
  "LDH" = 3,
  "BLH" = 4
)

labels = c(
  "Guangxi" = "China (Guangxi)",
  "Indonesia" = "Indonesia",
  "Yunnan" = "China (Yunnan)",
  "Thailand" = "Thailand",
  "Vietnam" = "Vietnam",
  "bankiva" = "GGb",
  "gallus" = "GGg",
  "jabouillei" = "GGj",
  "spadiceus_China" = "GGs (China)",
  "spadiceus_Thailand" = "GGs (Thailand)",
  "spadiceus_Singapore" = "GGs (Singapore)",
  "RJFizoo" = "RJF (izoo)",
  #"RJFabrc" = "RJF (ABRC)",
  "RJFkmt" = "RJF (Kumamoto)",
  "RJFtama" = "RJF (Tama)",
  "RJFmrym" = "RJF (Maruyama)",
  "LDH" = "RIR",
  "BLH" = "WL"
)
