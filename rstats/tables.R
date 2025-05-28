library(conflicted)
library(tidyverse)

source("rstats/common_settings.R")

.flagstat = readr::read_tsv("out/flag_summary.tsv")
mapping_stat = readr::read_tsv("out/mapping_summary.tsv") |>
  dplyr::select(sample, meandp, meancov) |>
  dplyr::full_join(.flagstat, by = "sample")

sheet01 = dplyr::tibble(
  "Sheet" = c("tableS1", "tableS2", "tableS3"),
  "Description" = c(
    "Sample information of newly sequenced samples",
    "Sample information of additional sequences",
    "Genes in high Fst regions"
))

## Table S1. Sample information of newly sequenced samples ---------------------

zoo_samples = sample2group |>
  dplyr::filter(group == "Zoo") |>
  dplyr::mutate(sex = dplyr::case_when(
    stringr::str_detect(IID, "M$") ~ "Male",
    IID %in% c("RJF5903", "RJF5910", "RJF5911") ~ "Male",
    .default = "Female"
  )) |>
  dplyr::mutate(zoo = dplyr::case_when(
    FID == "RJFmrym" ~ "Maruyama Zoo",
    FID == "RJFizoo" ~ "iZoo",
    FID == "RJFkmt" ~ "Kumamoto City Zoo",
    FID == "RJFtama" ~ "Tama Zoological Park"
  )) |>
  dplyr::select(sample_id, sex, zoo) |>
  dplyr::left_join(mapping_stat, by = dplyr::join_by(sample_id == sample)) |>
  dplyr::rename(
    "mean_read_depth" = meandp,
    "average_coverage" = meancov,
    "number_of_reads" = num_reads,
    "proportion_of_properly_mapped_reads" = prop_mapped
  )


## Table S2. Sample information of additional sequences ------------------------

ext_samples = sample2group |>
  dplyr::filter(group != "Zoo") |>
  dplyr::left_join(mapping_stat, by = dplyr::join_by(IID == sample)) |>
  dplyr::select(sample_id, FID, group, meandp, meancov, num_reads, prop_mapped, IID) |>
  dplyr::rename(
    "population" = FID,
    "mean_read_depth" = meandp,
    "average_coverage" = meancov,
    "number_of_reads" = num_reads,
    "proportion_of_properly_mapped_reads" = prop_mapped,
    "run" = IID
  )

## Table S3. Gene list of Fst and Tajima'D -------------------------------------

fst = readr::read_tsv("out/fst/zoo_vs_spadiceus.windowed.weir.fst") |>
  dplyr::full_join(acc2chr, by = "CHROM") |>
  dplyr::mutate(chr = forcats::fct_relevel(chr, chr_levels)) |>
  dplyr::filter(N_VARIANTS > 10) |>
  dplyr::filter(!chr %in% c("Un", "MT", "W", "Z"))
.fst_threshold = mean(fst$MEAN_FST) + 3 * sd(fst$MEAN_FST)

tajimaD = readr::read_tsv("out/fst/zoo.Tajima.D") |>
  dplyr::full_join(acc2chr, by = "CHROM") |>
  dplyr::mutate(
    chr = forcats::fct_relevel(chr, chr_levels),
    BIN_START= BIN_START + 1
  ) |>
  dplyr::filter(N_SNPS > 10) |>
  dplyr::filter(!chr %in% c("Un", "MT", "W", "Z"))
.tajima_up = mean(tajimaD$TajimaD) + 3 * sd(tajimaD$TajimaD)
.tajima_low = mean(tajimaD$TajimaD) - 3 * sd(tajimaD$TajimaD)

fst_tajimaD_genes = tajimaD |> 
  dplyr::select(CHROM, BIN_START, TajimaD) |>
  dplyr::inner_join(fst, by = c("CHROM", "BIN_START")) |>
  dplyr::filter(MEAN_FST > .fst_threshold) |>
  myrrr::annotate(seqnames = "CHROM", start = "BIN_START", end = "BIN_END") |>
  dplyr::filter(!is.na(gene) & !stringr::str_detect(gene, "^LOC")) |>
  dplyr::distinct(gene, .keep_all = TRUE) |>
  dplyr::rename(Fst = MEAN_FST) |>
  dplyr::select(CHROM, chr, BIN_START, BIN_END, TajimaD, Fst, gene)


## Write xlsx ------------------------------------------------------------------

writexl::write_xlsx(
  list("Sheet01" = sheet01, "tableS1" = zoo_samples, "tableS2" = ext_samples, "tableS3" = fst_tajimaD_genes),
  path = "docs/supplementary_tables.xlsx",
  format_headers = FALSE
)
