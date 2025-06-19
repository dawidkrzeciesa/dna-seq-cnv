log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")
rlang::global_entrace()

library(multidplyr)
library(tidyverse)

multi_core_cluster <- new_cluster(snakemake@threads)

partitioned_input = read_tsv(snakemake@input[["table"]]) |>
  arrange(
    exon_chrom_start
  ) |>
  group_by(
    external_gene_name,
    ensembl_transcript_id,
    chromosome_name,
    strand,
    transcript_start,
    transcript_end
  ) |>
  partition(
    multi_core_cluster
  )

ref_flat <- partitioned_input |>
  mutate(
    # make cds positions absolute on the chromosome, vs. relative to the transcript
    across(c(cds_start, cds_end), ~ .x + exon_chrom_start - 1)
  ) |>
  summarize(
    cds_start = min(cds_start, na.rm = TRUE),
    cds_end = max(cds_end, na.rm = TRUE),
    num_exons = dplyr::n(),
    exon_starts = stringr::str_flatten(exon_chrom_start, collapse = ","),
    exon_ends = stringr::str_flatten(exon_chrom_end, collapse = ",")
  ) |>
  mutate(
    strand = dplyr::case_match(
      strand,
      1 ~ "+",
      -1 ~ "-"
    ),
    # to mirror the refFlat.txt format, we add trailing commas to
    # lists of exon start and end positions
    across(c(exon_starts, exon_ends), ~ stringr::str_c(.x, ",")),
    # always put something into the cds start and end:
    # * the min() and max() return Inf / -Inf for entries with only NA
    # * refFlat.txt files seem to have the transcript_end position in
    # both (without something in there, the cnvkit batch format sniffing
    # regex for this fails)
    across(
      c(cds_start, cds_end),
      ~ dplyr::if_else(is.infinite(.x), transcript_end, .x)
    )
  ) |>
  select(
    external_gene_name,
    ensembl_transcript_id,
    chromosome_name,
    strand,
    transcript_start,
    transcript_end,
    cds_start,
    cds_end,
    num_exons,
    exon_starts,
    exon_ends
  ) |>
  collect()

write_tsv(
  ref_flat,
  file = snakemake@output[["table"]],
  col_names = FALSE,
  na = ""
)
