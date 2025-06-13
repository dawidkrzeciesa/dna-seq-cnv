log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")
rlang::global_entrace()

library(tidyverse)

ref_flat = read_tsv(snakemake@input[["table"]]) |>
  group_by(
    external_gene_name,
    ensembl_transcript_id,
    chromosome_name,
    strand,
    transcript_start,
    transcript_end
  ) |>
  summarize(
    cds_start = min(cds_start, na.rm = TRUE),
    cds_end = max(cds_end, na.rm = TRUE),
    num_exons = n(),
    exon_starts = str_flatten(exon_chrom_start, collapse = ","),
    exon_ends = str_flatten(exon_chrom_end, collapse = ",")
  ) |>
  mutate(
    strand = case_match(
      strand,
      1 ~ "+",
      -1 ~ "-"
    ),
    # to mirror the refFlat.txt format, we add trailing commas to
    # lists of exon start and end positions
    across(c(exon_starts, exon_ends), ~ str_c(.x, ",")),
    # if no cds exists for a (non-coding) transcript, we want an NA
    # value that will be written as an empty value below
    cds_start = na_if(cds_start, Inf),
    cds_end = na_if(cds_end, -Inf)
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
  )

write_tsv(
  ref_flat,
  file = snakemake@output[["table"]],
  col_names = FALSE,
  na = ""
)
