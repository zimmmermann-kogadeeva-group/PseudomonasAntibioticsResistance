#!/usr/bin/env Rscript

library(tidyverse)
library(DESeq2)
library(tximport)
library(here)

# path_to_metadata <- "Data/PublicDatasets/Khaledi/phenotypes.txt"
df_meta <- read_tsv(snakemake@input[["metadata"]]) %>%
  rename(all_of(c(final_all = "isolate", `Meropenem_S-vs-R` = "mem"))) %>%
  select(isolate, mem) %>%
  drop_na() %>%
  mutate(mem = as.factor(mem))

# path_to_data <- "Data/PublicDatasets/Khaledi/GSE123544_ProcessedDataMatrix_clinical_isolates.xlsx"
df <- readxl::read_excel(snakemake@input[["counts"]]) %>%
  select(locus, all_of(df_meta %>% pull(isolate))) %>%
  separate_wider_delim(
    cols = "locus",
    delim = ",",
    names = c("gene_id", "gene_name"),
    too_few = "align_start"
  ) %>%
  separate_longer_delim(cols = "gene_id", delim = "-")

dds <- DESeqDataSetFromMatrix(
  countData = df %>%
    select(-gene_name) %>%
    column_to_rownames(var = "gene_id"),
  colData = df_meta,
  design = ~mem
)
dds <- DESeq(dds)

df_orthologs <- read_csv(snakemake@input[["orthologs"]]) %>%
  select(`Locus Tag (Query)`, `Locus Tag (Hit)`) %>%
  rename(all_of(c(
    `Locus Tag (Query)` = "new_gene_id",
    `Locus Tag (Hit)` = "gene_id"
  )))


df_lfc <- results(dds) %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene_id") %>%
  dplyr::rename(log2fc = log2FoldChange) %>%
  mutate(regulation = case_when(
    padj < 0.05 & log2fc > 1 ~ "Up",
    padj < 0.05 & log2fc < -1 ~ "Down",
    TRUE ~ "None"
  )) %>%
  left_join(df_orthologs, by = "gene_id")

df_lfc %>% write_csv(snakemake@output[[1]])
