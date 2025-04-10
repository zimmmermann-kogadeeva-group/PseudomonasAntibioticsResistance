#!/usr/bin/env Rscript

# Data from: https://www.ebi.ac.uk/ena/browser/view/PRJNA1066021

library(tidyverse)
library(DESeq2)
library(tximport)


# quant_files <- list.files(
#   here("Output/PublicDatasets/PRJNA1066021/Quants"),
#   pattern = "quant.sf",
#   recursive = TRUE,
#   full.names = TRUE
# )
quant_files <- snakemake@input[["quants"]]

samples <- quant_files %>%
  dirname() %>%
  basename()

df_samples <- data.frame(
  row.names = samples,
  sample = samples,
  files = quant_files
) %>%
  filter(str_detect(sample, "MEM|control")) %>%
  separate_wider_delim(
    cols = "sample",
    delim = "_",
    names = c("species", "cond", "repl"),
    cols_remove = FALSE
  )

tx2gene <- read_tsv(snakemake@input[["annot"]]) %>%
  select(transcript_stable_id, gene_stable_id)

# Create metadata table
meta <- df_samples %>%
  mutate(cond = factor(cond)) %>%
  select(sample, cond) %>%
  column_to_rownames(var = "sample")

# Import all the salmon outputs
txi <- tximport(
  df_samples %>% pull(files),
  type = "salmon",
  tx2gene = tx2gene,
  countsFromAbundance = "lengthScaledTPM"
)

# Run DESeq2 and save to rds
dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~cond)
dds <- DESeq(dds)

df_ids <- read_csv(snakemake@input[["id_map"]]) %>%
  select(gene_id, gene_name)

df_lfc <- results(dds) %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene_id") %>%
  dplyr::rename(log2fc = log2FoldChange) %>%
  mutate(regulation = case_when(
    padj < 0.05 & log2fc > 1 ~ "Up",
    padj < 0.05 & log2fc < -1 ~ "Down",
    TRUE ~ "None"
  )) %>%
  left_join(df_ids, by = "gene_id")

df_lfc %>% write_csv(snakemake@output[[1]])
