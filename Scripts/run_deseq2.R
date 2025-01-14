#!/usr/bin/env Rscript

library(tidyverse)
library(tximport)
library(DESeq2)

get_contrast <- function(strain, comp) {
    return(c(
        "group", 
        paste0(strain, "_", str_sub(comp, 1, 1)), 
        paste0(strain, "_", str_sub(comp, 4, 4))
    ))
}

main <- function(quant_files, comparisons, annot_path, output_dir) {
    samples <- basename(dirname(quant_files))

    df_samples <- data.frame(
        row.names=samples, 
        sample=samples, 
        files=quant_files
    ) %>%
    separate_wider_delim(
        cols="sample", 
        delim="_", 
        names=c("strain", "condition", "replicate"), 
        cols_remove=F
    )

    selected_strains <- df_samples %>% pull(strain) %>% unique()

    # Read in the annotation file
    tx2gene <- read_tsv(annot_path, show_col_types=F) %>% 
    select(transcript_stable_id, gene_stable_id)

    # Create metadata table
    meta <- df_samples %>%
    mutate(group=factor(paste(strain, condition, sep="_"))) %>%
    select(sample, condition, group) %>% 
    column_to_rownames(var="sample")

    # Import all the salmon outputs
    txi <- tximport(
        df_samples[["files"]], 
        type="salmon", 
        tx2gene=tx2gene, 
        countsFromAbundance="lengthScaledTPM"
    )

    # Run DESeq2 and save to rds
    dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ group)
    dds <- DESeq(dds)
    saveRDS(dds, paste0(output_dir, "/dds.rds"))

    # Get raw counts and save to csv
    dds %>%
    counts() %>% 
    as.data.frame() %>% 
    rownames_to_column(var="gene_id") %>% 
    write_csv(file=paste0(output_dir, "/counts.csv"))
    
    # Get normalised counts and save to csv
    dds %>%
    counts(normalized=TRUE) %>% 
    as.data.frame() %>% 
    rownames_to_column(var="gene_id") %>% 
    write_csv(file=paste0(output_dir, "/counts_normalised.csv"))

    # Get variance stabilising transformation of counts
    dds %>% 
    rlog() %>% 
    assay() %>% 
    as.data.frame() %>%
    rownames_to_column(var="gene_id") %>%
    write_csv(file=paste0(output_dir, "/rlog.csv"))

    # Get variance stabilising transformation of counts
    dds %>% 
    vst() %>% 
    assay() %>% 
    as.data.frame() %>%
    rownames_to_column(var="gene_id") %>%
    write_csv(file=paste0(output_dir, "/vst.csv"))

    expand_grid(
        strain=selected_strains, 
        comparison=comparisons
    ) %>% 
    filter(!(strain == "090.3" & str_detect(comparison, "M"))) %>%
    rowwise() %>%
    mutate(both = list(c_across(everything()))) %>%
    pull(both) %>%
    map(
        ~{
            results(dds, contrast=get_contrast(.x[[1]], .x[[2]])) %>%
            as.data.frame() %>%
            rownames_to_column(var="gene_id") %>%
            mutate(strain=.x[[1]], comparison=.x[[2]])
        }
    ) %>%
    bind_rows() %>%
    dplyr::rename(log2fc=log2FoldChange) %>%
    write_csv(paste0(output_dir, "/dds_results.csv"))

}

main(
     snakemake@input[["quant_files"]], 
     snakemake@config[["comparisons"]],
     snakemake@input[["annot"]],
     snakemake@params[["output_dir"]]
)

