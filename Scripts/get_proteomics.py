#!/usr/bin/env python3

from functools import reduce
import numpy as np
import pandas as pd
import re


def separate(data, col, new_cols, sep=",", **kwargs):
    new_data = data.copy()
    new_data[new_cols] = data[col].str.split(sep, **kwargs)
    return new_data


def rename_lfc_cols(filename):
    # Rename log2FC, P.Value and adj.P.Val columns to avoid conflicts when
    # merging tables. `zfill` pads with zeroes to make strain label consistant
    comparison = re.sub(
        r".*_report_([0-9]*)([MCP])([0-9])_vs_.*([MCP]).*_fc1_.*",
        "\\1.\\3_\\2vs\\4",
        filename,
    ).zfill(10)
    return {x: f"{x}_{comparison}" for x in ("log2FC", "P.Value", "adj.P.Val")}


def merge_lfc(x, y):
    cols = "^(protein_name|log2FC|P.Value|adj.P.Val)"
    return pd.merge(
        x.filter(regex=cols), y.filter(regex=cols), on="protein_name", how="outer"
    )


def rename_raw_cols(col):
    if col.endswith(".raw"):
        return re.sub(
            ".*_.*_([0-9]*)([MPC])([0-9])_([A-D]).*\.raw",
            "\\1.\\3_\\2_\\4.raw",
            col,
        ).zfill(13)
    else:
        return col


def merge_raw(x, y):
    x = x.filter(regex="\.raw$|^protein_name").copy()
    y = y.filter(regex="\.raw$|^protein_name").copy()

    df_merged = pd.merge(x, y, on="protein_name", how="outer")
    for col in set(x.columns).intersection(y.columns) - {"protein_name"}:
        df_merged[col] = df_merged[col + "_x"].fillna(df_merged[col + "_y"])
        df_merged.drop(columns=[col + "_x", col + "_y"], inplace=True)
    return df_merged


def quantile_norm(data):
    rank_median = (
        pd.DataFrame(
            {
                name: col.sort_values(na_position="first").values
                for name, col in data.items()
            }
        )
        .median(axis=1)
        .tolist()
    )
    return data.apply(
        lambda x: [
            np.nanpercentile(rank_median, i * 100) if ~np.isnan(i) else np.nan
            for i in x.rank(pct=True, method="max")
        ]
    )


def main(input_files, output_dir):
    all_files = [
        pd.read_csv(filename, sep="\t", index_col=0)
        .rename(columns=rename_lfc_cols(filename))
        .rename(columns=rename_raw_cols)
        .rename(columns={"ProteinName": "protein_name"})
        for filename in input_files
    ]

    df_raw = (
        reduce(merge_raw, all_files)
        .rename(columns=lambda x: x.removesuffix(".raw"))
        .assign(protein_name=lambda x: x.protein_name.str.split(";"))
        .explode("protein_name")
        .set_index("protein_name")
        .filter(regex="^(?!090.3_M).*")
    )
    df_raw.to_csv(f"{output_dir}/raw.csv")

    df_norm = df_raw.pipe(np.log10).pipe(quantile_norm)
    df_norm.to_csv(f"{output_dir}/normalised.csv")

    df_lfc = (
        reduce(merge_lfc, all_files)
        .pipe(lambda x: x.assign(**(-x.filter(regex="log2FC.*MvsC"))))
        .rename(columns=lambda x: x.replace("MvsC", "CvsM"))
        .pipe(
            pd.wide_to_long,
            stubnames=("log2FC", "P.Value", "adj.P.Val"),
            i="protein_name",
            j="comparison",
            sep="_",
            suffix=".*",
        )
        .reset_index()
        .rename(columns={"log2FC": "log2fc", "P.Value": "pvalue", "adj.P.Val": "padj"})
        .assign(protein_name=lambda x: x.protein_name.str.split(";"))
        .explode("protein_name")
        .pipe(
            separate,
            "protein_name",
            ["rev", "prot_id", "protein_name"],
            sep="|",
            expand=True,
        )
        .pipe(separate, "comparison", ["strain", "comparison"], sep="_", expand=True)
        .query("~(strain == '090.3' and comparison.str.contains('M'))")
    )
    df_lfc.to_csv(f"{output_dir}/lfc.csv", index=False)


if __name__ == "__main__":
    main(snakemake.input, snakemake.params["output_dir"])
