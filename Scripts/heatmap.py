#!/usr/bin/env python3

from itertools import product

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def get_mut_idxs(mut, gene_names, comp_strain_omics, add_cols=None):
    if add_cols is None:
        add_cols = []

    return (
        pd.DataFrame(
            [
                (i + 0.5, j + 1, gene, comp, strain, omics)
                for (i, gene), (j, (comp, strain, omics)) in product(
                    enumerate(gene_names), enumerate(comp_strain_omics)
                )
                if omics == "t"
            ],
            columns=[
                "pos_y",
                "pos_x",
                "gene_name",
                "condition",
                "strain",
                "omics",
            ],
        )
        .merge(mut, how="inner", on=["gene_name", "condition", "strain"])
        .get(["pos_y", "pos_x", "gene_name", "condition", "strain", "type", *add_cols])
        .drop_duplicates()
    )


def get_label_pos(index):
    return (
        index.to_series()
        .groupby(level=0)
        .count()
        .reindex(index.unique())
        .reset_index(name="length")
        .assign(
            start=lambda x: x.length.shift(1, fill_value=0).cumsum(),
            middle=lambda x: x.start + x.length / 2,
            end=lambda x: x.length.cumsum(),
        )
    )


# MIC array barplot
def plot_mic(data, ax):
    num_rows, num_cols = data.shape
    mic_arr = data.to_numpy().reshape(-1)
    x0 = np.concatenate([np.arange(num_cols)] * num_rows)
    y = np.repeat(np.arange(num_rows), num_cols)

    ax.barh(y, mic_arr / 256, left=x0, height=1, color="green")
    ax.tick_params(bottom=False, labelbottom=False)
    ax.set(
        yticks=np.arange(num_rows),
        yticklabels=data.index.tolist(),
        xlim=(0, num_cols),
    )

    for x_pos, y_pos, val in zip(x0, y, mic_arr):
        ax.text(
            x_pos + 0.5,
            y_pos,
            f"{val:.0f}",
            ha="center",
            va="center",
        )


def plot_mut(data, ax, legend_kwargs=None, **kwargs):
    legend_opts = dict(loc="upper left", bbox_to_anchor=(1.02, 1.07))
    if legend_kwargs is not None:
        legend_opts.update(legend_kwargs)

    if not data.empty:
        main_opts = dict(
            x="pos_x",
            y="pos_y",
            style="Mutation type",
            color="black",
            ax=ax,
            style_order=("snp", "del"),
        )
        main_opts.update(kwargs)

        # Mutations
        sns.scatterplot(data, **main_opts)
        sns.move_legend(ax, **legend_opts)
        ax.set(xlabel=None, ylabel=None)


def plot_heatmap(
    fig,
    data_lfc,
    data_mic=None,
    data_mut=None,
    mut_kwargs=None,
    mut_legend_kwargs=None,
    x_label_kwargs=None,
    y_label_kwargs=None,
    cbar_pad=0.04,
    **kwargs,
):
    x_label_opts = dict(
        line_pos=(1.15, 1.2),
        text_pos=-5.5,
        pad=20,
        rotation=0,
        lw=1.5,
    )
    if x_label_kwargs is not None:
        x_label_opts.update(x_label_kwargs)

    y_label_opts = dict(
        line_pos=(0.15, 0.8),
        text_pos=-12,
        pad=5,
        rotation=0,
        lw=1.5,
    )
    if y_label_kwargs is not None:
        y_label_opts.update(y_label_kwargs)

    heatmap_opts = dict(
        cmap="coolwarm",
        center=0,
        vmin=-2,
        vmax=2,
    )
    heatmap_opts.update(kwargs)

    axs = fig.get_axes()
    assert len(axs) == 2

    if "cbar_ax" not in heatmap_opts:
        cax = axs[1].inset_axes([1.0 + cbar_pad, 0.0, 0.05, 1.0])
        heatmap_opts["cbar_ax"] = cax

    sns.heatmap(data_lfc, **heatmap_opts, ax=axs[1])
    axs[1].tick_params(
        "x",
        top=False,
        labeltop=True,
        bottom=False,
        labelbottom=False,
        pad=x_label_opts["pad"],
        rotation=x_label_opts["rotation"],
    )
    axs[1].set(
        xlabel="",
        ylabel="",
        yticks=np.arange(0, data_lfc.shape[0]) + 0.5,
        yticklabels=data_lfc.index.get_level_values(-1),
        xticks=np.arange(0, data_lfc.shape[1], 2) + 1,
        xticklabels=data_lfc.columns.get_level_values(1)[::2],
    )

    # Vertical lines and lables
    for x_val in range(0, data_lfc.shape[1] + 2, 2):
        line = axs[1].axvline(
            x_val,
            ymax=x_label_opts["line_pos"][0],
            lw=x_label_opts["lw"],
            c="k",
            clip_on=False,
        )

    for idx, col in get_label_pos(data_lfc.columns.get_level_values(0)).iterrows():
        axs[1].text(
            col.middle,
            x_label_opts["text_pos"],
            col.comparison,
            ha="center",
            va="center",
        )
        axs[1].axvline(
            col.start,
            x_label_opts["line_pos"][0],
            x_label_opts["line_pos"][1],
            lw=x_label_opts["lw"],
            c="k",
            clip_on=False,
        )
    axs[1].axvline(
        data_lfc.shape[1],
        x_label_opts["line_pos"][0],
        x_label_opts["line_pos"][1],
        lw=x_label_opts["lw"],
        c="k",
        clip_on=False,
    )

    # Horizontal lines and labels
    for y_val in range(0, data_lfc.shape[0] + 1, 1):
        axs[1].axhline(
            y_val,
            xmin=-y_label_opts["line_pos"][0],
            lw=y_label_opts["lw"],
            c="k",
            clip_on=False,
        )

    if data_lfc.index.nlevels > 1:
        for idx, row in get_label_pos(data_lfc.index.get_level_values(-2)).iterrows():
            axs[1].text(
                y_label_opts["text_pos"],
                row.middle,
                row.group,
                ha="center",
                va="center_baseline",
            )
            axs[1].axhline(
                row.start,
                -y_label_opts["line_pos"][1],
                -y_label_opts["line_pos"][0],
                lw=y_label_opts["lw"],
                c="k",
                clip_on=False,
            )
        axs[1].axhline(
            data_lfc.shape[0],
            -y_label_opts["line_pos"][1],
            -y_label_opts["line_pos"][0],
            lw=y_label_opts["lw"],
            c="k",
            clip_on=False,
        )

    plot_mic(data_mic, axs[0])
    mut_kwargs = mut_kwargs or {}
    plot_mut(data_mut, axs[1], legend_kwargs=mut_legend_kwargs, **mut_kwargs)
