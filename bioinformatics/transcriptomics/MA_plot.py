from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from adjustText import adjust_text
from matplotlib.axes import Axes
from matplotlib.lines import Line2D
from pandas import DataFrame
import requests


def get_tf_list() -> list:
    taxon_id = "9606"  # Human
    """
    ════════════════════════════════════════════════
    example go terms:
    - 0009986 Cell surface
    - 0048870 Cell motility
    - 0003700 transcription factors
    - 0030217 T cell
    - 0016020 cell membrane
    - 0005576 ECM
    - 0003700 transcription factors
    - 0005737 cytoplasm
    - 0005840 ribosome
    - 0005739 mitochondria
    - 0031047 negative regulation of gene expression
    ════════════════════════════════════════════════
    """
    # specify the go id you want here
    go_terms = ["0003700"]
    descriptions = ["TF"]
    df_out = None
    for i, gt in enumerate(go_terms):
        try:
            url_genes_by_go = f"https://api.geneontology.org/api/bioentity/function/GO%3A{gt}/genes?taxon=NCBITaxon%3A{taxon_id}&relationship_type=involved_in&rows=999999"
            response = requests.get(url_genes_by_go).json()["associations"]
            print(f"Got {len(response)} entries.")
            print(response[0])
            df_out_next = pd.json_normalize(response)
            df_out_next["go_name"] = descriptions[i]
            if df_out is None:
                df_out = df_out_next
            else:
                df_out = pd.concat([df_out, df_out_next], axis=0)
        except requests.exceptions.HTTPError as error:
            print(error)
    return list(np.unique(df_out["subject.label"]))


def volcano_plot_differential_expression(
    df: DataFrame,
    LFC_key: str,
    xval: str,
    annot_key: str,
    min_lfc: float = 1.0,
    min_exp: float = 0.0,
    y_annot_up: str = "",
    y_annot_down: str = "",
    y_annot_xpos: float = 0.75,
    y_annot_ypos: float = 0.75,
    y_max: float = np.inf,
    legend_loc: str = "upper right",
    legend_fontsize: int = 8,
    gene_domain: list[str] = None,
    ax: Optional[Axes] = None,
    **kwargs,
) -> Axes:
    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (5, 5)))
    np.random.seed(42)
    df = df.iloc[np.random.permutation(df.shape[0]), :]
    df["scatter_color"] = "#000000"
    df_sig = df.query(
        f"{xval} > {min_exp} & ({LFC_key} >= {min_lfc} | {LFC_key} <= -{min_lfc})"
    ).reset_index()

    ymax = max(np.max(df[LFC_key]), np.max(-df[LFC_key]))
    ymax = max(ymax, y_max)
    ymin = -ymax * 1.1
    y_clip = df[LFC_key].__array__()
    y_clip[y_clip > ymax] = ymax
    y_clip[y_clip < -ymax] = -ymax
    y_sig_clip = df_sig[LFC_key].__array__()
    y_sig_clip[y_sig_clip > ymax] = ymax
    y_sig_clip[y_sig_clip < -ymax] = -ymax

    ax.scatter(
        df[xval],
        y_clip,
        c=df["scatter_color"],
        s=0.5,
        alpha=0.5,
    )
    ax.scatter(
        df_sig[xval],
        y_sig_clip,
        color="red",
        s=4,
        alpha=0.5,
    )
    texts = []

    for i, row in df_sig.iterrows():
        if gene_domain is not None:
            if row[annot_key] not in gene_domain:
                continue
        texts.append(
            ax.text(
                x=row[xval],
                y=y_sig_clip[i],
                s=row[annot_key],
                fontsize=8,
                alpha=1.0,
                ha="left",
                va="bottom",
                bbox=dict(
                    boxstyle="round,pad=0.1", facecolor="white", edgecolor="black"
                ),
            )
        )

    xlims = ax.get_xlim()
    xmin = -0.5
    xmax = xlims[1] * 1.1

    ax.text(
        x=xmax * y_annot_xpos,
        y=ymax * y_annot_ypos,
        s=y_annot_up,  # assuming gene names are in 'names' column
        fontsize=16,
        alpha=1.0,
        ha="left",  # horizontal alignment
        va="bottom",  # vertical alignment
    )
    ax.text(
        x=xmax * y_annot_xpos,
        y=-ymax * y_annot_ypos,
        s=y_annot_down,  # assuming gene names are in 'names' column
        fontsize=16,
        alpha=1.0,
        ha="left",  # horizontal alignment
        va="top",  # vertical alignment
    )
    ax.annotate(
        "",
        xy=(xmax * y_annot_xpos, ymax * y_annot_ypos * 0.95),
        xytext=(xmax * y_annot_xpos, ymax * y_annot_ypos * 0.6),
        arrowprops=dict(arrowstyle="->", color="black", lw=2),
    )
    ax.annotate(
        "",
        xy=(xmax * y_annot_xpos, -ymax * y_annot_ypos * 0.95),
        xytext=(xmax * y_annot_xpos, -ymax * y_annot_ypos * 0.6),
        arrowprops=dict(arrowstyle="->", color="black", lw=2),
    )
    legend_elements = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor="black",
            markersize=6,
            label="not differentially expressed",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor="red",
            markersize=6,
            label="unique markers",
        ),
    ]
    ax.legend(
        handles=legend_elements,
        loc=legend_loc,
        frameon=True,
        fancybox=True,
        shadow=True,
        fontsize=legend_fontsize,
    )
    ax.set_xticks(np.round(np.linspace(0, xmax, 5)))
    yticks = np.round(np.linspace(0, ymax, 4), 1)
    yticks = np.sort(np.unique(np.concatenate([-yticks, yticks])))
    yticks_lab = yticks.astype(str)
    yticks_lab[0] = f"<= -{np.round(ymax, 1)}"
    yticks_lab[-1] = f">= {np.round(ymax, 1)}"
    ax.set_yticks(yticks, labels=yticks_lab)
    ax.set_ylabel(r"Log$_{2}$ Fold Change")
    ax.set_xlabel("Mean Expression")

    ax.set_ylim([ymin, ymax * 1.1])
    ax.set_xlim([xmin, xmax])

    adjust_text(
        texts,
        arrowprops=dict(
            arrowstyle="-",
            color="black",
            lw=1,
            alpha=0.5,
            shrinkA=10,
            shrinkB=0,
        ),
        expand=(1.25, 1.25),
        avoid_text=True,
        avoid_point=True,
        prevent_crossings=True,
        min_arrow_len=1,
        max_move=[15, 15],
    )
    return ax


def main():
    TSV_PATH = "88NFMB-expression-matrix.tsv"
    OUTPUT_PATH_1 = "MA_plot_BMP4_4h_pseudocount_1.png"
    OUTPUT_PATH_2 = "MA_plot_BMP4_4h_pseudocount_p01.png"
    CTRL_COL = "88NFMB_1_cpm"
    TREAT_COL = "88NFMB_2_cpm"
    TF_LIST = get_tf_list()
    # here I usually use 1 to shrink high fc of lowly expressed genes 
    # (e.g. (0.5 + 1) / (0.1 + 1) = 1.36 whereas 0.5 / 0.1 = 5). 
    # Also try to be consistent with x, which is natural log(1 + mean cpm)
    PSEUDOCOUNT = 1
    df = pd.read_csv(TSV_PATH, sep="\t")
    df["log2fc"] = np.log2((df[TREAT_COL] + PSEUDOCOUNT) / (df[CTRL_COL] + PSEUDOCOUNT))
    df["mean_cpm"] = (df[CTRL_COL] + df[TREAT_COL]) / 2
    df["x"] = np.log(df["mean_cpm"] + 1)
    # add IDs
    TF_LIST += list(df["gene_name"][df["gene_name"].str.contains("^ID", na=False)])
    ax = volcano_plot_differential_expression(
        df=df,
        xval="x",
        LFC_key="log2fc",
        annot_key="gene_name",
        min_lfc=1.3,  # these are not for filtering, only for highlighting differentially expressed genes
        min_exp=0.5,
        gene_domain=TF_LIST,
        legend_loc="lower left",
        legend_fontsize=10,
        y_max=6.2,
        y_annot_up="BMP4",
        y_annot_down="Control",
        y_annot_xpos=0.67,
        y_annot_ypos=0.68,
    )
    ax.set_xlim(ax.get_xlim()[0], 10)
    ax.set_xticks([0, 2, 3, 4, 6, 8, 10])
    ax.set_xlabel("log(mean_CPM + 1)", fontsize=14, fontweight="bold")
    ax.set_ylabel(ax.get_ylabel(), fontsize=14, fontweight="bold")
    ax.set_title("RNAseq (Plasmidsaurus)", fontsize=14, fontweight="bold")
    plt.savefig(OUTPUT_PATH_1, dpi=300, bbox_inches="tight")

    # Here repeat Minja's result with 0.01 pseudocount
    PSEUDOCOUNT = 0.01
    df["log2fc"] = np.log2((df[TREAT_COL] + PSEUDOCOUNT) / (df[CTRL_COL] + PSEUDOCOUNT))
    # add IDs
    TF_LIST += list(df["gene_name"][df["gene_name"].str.contains("^ID", na=False)])
    ax = volcano_plot_differential_expression(
        df=df,
        xval="x",
        LFC_key="log2fc",
        annot_key="gene_name",
        min_lfc=1.3,
        min_exp=0.5,
        gene_domain=TF_LIST,
        legend_loc="lower left",
        legend_fontsize=10,
        y_max=6.2,
        y_annot_up="BMP4",
        y_annot_down="Control",
        y_annot_xpos=0.67,
        y_annot_ypos=0.68,
    )
    ax.set_xlim(ax.get_xlim()[0], 10)
    ax.set_xticks([0, 2, 3, 4, 6, 8, 10])
    ax.set_xlabel("log(mean_CPM + 1)", fontsize=14, fontweight="bold")
    ax.set_ylabel(ax.get_ylabel(), fontsize=14, fontweight="bold")
    ax.set_title("RNAseq (Plasmidsaurus)", fontsize=14, fontweight="bold")
    plt.savefig(OUTPUT_PATH_2, dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main()
