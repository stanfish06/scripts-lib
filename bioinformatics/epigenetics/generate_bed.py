#!/usr/bin/env python
import pandas as pd

genes_to_exclude = ["ENSG00000308881"]
path_to_csv = "./hg19.csv"
output_name = "./hg19.bed"


def main():
    ucsc_df = pd.read_csv(path_to_csv)
    ucsc_df = ucsc_df[[g not in genes_to_exclude for g in ucsc_df["geneName"]]]
    bed12_fields = [
        "chrom",
        "chromStart",
        "chromEnd",
        "name",
        "score",
        "strand",
        "thickStart",
        "thickEnd",
        "itemRgb",
        "blockCount",
        "blockSizes",
        "blockStarts",
    ]
    ucsc_df.loc[:, "chrom"] = ucsc_df.iloc[:, 0]
    ucsc_df.loc[:, "score"] = 500
    ucsc_df.loc[:, "itemRgb"] = "0,0,255"
    ucsc_df.loc[:, "blockStarts"] = ucsc_df.loc[:, "chromStarts"]
    column_order = [
        "chrom",
        "chromStart",
        "chromEnd",
        "geneName",
        "score",
        "strand",
        "thickStart",
        "thickEnd",  # thickStart is TSS
        "itemRgb",
        "blockCount",
        "blockSizes",
        "blockStarts",
    ]
    ucsc_df = ucsc_df[column_order]
    ucsc_df.columns = bed12_fields
    # sort by first column alphabetically
    # sort second and third columns numerically
    ucsc_df = ucsc_df.sort_values(
        by=[bed12_fields[0], bed12_fields[1], bed12_fields[2]],
        ascending=[True, True, True],
    )
    ucsc_df.to_csv(
        output_name,
        sep="\t",
        header=False,
        index=False,
    )


if __name__ == "__main__":
    main()
