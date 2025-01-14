#!/usr/bin/env python3

import argparse
import pandas as pd
import urllib


def main(output):
    # URL of the UNIPROT REST API
    url = "https://rest.uniprot.org/uniprotkb/stream?"

    # Define columns in the table to be downloaded from UNIPROT
    cols = (
        "accession",
        "reviewed",
        "id",
        "protein_name",
        "gene_names",
        "organism_name",
        "length",
    )

    # Convert query parameters to a string
    query_params = urllib.parse.urlencode(
        {
            "compressed": "false",
            "download": "true",
            "fields": ",".join(cols),
            "format": "tsv",
            "query": "((taxonomy_id:208964))",
        }
    )

    # Download the table and process it before saving it as csv
    df_uniprot = (
        pd.read_csv(url + query_params, sep="\t")
        .rename(columns={"Gene Names": "gene_names"})
        .dropna(subset="gene_names")
        .assign(
            gene_id=lambda x: x.gene_names.str.extract(r"^(?:\S+\s+)*(\S+)$"),
            gene_name=lambda x: x.gene_names.str.extract(r"^(\S+)(?:\s+\S+)*$"),
        )
        .rename(columns=lambda x: x.lower().replace(" ", "_"))
        .rename(columns={"entry": "prot_id", "entry_name": "prot_name"})
    )

    # Additional values to be added to the table
    df_uniprot.loc[df_uniprot.gene_id == "PA3047", "gene_name"] = "dacB"
    df_uniprot.loc[df_uniprot.gene_id == "PA5159", "gene_name"] = "emrA"
    df_uniprot.loc[df_uniprot.gene_id == "PA5160", "gene_name"] = "emrB"

    df_uniprot.to_csv(output, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "get_uniprot_mapping.py", "Download uniprot mapping"
    )
    parser.add_argument("output", help="filepath")
    args = parser.parse_args()

    main(args.output)
