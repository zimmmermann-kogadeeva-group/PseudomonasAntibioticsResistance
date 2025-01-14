#!/usr/bin/env python3

from datetime import datetime
from functools import wraps
import hashlib
from io import StringIO
import inspect
import json
import pandas as pd
import re
import requests
from time import sleep


def correct_term_col(value):
    if type(value) is dict:
        if "id" in value:
            return f"{value['label']} ({value['id']})"
        else:
            return value["label"]
    else:
        return value


def get_ora(
    genes,
    organism,
    annotDataSet="GO:0008150",
    enrichmentTestType="FISHER",
    correction="FDR",
    cache=None,
):
    if not isinstance(genes, str):
        genes = ",".join(genes)

    key = genes + str(organism) + annotDataSet + enrichmentTestType + correction
    if cache is not None and cache.get(key) is not None:
        return cache.get(key)

    # The default value for annotDataSet parameter means biological process GO terms
    query_params = {
        "geneInputList": genes,
        "organism": organism,
        "annotDataSet": annotDataSet,
        "enrichmentTestType": enrichmentTestType,
        "correction": correction,
    }

    # Send a GET request to GO API
    response = requests.post(
        "https://pantherdb.org/services/oai/pantherdb/enrich/overrep",
        headers={"accept": "application/json"},
        data=query_params,
    )

    df_res = pd.DataFrame()
    if response and "results" in response.json():
        # Convert to a dataframe, then correct term column and append
        # additional columns
        df_res = pd.DataFrame(response.json()["results"]["result"]).assign(
            term=lambda x: x.term.apply(correct_term_col)
        )

    if cache is not None:
        cache.set(key, df_res)

    sleep(2)

    return df_res


# TODO: check that get_gsea still works
def get_gsea(
    genes,
    organism,
    annotDataSet="GO:0008150",
    correction="FDR",
    cache=None,
):
    key = genes + str(organism) + annotDataSet + correction
    if cache is not None and cache.get(key) is not None:
        return cache.get(key)
    # genes param should be a dataframe with ranks
    # The default value for annotDataSet parameter means biological process GO terms

    # Set up the dataframe for the request
    tsv_data = StringIO()
    genes.to_csv(tsv_data, sep="\t", index=False)
    tsv_data.seek(0)  # Reset the stream position to the beginning
    files = {"geneExp": ("data.tsv", tsv_data, "text/tab-separated-values")}

    response = requests.post(
        "https://pantherdb.org/services/oai/pantherdb/enrich/statenrich",
        headers={"accept": "application/json"},
        data={
            "organism": organism,
            "annotDataSet": annotDataSet,
            "correction": correction,
        },
        files=files,
    )

    df_res = pd.DataFrame()
    if response and len(response.content) and "results" in response.json():
        df_res = pd.DataFrame(response.json()["results"]["result"]).assign(
            term=lambda x: x.term.apply(correct_term_col)
        )

    if cache is not None:
        cache.set(key, df_res)

    sleep(2)

    return df_res
