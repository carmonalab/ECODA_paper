import os
from pathlib import Path

import pandas as pd

_ENSEMBL105_MAP = None

def _load_ensembl105_map():
    global _ENSEMBL105_MAP
    if _ENSEMBL105_MAP is not None:
        return _ENSEMBL105_MAP
    project_root = Path(__file__).resolve().parents[3]
    filename = "EnsemblGenes105_Hsa_GRCh38.p13.txt.gz"
    if "ECODA_AUX_ROOT" in os.environ:
        aux_root = os.environ["ECODA_AUX_ROOT"]
        path = Path(aux_root) / filename
        if not aux_root or not path.is_file() or not os.access(path, os.R_OK):
            raise FileNotFoundError(
                f"ECODA_AUX_ROOT Ensembl map is missing or unreadable: {path}"
            )
    else:
        path = project_root / "aux" / filename
    df = pd.read_csv(path, sep="\t")

    stable_ids = df[["Gene stable ID", "Gene name"]].copy()
    stable_ids.columns = ["key", "value"]
    stable_ids["key"] = (
        stable_ids["key"].astype("string").str.replace(r"\.[0-9]+$", "", regex=True)
    )
    stable_ids = stable_ids.dropna(subset=["key", "value"])
    stable_ids = stable_ids[
        (stable_ids["key"] != "") & (stable_ids["value"] != "")
    ]

    identity = df[["Gene name", "Gene name"]].drop_duplicates()
    identity.columns = ["key", "value"]

    aliases = df[["Gene Synonym", "Gene name"]].copy()
    aliases.columns = ["key", "value"]
    aliases = aliases.dropna(subset=["key"])
    aliases = aliases[aliases["key"] != ""]

    combined = pd.concat([stable_ids, identity, aliases], ignore_index=True)
    combined = combined[~combined["key"].duplicated(keep="first")]

    _ENSEMBL105_MAP = dict(zip(combined["key"], combined["value"]))
    return _ENSEMBL105_MAP


def standardize_gene_symbols(adata):
    gene_map = _load_ensembl105_map()
    standardized = []
    for gene in adata.var_names:
        text = str(gene)
        base = text.split(".", 1)[0]
        standardized.append(gene_map.get(text, gene_map.get(base, text)))
    adata.var_names = standardized
