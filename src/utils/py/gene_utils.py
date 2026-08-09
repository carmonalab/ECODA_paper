import pandas as pd
from pathlib import Path


_ENSEMBL105_MAP = None

def _load_ensembl105_map():
    global _ENSEMBL105_MAP
    if _ENSEMBL105_MAP is not None:
        return _ENSEMBL105_MAP
    project_root = Path(__file__).resolve().parents[3]
    path = project_root / "aux" / "EnsemblGenes105_Hsa_GRCh38.p13.txt.gz"
    df = pd.read_csv(path, sep="\t")

    identity = df[["Gene name", "Gene name"]].drop_duplicates()
    identity.columns = ["key", "value"]

    aliases = df[["Gene Synonym", "Gene name"]].copy()
    aliases.columns = ["key", "value"]
    aliases = aliases.dropna(subset=["key"])
    aliases = aliases[aliases["key"] != ""]

    combined = pd.concat([identity, aliases], ignore_index=True)
    combined = combined[~combined["key"].duplicated(keep="first")]

    _ENSEMBL105_MAP = dict(zip(combined["key"], combined["value"]))
    return _ENSEMBL105_MAP


def standardize_gene_symbols(adata):
    gene_map = _load_ensembl105_map()
    adata.var_names = [gene_map.get(g, g) for g in adata.var_names]
