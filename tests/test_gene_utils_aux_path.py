#!/usr/bin/env python3
"""Exercise ECODA auxiliary-map selection and fail-closed overrides."""
from __future__ import annotations

import gzip
import importlib.util
import os
import shutil
import tempfile
from contextlib import contextmanager
from pathlib import Path
from types import ModuleType

ROOT = Path(__file__).resolve().parents[1]
SOURCE_GENE_UTILS = ROOT / "src/utils/py/gene_utils.py"
MAP_NAME = "EnsemblGenes105_Hsa_GRCh38.p13.txt.gz"
_MISSING = object()


@contextmanager
def auxiliary_root(value: Path | None):
    previous = os.environ.get("ECODA_AUX_ROOT", _MISSING)
    if value is None:
        os.environ.pop("ECODA_AUX_ROOT", None)
    else:
        os.environ["ECODA_AUX_ROOT"] = str(value)
    try:
        yield
    finally:
        if previous is _MISSING:
            os.environ.pop("ECODA_AUX_ROOT", None)
        else:
            os.environ["ECODA_AUX_ROOT"] = previous


def load_gene_utils(module_path: Path, module_name: str) -> ModuleType:
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    if spec is None or spec.loader is None:
        raise AssertionError(f"could not load gene_utils from {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_ensembl_map(path: Path, symbol: str, synonym: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write("Gene stable ID\tGene name\tGene Synonym\n")
        handle.write(f"ENSG000001.1\t{symbol}\t{synonym}\n")


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="ecoda-gene-utils-aux-") as raw:
        root = Path(raw)
        module_path = root / "src/utils/py/gene_utils.py"
        module_path.parent.mkdir(parents=True)
        shutil.copy2(SOURCE_GENE_UTILS, module_path)

        fallback_aux = root / "aux"
        override_aux = root / "override-aux"
        write_ensembl_map(fallback_aux / MAP_NAME, "FALLBACK_SYMBOL", "fallback_alias")
        write_ensembl_map(override_aux / MAP_NAME, "OVERRIDE_SYMBOL", "override_alias")

        with auxiliary_root(override_aux):
            override_module = load_gene_utils(module_path, "gene_utils_aux_override")
            override_map = override_module._load_ensembl105_map()
            assert override_map["ENSG000001"] == "OVERRIDE_SYMBOL"
            assert override_map["override_alias"] == "OVERRIDE_SYMBOL"

        with auxiliary_root(None):
            fallback_module = load_gene_utils(module_path, "gene_utils_aux_fallback")
            fallback_map = fallback_module._load_ensembl105_map()
            assert fallback_map["ENSG000001"] == "FALLBACK_SYMBOL"
            assert fallback_map["fallback_alias"] == "FALLBACK_SYMBOL"

        missing_override = root / "missing-aux"
        with auxiliary_root(missing_override):
            invalid_module = load_gene_utils(module_path, "gene_utils_aux_invalid")
            try:
                invalid_module._load_ensembl105_map()
            except FileNotFoundError as exc:
                assert "ECODA_AUX_ROOT" in str(exc)
                assert str(missing_override / MAP_NAME) in str(exc)
            else:
                raise AssertionError("invalid ECODA_AUX_ROOT was accepted via legacy fallback")

    print("gene utils auxiliary path: OK")


if __name__ == "__main__":
    main()
