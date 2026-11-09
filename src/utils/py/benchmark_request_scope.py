"""Exact benchmark requests, normalized family scope, and standalone artifacts.

The public TSV has three fields (whole family) or four (one parameter token).
The original request stays immutable; selection.tsv retains its existing
three-column family schema. Parameter requests never authorize family expansion.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

try:
    from .derivative_artifacts import active_binding, configured_files, context_root
except ImportError:
    from derivative_artifacts import active_binding, configured_files, context_root

BASELINE = ("gloscope", "mofa", "pseudobulk", "composition", "scitd", "mrvi", "scpoli", "pilot", "qot", "pilotgm")
METHODS = (*BASELINE, "prepare_pseudobulk", "trans", "zeroimp")
PB = ("schvg2000", "hvg2000", "hvg500", "hvg2000_bl", "hvg1000", "hvg3000")
CT_HR = ("CT_HR_hvg2000", "CT_HR_hvg500")
AUTHOR_HR = ("ECODA_authors_HR", "ECODA_authors_HR_NULL", "GloProp", "Freq_highres",
    *(f"ECODA_authors_HR_top_varexp{x}" for x in ("0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9")),
    "ECODA_authors_HR_3most_varcts", "ECODA_authors_HR_2least_varcts", "ECODA_authors_HR_3least_varcts",
    *(f"ECODA_authors_HR_{n}_PCA_dims" for n in (2, 3, 5, 10, 15)))
TOKENS = {
    "prepare_pseudobulk": PB,
    "pseudobulk": (*PB, *CT_HR),
    "composition": AUTHOR_HR,
    "gloscope": ("hvg2000_pcadims10", "hvg2000_pcadims30", "hvg2000_pcadims50", "hvg1000_pcadims30", "hvg3000_pcadims30"),
    "mrvi": ("hvg1000", "hvg2000", "hvg3000"),
    "scpoli": ("hvg2000_lowres_dims15", *(f"hvg2000_highres_dims{n}" for n in (2, 3, 5, 10, 15)), "hvg1000_highres_dims15", "hvg3000_highres_dims15"),
    "pilot": ("hvg2000_lowres", "hvg1000_highres", "hvg2000_highres", "hvg3000_highres"),
    "qot": ("hvg2000_lowres", "hvg1000_highres", "hvg2000_highres", "hvg3000_highres"),
    "pilotgm": tuple(f"hvg{n}_highres" for n in (1000, 2000, 3000)),
}


def read_request(path: Path, config: dict) -> dict[tuple[str, str, str], tuple[str, ...] | None]:
    scope = {}
    for number, line in enumerate(path.read_text().splitlines(), 1):
        fields = line.split("\t")
        if len(fields) not in (3, 4) or any(not field for field in fields):
            raise ValueError(f"request row {number} requires three or four nonempty fields")
        dataset, view, method = fields[:3]
        configured_files(config, dataset, view, "stage5")
        if view != "benchmark_analysis":
            raise ValueError("parameter requests are supported only for benchmark_analysis")
        expanded = BASELINE if method == "all" and len(fields) == 3 else (method,)
        for family in expanded:
            if family not in METHODS:
                raise ValueError(f"unknown requested method: {family}")
            key = (dataset, view, family)
            if len(fields) == 3:
                if key in scope:
                    raise ValueError(f"overlapping whole-family request: {key}")
                scope[key] = None
                continue
            token = fields[3]
            if token not in TOKENS.get(family, ()):
                raise ValueError(f"unknown parameter token for {family}: {token}")
            if key in scope and scope[key] is None:
                raise ValueError(f"whole-family/parameter overlap: {key}")
            previous = scope.get(key, ())
            if token in previous:
                raise ValueError(f"duplicate parameter token: {key}/{token}")
            scope[key] = (*previous, token)
    if not scope:
        raise ValueError("empty benchmark request")
    for dataset, view in dict.fromkeys(key[:2] for key in scope):
        selected = scope.get((dataset, view, "pseudobulk"))
        if selected is not None:
            prepared = preparation(scope, dataset, view)
            overlap = set(selected).intersection(PB if prepared is None else prepared)
            if overlap:
                raise ValueError(f"selected pseudobulk/preparation cache overlap: {dataset}/{view}/{sorted(overlap)}")
    return scope


def preparation(scope: dict, dataset: str, view: str) -> tuple[str, ...] | None:
    """Only whole-family dependents or explicit preparation require a prep job.

    Selected pseudobulk variants prepare their own exact caches in the method
    worker, isolating them from an unrelated full-cohort preparation failure.
    """
    variants = set()
    explicit = (dataset, view, "prepare_pseudobulk")
    if explicit in scope:
        if scope[explicit] is None:
            return None
        variants.update(scope[explicit])
    for family, required in (("pseudobulk", PB), ("mofa", ("hvg1000", "hvg2000", "hvg3000")), ("composition", ("hvg2000",))):
        key = (dataset, view, family)
        if key in scope and scope[key] is None:
            variants.update(required)
    if variants == set(PB):
        return None
    return tuple(token for token in PB if token in variants)


def tokens_for(scope: dict, dataset: str, view: str, method: str):
    if method == "prepare_pseudobulk":
        value = preparation(scope, dataset, view)
        return value if value is None or value else "ABSENT"
    return scope.get((dataset, view, method), "ABSENT")


def normalized_rows(scope: dict) -> list[tuple[str, str, str]]:
    rows = []
    for dataset, view in dict.fromkeys(key[:2] for key in scope):
        if tokens_for(scope, dataset, view, "prepare_pseudobulk") != "ABSENT":
            rows.append((dataset, view, "prepare_pseudobulk"))
        rows.extend(key for key in scope if key[:2] == (dataset, view) and key[2] != "prepare_pseudobulk")
    return rows


def selected_files(dataset: str, method: str, tokens: tuple[str, ...]) -> list[str]:
    files = []
    if method == "composition":
        files.append(f"results/{dataset}_metadata.rds")
    for token in tokens:
        if method in ("prepare_pseudobulk", "pseudobulk"):
            if token in PB:
                files.append(f"pseudobulks/{dataset}_pseudobulk_{token}.rds")
            if method == "pseudobulk":
                files.append(f"results/{dataset}_Pseudobulk_{token}.rds")
        elif method == "composition":
            files.append(f"results/{dataset}_{token}.rds")
        elif method == "gloscope":
            files.append(f"results/{dataset}_GloScope_{token}.rds")
        elif method == "mrvi":
            files.append(f"embeddings/{dataset}_{token}_mrvi_dists.feather")
        elif method == "scpoli":
            hvg, resolution, dims = token.split("_")
            files.append(f"embeddings/{dataset}_{hvg}_{resolution}_scpoli_{dims}_embs.feather")
        elif method in ("pilot", "qot", "pilotgm"):
            files.append(f"embeddings/{dataset}_{token}_{method}_dists.feather")
        else:
            raise ValueError(f"method has no standalone parameter artifacts: {method}")
    return files


def active_scope(config_path: Path) -> dict:
    value, _ = active_binding(config_path)
    request = Path(value["run_root"]) / "manifests/requested_selection.tsv"
    if value.get("requested_selection") != str(request):
        raise ValueError("requested parameter selection is not bound to the active run")
    config = json.loads(config_path.read_text())
    scope = read_request(request, config)
    expected = "".join("\t".join(row) + "\n" for row in normalized_rows(scope))
    if Path(value["selection"]).read_text() != expected:
        raise ValueError("normalized selection differs from the exact benchmark request")
    return scope


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("normalize", "tokens", "files", "check", "worker"))
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--request", type=Path)
    parser.add_argument("--dataset")
    parser.add_argument("--view", default="benchmark_analysis")
    parser.add_argument("--method")
    parser.add_argument("--token", default="")
    parser.add_argument("--dispatch", type=Path)
    args = parser.parse_args()
    if args.command == "normalize":
        scope = read_request(args.request, json.loads(args.config.read_text()))
        for row in normalized_rows(scope):
            print("\t".join(row))
        return
    scope = active_scope(args.config)
    if args.command == "check":
        if args.dispatch is not None:
            value, _ = active_binding(args.config)
            if args.dispatch != Path(value["run_root"]) / "manifests/dispatch_selection.tsv":
                raise ValueError("dispatch manifest is not owned by the active run")
            rows = [tuple(line.split("\t")) for line in args.dispatch.read_text().splitlines()]
            expected = normalized_rows(scope)
            if len(rows) != len(set(rows)) or set(rows) != set(expected):
                raise ValueError("dispatch rows differ from the exact benchmark request")
        print("BENCHMARK_REQUEST_VALIDATED")
        return
    tokens = tokens_for(scope, args.dataset, args.view, args.method)
    if args.command == "worker":
        if tokens == "ABSENT":
            raise ValueError("worker method is outside the exact benchmark request")
        if tokens is not None:
            expected = ",".join(tokens)
            if args.method in ("prepare_pseudobulk", "pseudobulk", "composition"):
                valid = args.token == expected
            else:
                valid = args.token in tokens
            if not valid:
                raise ValueError("worker parameter scope differs from the exact benchmark request")
        print("BENCHMARK_WORKER_SCOPE_VALIDATED")
        return
    if args.command == "tokens":
        print("FULL" if tokens is None else tokens if tokens == "ABSENT" else ",".join(tokens))
    else:
        if tokens is None or tokens == "ABSENT":
            raise ValueError("standalone artifact expansion requires selected parameter tokens")
        value, _ = active_binding(args.config)
        root = context_root(value, "benchmark")
        if os.environ.get("ANALYSIS_ROOT") != str(root):
            raise ValueError("standalone artifact root differs from the active binding")
        for relative in selected_files(args.dataset, args.method, tokens):
            print(root / relative)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, KeyError, TypeError, OSError) as exc:
        raise SystemExit(f"ERROR: {exc}")
