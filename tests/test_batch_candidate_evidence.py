#!/usr/bin/env python3
"""Tiny synthetic evidence run covering integrity and applicability gates."""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
DATASETS = json.loads((ROOT / "datasets.json").read_text())
SPECS_NS = {}
exec((ROOT / "notebooks/dataset_onboarding/dataset_specs.py").read_text(), SPECS_NS)
ORDER = SPECS_NS["BATCH_EFFECT_DATASET_ORDER"]
CANDIDATES = SPECS_NS["BATCH_EFFECT_SPECS"]
LABELS = {
    "Alzheimer": "Cognitive status",
    "Breast_cancer": "disease",
    "Covid19_PBMC": "CoVID-19 severity",
    "Kidney_KPMP": "condition.l1",
    "Myocardial_infarction": "patient_group",
    "Diabetes": "disease",
    "Lupus_PBMC": "Status",
    "Lung": "disease",
    "Parkinson": "disease",
    "Joanito": "sample.origin",
    "Stephenson": "Status",
    "CombinedPBMC": "cond",
}
COMPOSITION = [
    "ECODA_authors_HR",
    "ECODA_HiTME_HR_layer2",
    "ECODA_scATOMIC_HR",
    "ECODA_seuratres_2",
    "ECODA_authors_HR_NULL",
]


def sidecar(path: Path) -> None:
    path.with_name(f"{path.name}.md5").write_text(
        f"MD5={hashlib.md5(path.read_bytes()).hexdigest()}\n"
        f"SIZE={path.stat().st_size}\nPATH={path}\n"
    )


def make_h5ads(input_root: Path) -> None:
    for ds in ORDER:
        output_name = DATASETS[ds]["views"]["batch_effect_uncorrected"]["output_file_name"]
        path = input_root / ds / "output" / output_name
        path.parent.mkdir(parents=True)
        samples = [f"s{i}" for i in range(1, 7)]
        labels = ["A", "A", "B", "B", "C", "C"]
        obs: dict[str, object] = {"Sample": samples, LABELS[ds]: labels}
        for candidate in CANDIDATES[ds]:
            values: list[object]
            if ds == "Alzheimer" and candidate == "assay":
                continue  # missing candidate row
            if ds == "Breast_cancer" and candidate == "assay":
                values = ["A", "A", "B", "B", "C", None]  # incomplete
            elif ds == "Covid19_PBMC" and candidate == "City":
                values = ["A"] * 6  # constant
            elif ds == "Kidney_KPMP" and candidate == "experiment":
                values = [f"u{i}" for i in range(6)]  # sample-unique
            elif ds == "Myocardial_infarction" and candidate == "batch":
                values = labels[:]  # perfectly confounded
            else:
                values = ["A", "A", "B", "B", "A", "B"]
            obs[candidate] = values
        matrix = sparse.csr_matrix(np.ones((6, 2), dtype="float32"))
        data = ad.AnnData(
            X=matrix,
            obs=pd.DataFrame(obs, index=[f"c{i}" for i in range(1, 7)]),
            var=pd.DataFrame(index=["g1", "g2"]),
        )
        data.layers["counts"] = matrix.copy()
        data.write_h5ad(path)
        sidecar(path)


def make_feathers(analysis_root: Path) -> None:
    embeddings = analysis_root / "embeddings"
    embeddings.mkdir(parents=True, exist_ok=True)
    ids = [f"s{i}" for i in range(1, 7)]
    matrix = np.abs(np.subtract.outer(np.arange(6, dtype=float), np.arange(6, dtype=float)))
    for ds in ORDER:
        stem = f"{ds}_batch_effect_uncorrected"
        for suffix in ("pilot", "pilotgm", "mrvi", "qot"):
            path = embeddings / f"{stem}_hvg2000_highres_{suffix}_dists.feather"
            frame = pd.DataFrame(matrix, index=ids, columns=ids)
            frame.to_feather(path)
            sidecar(path)


def make_rds(analysis_root: Path) -> None:
    script = analysis_root.parent / "make_rds.R"
    script.write_text(
        """
        root <- Sys.getenv("ANALYSIS_ROOT")
        order <- strsplit(Sys.getenv("ORDER"), "|", fixed=TRUE)[[1]]
        excluded <- c("Alzheimer", "Diabetes", "Parkinson")
        ids <- paste0("s", 1:6)
        mat <- as.matrix(dist(seq_along(ids))); dimnames(mat) <- list(ids, ids)
        write_one <- function(path, names) {
          b <- setNames(lapply(names, function(n) list(dist_mat=mat)), names)
          saveRDS(b, path)
          d <- unname(tools::md5sum(path)); s <- file.info(path)$size
          writeLines(c(paste0("MD5=", d), paste0("SIZE=", s), paste0("PATH=", path)), paste0(path, ".md5"))
        }
        dir.create(file.path(root, "results"), recursive=TRUE, showWarnings=FALSE)
        for (ds in order) {
          stem <- paste0(ds, "_batch_effect_uncorrected")
          names <- COMPOSITION <- c("ECODA_authors_HR", "ECODA_HiTME_HR_layer2", "ECODA_scATOMIC_HR", "ECODA_seuratres_2", "ECODA_authors_HR_NULL")
          if (ds %in% excluded) names <- setdiff(names, c("ECODA_HiTME_HR_layer2", "ECODA_scATOMIC_HR"))
          write_one(file.path(root, "results", paste0(stem, "_composition.rds")), names)
          write_one(file.path(root, "results", paste0(stem, "_pseudobulk.rds")), "Pseudobulk_hvg2000")
          write_one(file.path(root, "results", paste0(stem, "_gloscope.rds")), "GloScope_hvg2000_pcadims30")
        }
        """
    )
    env = os.environ.copy()
    env.update({"ANALYSIS_ROOT": str(analysis_root), "ORDER": "|".join(ORDER)})
    subprocess.run(["pixi", "run", "Rscript", "--vanilla", str(script)], cwd=ROOT, env=env, check=True)

def run_builder(input_root: Path, analysis_root: Path, output_root: Path) -> subprocess.CompletedProcess[str]:
    selection = input_root.parent / "selection.tsv"
    command = [
        "pixi", "run", "Rscript", "--vanilla",
        str(ROOT / "notebooks/dataset_onboarding/build_batch_candidate_evidence.R"),
        "--selection-file", str(selection),
        "--analysis-root", str(analysis_root),
        "--input-root", str(input_root),
        "--output-dir", str(output_root),
    ]
    env = os.environ.copy()
    env["RETICULATE_PYTHON"] = str(ROOT / ".pixi/envs/default/bin/python")
    return subprocess.run(command, cwd=ROOT, text=True, capture_output=True, env=env)


def write_selection(input_root: Path) -> None:
    selection = input_root.parent / "selection.tsv"
    selection.write_text("".join(f"{ds}\tbatch_effect_uncorrected\n" for ds in ORDER))
    sidecar(selection)


def expect_failure(input_root: Path, analysis_root: Path, output_root: Path, text: str) -> None:
    result = run_builder(input_root, analysis_root, output_root)
    assert result.returncode != 0, result.stdout + result.stderr
    assert text.lower() in (result.stdout + result.stderr).lower(), result.stdout + result.stderr


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="ecoda-evidence-") as temporary:
        root = Path(temporary)
        input_root = root / "input"
        analysis_root = root / "analysis"
        output_root = root / "evidence"
        make_h5ads(input_root)
        make_feathers(analysis_root)
        make_rds(analysis_root)
        missing_selection = run_builder(input_root, analysis_root, output_root)
        assert missing_selection.returncode != 0, missing_selection.stdout + missing_selection.stderr
        missing_text = missing_selection.stdout + missing_selection.stderr
        assert "selection" in missing_text.lower(), missing_text
        assert not output_root.exists()
        write_selection(input_root)

        first = run_builder(input_root, analysis_root, output_root)
        assert first.returncode == 0, first.stdout + first.stderr
        review = pd.read_csv(output_root / "batch_candidate_review.csv")
        assert len(list(output_root.glob("*_batch_candidate_evidence.csv"))) == 12
        assert (output_root / "batch_candidate_review.csv.md5").is_file()
        assert bool(review.loc[(review.dataset == "Alzheimer") & (review.candidate == "assay"), "present"].eq(False).any())
        assert bool(review.loc[(review.dataset == "Breast_cancer") & (review.candidate == "assay"), "warnings"].str.contains("incomplete").any())
        assert bool(review.loc[(review.dataset == "Covid19_PBMC") & (review.candidate == "City"), "constant_candidate"].eq(True).any())
        assert bool(review.loc[(review.dataset == "Kidney_KPMP") & (review.candidate == "experiment"), "sample_unique_candidate"].eq(True).any())
        assert bool(review.loc[(review.dataset == "Myocardial_infarction") & (review.candidate == "batch"), "perfect_confounded"].eq(True).any())
        unavailable = review[(review.dataset == "Alzheimer") & (review.method == "ECODA_HiTME_HR_layer2")]
        assert not unavailable.empty and unavailable.method_available.eq(False).all()
        assert unavailable.method_reason.eq("not_suitable_for_auto_annotation").all()

        selection = input_root.parent / "selection.tsv"
        selection_lines = selection.read_text().splitlines()
        mutations = {
            "reordered": [selection_lines[1], selection_lines[0], *selection_lines[2:]],
            "duplicate": [selection_lines[0], selection_lines[0], *selection_lines[2:]],
            "wrong-view": [
                selection_lines[0].replace("batch_effect_uncorrected", "benchmark_analysis", 1),
                *selection_lines[1:],
            ],
            "extra-row": [*selection_lines, "Extra\tbatch_effect_uncorrected"],
        }
        for name, mutated_lines in mutations.items():
            selection.write_text("\n".join(mutated_lines) + "\n")
            sidecar(selection)
            invalid = run_builder(input_root, analysis_root, root / f"selection-{name}")
            assert invalid.returncode != 0, invalid.stdout + invalid.stderr
            assert not (root / f"selection-{name}").exists()
        selection.write_text("\n".join(selection_lines) + "\n")
        sidecar(selection)

        missing = analysis_root / "embeddings" / "Alzheimer_batch_effect_uncorrected_hvg2000_highres_pilot_dists.feather"
        missing.unlink(); Path(f"{missing}.md5").unlink()
        expect_failure(input_root, analysis_root, root / "missing-output", "missing applicable")
        make_feathers(analysis_root)

        h5ad = next((input_root / "Alzheimer" / "output").glob("*.h5ad"))
        h5ad.with_name(f"{h5ad.name}.md5").write_text("MD5=00000000000000000000000000000000\nSIZE=1\nPATH=" + str(h5ad) + "\n")
        expect_failure(input_root, analysis_root, root / "corrupt-output", "MD5 mismatch")
        sidecar(h5ad)

        duplicate = analysis_root / "embeddings" / "Alzheimer_batch_effect_uncorrected_hvg2000_highres_pilot_dists.feather"
        frame = pd.read_feather(duplicate)
        frame.index = ["s1", "s1", "s3", "s4", "s5", "s6"]
        frame.to_feather(duplicate); sidecar(duplicate)
        expect_failure(input_root, analysis_root, root / "duplicate-output", "duplicated")
        make_feathers(analysis_root)

        mismatch = analysis_root / "embeddings" / "Alzheimer_batch_effect_uncorrected_hvg2000_highres_pilot_dists.feather"
        frame = pd.read_feather(mismatch)
        frame.index = list(reversed(frame.index.tolist()))
        frame.to_feather(mismatch); sidecar(mismatch)
        expect_failure(input_root, analysis_root, root / "order-output", "sample-order mismatch")
    print("batch candidate evidence contracts: OK")


if __name__ == "__main__":
    main()
