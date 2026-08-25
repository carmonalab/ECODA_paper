"""Authoritative metadata roles for the nine onboarding cohorts.

The registry is deliberately not inferred from column-name heuristics. Candidate
roles come from the full-file metadata audits; declared roles are user-confirmed
and the audits retain heuristic rankings as evidence only.
"""

from __future__ import annotations


_DATASET_KEYS = (
    "Alzheimer",
    "Breast_cancer",
    "Covid19_PBMC",
    "Diabetes",
    "Kidney_KPMP",
    "Lupus_PBMC",
    "Lung",
    "Myocardial_infarction",
    "Parkinson",
)


_ROLE_KEYS = {
    "sample",
    "label",
    "cell_type_low_res",
    "cell_type_high_res",
    "annotation_source",
}
_ANNOTATION_ROLES = {"low": "cell_type_low_res", "high": "cell_type_high_res"}


def _source(
    cells: int | None,
    independent_units: int | None,
    unit_kind: str,
    display_name: str,
    tissue: str,
    normal_tissue: bool,
    source_note: str,
) -> dict:
    return {
        "cells": cells,
        "independent_units": independent_units,
        "unit_kind": unit_kind,
        "display_name": display_name,
        "tissue": tissue,
        "normal_tissue": normal_tissue,
        "source_note": source_note,
    }


# Keep the order of candidates source-facing: it is used only for deterministic
# reporting and alias tie-breaking. It never overrides a failed gate or a
# user-confirmed registry role.
DATASET_SPECS: dict[str, dict] = {
    "Alzheimer": {
        "key": "Alzheimer",
        "file_name": "SEAAD_Alzheimer.h5ad",
        "expected_source": _source(
            1_395_601, 83, "donor", "Alzheimer (SEA-AD)", "Brain", False,
            "SEA-AD full-file metadata; source report 1,395,601 cells from 83 donors",
        ),
        "sample_candidates": ["donor_id", "Specimen ID"],
        "sample_stable_cols": ["Cognitive status", "sex", "disease", "tissue", "donor_id"],
        "bio_col": "Cognitive status",
        "batch_cols": ["assay", "tissue_type", "PMI"],
        "cell_type_candidates": ["cell_type", "Class", "Subclass", "Supertype"],
        "initial_registry_mode": "two_pass_batch_effect",
        "registry_roles": {
            "sample": "donor_id",
            "label": "Cognitive status",
            "cell_type_low_res": "Subclass",
            "cell_type_high_res": "Supertype",
            "annotation_source": {"low": "author", "high": "author"},
        },
        "decision_notes": [
            "User-selected Subclass is the canonical 24-label author tier between Class and Supertype.",
            "The prior Class selection is retained as heuristic evidence, not authority.",
        ],
        "not_suitable_for_auto_annotation": True,
        "subset_vars": {},
    },
    "Breast_cancer": {
        "key": "Breast_cancer",
        "file_name": "BreastCncr_processed.h5ad",
        "expected_source": _source(
            714_331, 167, "sample", "Breast cancer", "Breast tumor", False,
            "Kumar 2023 source report; user-confirmed sample_id unit has 167 samples",
        ),
        "sample_candidates": ["donor_id", "sample_id", "accSample"],
        "sample_stable_cols": ["disease", "sex", "tissue", "donor_id"],
        "bio_col": "disease",
        "batch_cols": [
            "assay", "sequencing_platform", "sample_preservation_method",
            "suspension_dissociation_time", "suspension_dissociation_reagent",
        ],
        "cell_type_candidates": ["broad_cell_type", "cell_type", "author_cell_type"],
        "initial_registry_mode": "two_pass_batch_effect",
        "registry_roles": {
            "sample": "sample_id",
            "label": "disease",
            "cell_type_low_res": "broad_cell_type",
            "cell_type_high_res": "author_cell_type",
            "annotation_source": {"low": "author", "high": "author"},
        },
        "decision_notes": [
            "Use the user-confirmed sample_id unit (167 samples), not the nested donor_id unit.",
            "The donor/sample nesting audit remains explanatory evidence.",
        ],
        "not_suitable_for_auto_annotation": False,
        "subset_vars": {},
    },
    "Covid19_PBMC": {
        "key": "Covid19_PBMC",
        "file_name": "Covid19_Ren2021.h5ad",
        "expected_source": _source(
            993_171, 151, "patient", "Covid-19 PBMC", "Blood", True,
            "Ren 2021 source report; PatientID is retained as an audit alternative",
        ),
        "sample_candidates": ["PatientID", "sampleID"],
        "sample_stable_cols": [
            "CoVID-19 severity", "Sex", "PatientID", "Outcome", "Sample time",
            "Sampling day (Days after symptom onset)",
        ],
        "bio_col": "CoVID-19 severity",
        "batch_cols": ["Single cell sequencing platform", "City", "datasets", "Sample type"],
        "cell_type_candidates": ["majorType", "celltype"],
        "initial_registry_mode": "two_pass_batch_effect",
        "registry_roles": {
            "sample": "sampleID",
            "label": "CoVID-19 severity",
            "cell_type_low_res": "majorType",
            "cell_type_high_res": "celltype",
            "annotation_source": {"low": "author", "high": "author"},
        },
        "decision_notes": [
            "Retain the user-confirmed sampleID and Covid-19 severity roles.",
            "The PatientID alternative remains in the audit as explanatory evidence.",
        ],
        "not_suitable_for_auto_annotation": False,
        "subset_vars": {},
    },
    "Diabetes": {
        "key": "Diabetes",
        "file_name": "diabetes.h5ad",
        "expected_source": _source(
            264_235, 52, "donor", "Diabetes (mouse)", "Pancreas", False,
            "Hrovatin 2023 source report; donor metadata conflicts are retained as warnings",
        ),
        "sample_candidates": ["donor_id", "dataset__design__sample"],
        "sample_stable_cols": ["disease", "sex", "donor_id", "strain", "diabetes_model", "age"],
        "bio_col": "disease",
        "batch_cols": ["batch_integration", "dataset", "design", "assay"],
        "cell_type_candidates": [
            "cell_type", "cell_type_reannotatedIntegrated", "cell_type_originalDataset",
        ],
        "initial_registry_mode": "two_pass_batch_effect",
        "registry_roles": {
            "sample": "donor_id",
            "label": "disease",
            "cell_type_low_res": "cell_type",
            "cell_type_high_res": "cell_type_reannotatedIntegrated",
            "annotation_source": {"low": "author", "high": "author"},
        },
        "decision_notes": [
            "Ten donor IDs mix sex metadata; retain this as a warning without changing the declared donor unit.",
            "The user decision overrides only the invariance selector, never ID or collision safety gates.",
        ],
        "not_suitable_for_auto_annotation": True,
        "subset_vars": {},
    },
    "Kidney_KPMP": {
        "key": "Kidney_KPMP",
        "file_name": "Kidney_KPMP.h5ad",
        "expected_source": _source(
            104_314, 47, "specimen", "Kidney (KPMP)", "Kidney", True,
            "KPMP source report; user-confirmed specimen unit has 47 specimens",
        ),
        "sample_candidates": ["donor_id", "specimen", "library"],
        "sample_stable_cols": [
            "condition.l1", "condition.l2", "condition.long", "sex", "donor_id",
            "region.l1", "region.l2", "tissue_type",
        ],
        "bio_col": "condition.l1",
        "batch_cols": [
            "experiment", "library", "tissue_type", "region.l1", "region.l2", "assay",
        ],
        "cell_type_candidates": [
            "major_cell_types", "subclass.l1", "subclass.l2", "subclass.l3",
            "subclass.full", "cell_type",
        ],
        "initial_registry_mode": "two_pass_batch_effect",
        "registry_roles": {
            "sample": "specimen",
            "label": "condition.l1",
            "cell_type_low_res": "subclass.l1",
            "cell_type_high_res": "subclass.l3",
            "annotation_source": {"low": "author", "high": "author"},
        },
        "decision_notes": [
            "Use the user-confirmed specimen unit (47 specimens) and the audited 16-to-67 hierarchy.",
            "Donor/library alternatives remain available for audit comparison only.",
        ],
        "not_suitable_for_auto_annotation": False,
        "subset_vars": {},
    },
    "Lupus_PBMC": {
        "key": "Lupus_PBMC",
        "file_name": "Lupus_Perez2022.h5ad",
        "expected_source": _source(
            1_263_676, 261, "sample", "Lupus PBMC", "Blood", True,
            "Perez 2022 source report; HiTME produces the configured annotation tiers",
        ),
        "sample_candidates": ["sampleID"],
        "sample_stable_cols": ["Status", "Sex", "Age", "pop_cov", "SLE_status"],
        "bio_col": "Status",
        "batch_cols": ["batch_cov", "Processing_Cohort", "ind_cov_batch_cov"],
        "cell_type_candidates": ["cell_types", "ct_cov", "layer1", "layer2"],
        "initial_registry_mode": "two_pass_batch_effect",
        "registry_roles": {
            "sample": "sampleID",
            "label": "Status",
            "cell_type_low_res": "layer1",
            "cell_type_high_res": "layer2",
            "annotation_source": {"low": "hitme", "high": "hitme"},
        },
        "decision_notes": [
            "HiTME must produce layer1/layer2; raw author cell_types/ct_cov are not fallbacks.",
            "Keep tissue exactly Blood so the established HiTME blood branch is used.",
            "Retain the two-sample Age conflict as a warning.",
        ],
        "not_suitable_for_auto_annotation": False,
        "subset_vars": {},
    },
    "Lung": {
        "key": "Lung",
        "file_name": "lungatlas.h5ad",
        "expected_source": _source(
            1_030_720, 336, "sample", "Lung atlas (10x subset)", "Lung", True,
            "Verified platform == 10x and assay-containing-10x masks are row-equivalent: 1,030,720 cells and 336 samples",
        ),
        "sample_candidates": ["donor_id", "sample"],
        "sample_stable_cols": [
            "origin", "origin_fine", "disease", "sex", "donor_id", "tissue",
        ],
        "bio_col": "disease",
        "batch_cols": ["dataset", "study", "platform", "assay", "origin_fine"],
        "cell_type_candidates": [
            "ann_coarse", "ann_fine", "cell_type_major", "cell_type", "cell_type_tumor",
        ],
        "initial_registry_mode": "two_pass_batch_effect",
        "registry_roles": {
            "sample": "sample",
            "label": "disease",
            "cell_type_low_res": "ann_coarse",
            "cell_type_high_res": "ann_fine",
            "annotation_source": {"low": "author", "high": "author"},
        },
        "decision_notes": [
            "The registry cohort is platform == 10x; verify filtered cell and sample counts before activation.",
            "ann_coarse to ann_fine must pass on the filtered cohort.",
            "The platform mask is authoritative if it differs from assay names containing 10x.",
            "The biological registry label is disease; origin remains a secondary field for tumor-versus-normal comparisons.",
        ],
        "not_suitable_for_auto_annotation": False,
        "subset_vars": {"platform": {"values": ["10x"], "op": "in"}},
    },
    "Myocardial_infarction": {
        "key": "Myocardial_infarction",
        "file_name": "Myocardial_Infarc_2.h5ad",
        "expected_source": _source(
            132_888, 24, "sample", "Myocardial infarction", "Heart", False,
            "Kuppe 2022 MI-2 source report; user-confirmed orig_ident has 24 identities",
        ),
        "sample_candidates": ["patient", "patient_region_id", "orig_ident", "accSample"],
        "sample_stable_cols": ["patient_group", "patient", "patient_region_id", "sampleType"],
        "bio_col": "patient_group",
        "batch_cols": ["batch", "sampleType", "dissociation_s1"],
        "cell_type_candidates": ["cell_type", "cell_subtype", "cell_subtype2"],
        "initial_registry_mode": "two_pass_batch_effect",
        "registry_roles": {
            "sample": "orig_ident",
            "label": "patient_group",
            "cell_type_low_res": "cell_type",
            "cell_type_high_res": "cell_subtype",
            "annotation_source": {"low": "author", "high": "author"},
        },
        "decision_notes": [
            "Use the user-confirmed orig_ident sample identity (24 identities), not patient_region_id.",
            "The prior patient-region ranking remains explanatory audit evidence.",
        ],
        "not_suitable_for_auto_annotation": False,
        "subset_vars": {},
    },
    "Parkinson": {
        "key": "Parkinson",
        "file_name": "Parkinson.h5ad",
        "expected_source": _source(
            2_096_155, 97, "donor", "Parkinson", "Brain", False,
            "Kamath 2022 source report; Leiden res-5 is produced downstream",
        ),
        "sample_candidates": ["donor_id"],
        "sample_stable_cols": ["disease", "sex", "tissue", "tissue_type", "donor_id"],
        "bio_col": "disease",
        "batch_cols": ["Brain_bank", "assay", "tissue_type", "PMI", "RIN"],
        "cell_type_candidates": [
            "derived_class2", "cell_type", "leiden_res_5_batch_effect_uncorrected_hvg2000",
        ],
        "initial_registry_mode": "two_pass_batch_effect",
        "registry_roles": {
            "sample": "donor_id",
            "label": "disease",
            "cell_type_low_res": "cell_type",
            "cell_type_high_res": "leiden_res_5_batch_effect_uncorrected_hvg2000",
            "annotation_source": {"low": "author", "high": "leiden"},
        },
        "decision_notes": [
            "Leiden res-5 is a produced high-resolution role; the corrected view substitutes its Harmony-qualified column.",
            "Ninety-six donors span more than one tissue value; retain this aggregation warning while disease remains the ground-truth label.",
        ],
        "not_suitable_for_auto_annotation": True,
        "subset_vars": {},
    },
}
DATASET_SPECS["Lung"]["expected_source"]["source_match_required"] = True


DEBUG_SPEC = {
    "key": "_debug",
    "file_name": "JoaI_2022_35773407_debug_5samples_batch_effect_analysis_uncorrected_ECODAprocessed.h5ad",
    "expected_source": {
        "cells": 6000,
        "independent_units": 12,
        "unit_kind": "sample",
        "display_name": "Debug (Joanito 5-sample subset)",
        "tissue": "Colon",
        "normal_tissue": False,
        "source_note": "Deterministic five-sample pipeline reference",
    },
    "sample_candidates": ["sample.ID"],
    "sample_stable_cols": ["sample.origin", "Gender", "sample.ID"],
    "bio_col": "sample.origin",
    "batch_cols": ["seqtec", "Site"],
    "cell_type_candidates": ["cell.type"],
    "initial_registry_mode": "debug_reference",
    "not_suitable_for_auto_annotation": False,
    "subset_vars": {},
}


def get_dataset_spec(name: str) -> dict:
    """Return a configured onboarding spec or fail closed."""
    if name == "_debug":
        return DEBUG_SPEC
    try:
        return DATASET_SPECS[name]
    except KeyError as exc:
        raise KeyError(
            f"Unknown onboarding dataset {name!r}; expected one of "
            f"{list(DATASET_SPECS)} or '_debug'"
        ) from exc


if set(DATASET_SPECS) != set(_DATASET_KEYS):
    raise RuntimeError("DATASET_SPECS keys do not match the nine onboarding cohorts")

for _name, _spec in DATASET_SPECS.items():
    _expected = {
        "key", "file_name", "expected_source", "sample_candidates", "sample_stable_cols",
        "bio_col", "batch_cols", "cell_type_candidates", "initial_registry_mode",
        "registry_roles", "decision_notes", "not_suitable_for_auto_annotation", "subset_vars",
    }
    if set(_spec) != _expected:
        raise RuntimeError(f"Unexpected fields in DATASET_SPECS[{_name!r}]: {set(_spec) ^ _expected}")
    if _spec["initial_registry_mode"] != "two_pass_batch_effect":
        raise RuntimeError(f"{_name}: initial_registry_mode must be two_pass_batch_effect")
    if set(_spec["registry_roles"]) != _ROLE_KEYS:
        raise RuntimeError(f"{_name}: registry_roles keys must be {_ROLE_KEYS}")
    sources = _spec["registry_roles"]["annotation_source"]
    if set(sources) != set(_ANNOTATION_ROLES):
        raise RuntimeError(f"{_name}: annotation_source must declare low and high sources")
    for role in ("sample", "label", "cell_type_low_res", "cell_type_high_res"):
        if not isinstance(_spec["registry_roles"][role], str) or not _spec["registry_roles"][role]:
            raise RuntimeError(f"{_name}: registry role {role} must be a non-empty column name")
