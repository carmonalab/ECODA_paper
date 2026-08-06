"""Python analog of src/utils/datasets_io.R::read_datasets_json().

Stdlib-only (json) so it can be imported from anywhere, including HPC
worker scripts that must not pull in R (rpy2) or heavy analysis deps.
"""

import json


def read_datasets_json(path="datasets.json", view=None):
    """Read datasets.json into a dict keyed by dataset name.

    Returns {ds_name: entry} with ALL datasets.json fields (future-proof).
    Entry summary fields (view_name, output_file, input_file, subset_vars)
    refer to the FIRST view matching `view` (with an output_file_name),
    mirroring the R version. `views` holds ALL matching views. Datasets
    without any matching view (e.g. Zhu, view-less raw source) are included
    with summary fields set to None/{} so raw inputs (file_names) are
    discoverable; the R version skips them (its callers use view filters).
    """
    with open(path) as f:
        datasets = json.load(f)

    result = {}
    for ds_name, ds in datasets.items():
        views = ds.get("views") or {}

        matched_views = {}
        for v_name, v in views.items():
            if view is not None and v_name != view:
                continue
            output_file = v.get("output_file_name")
            if not output_file:
                continue
            matched_views[v_name] = {
                "input_file": v.get("input_file_name"),
                "output_file": output_file,
                "subset_vars": v.get("subset_vars", {}),
            }

        first_name = next(iter(matched_views), None)
        first = matched_views.get(first_name)
        columns = ds.get("columns") or {}

        result[ds_name] = {
            # dataset-level fields (order mirrors datasets.json)
            "display_name": ds.get("display_name"),
            "file_names": ds.get("file_names"),
            "folder_name": ds.get("folder_name"),
            "tissue": ds.get("tissue"),
            "normal_tissue": ds.get("normal_tissue"),
            "use_for_benchmark": ds.get("use_for_benchmark"),
            "use_for_batch_effect": ds.get("use_for_batch_effect"),
            # columns (order mirrors datasets.json "columns")
            "sample_col": columns.get("sample"),
            "label_col": columns.get("label"),
            "batch_col": columns.get("batch"),
            "cell_type_low_res": columns.get("cell_type_low_res"),
            "cell_type_high_res": columns.get("cell_type_high_res"),
            "meta_cols_keep": ds.get("meta_cols_keep"),
            # first-matching-view summary (R-compat; mirrors view entry order)
            "view_name": first_name,
            "input_file": first["input_file"] if first else None,
            "output_file": first["output_file"] if first else None,
            "subset_vars": first["subset_vars"] if first else {},
            # all matching views
            "views": matched_views,
        }

    return result
