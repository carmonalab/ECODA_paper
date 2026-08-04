# Simplify path handling in `1.1.1_preprocess.py` `main()`

## Context
- Sole caller of the script: `src/3_scrnaseq_preprocessing/1.1_run_worker.sh:39-43` — sources `slurm_config.sh` (exports `PROJECT_ROOT`, `slurm_config.sh:13`) and passes all four CLI args **as absolute paths**:
  - `--config_path "${DATASETS_JSON_FILE}"` (= `${PROJECT_ROOT}/datasets.json`)
  - `--base_path "${DATA_DIR}"` (= `${HPC_SCRATCH_DIR}/${DS_NAME}/data`)
  - `--output_dir "${OUTPUT_DIR}"` (= `${SCRATCH_OUTPUT_DIR}/${DS_NAME}`)
- `resolve_path`'s relative branch is a no-op in production; `main()`'s path defaults and the `Path(__file__)` fallback are dead code.
- `ds_name=None` is the sentinel for "process all datasets" (`1.1.1_preprocess.py:195`); the worker always overrides it with `--ds_name "${DS_NAME}"`.

## Change 1 — remove `resolve_path` and `project_root` entirely
`src/3_scrnaseq_preprocessing/1.1.1_preprocess.py:176-188` becomes:

```python
def main(config_path, input_dir, output_dir, ds_name=None):
    os.makedirs(output_dir, exist_ok=True)
```

- `resolve_path` deleted (all callers pass absolute paths; no relative-path support needed).
- `project_root` line deleted — it was only used by `resolve_path`.
- `config_path is not None` guard gone with `resolve_path`.
- `os` import still needed (line 227, `SAMPLE_COLNAME`); `Path` still needed (line 9, `sys.path` insert). No import changes.
- Body update: line 224 `load_input(input_file_name, base_path, output_dir)` → `load_input(input_file_name, input_dir, output_dir)`.

## Change 2 — rename `base_path` → `input_dir`
- `1.1.1_preprocess.py`: `main()` signature (Change 1), argparse flag `--base_path` → `--input_dir` (help: "Directory with raw input files").
- `1.1_run_worker.sh:41`: `--base_path "${DATA_DIR}"` → `--input_dir "${DATA_DIR}"`.
- `src/utils/preprocess_utils.py`: param `base_path` → `input_dir` in `load_single_input()` (line 30) and `load_input()` (line 43), incl. `input_path = base_path / input_name`. Both callers (`1.1.1_preprocess.py:224`, `_create_combinedpbmc_dataset.py:58`) call positionally — no call-site changes.
- `docs/ARCHITECTURE.md:81`: `--config_path/--base_path/--output_dir` → `--config_path/--input_dir/--output_dir`.

## Change 3 — argparse: path args `required=True`, keep `ds_name=None`
`__main__` block (lines 264-277):
- `--config_path`, `--input_dir`, `--output_dir`: `required=True`, no defaults (defaults would now silently resolve relative to CWD since `resolve_path` is gone). Drop "relative to PROJECT_ROOT unless absolute" wording.
- `--ds_name`: keep `default=None` (semantics: process all datasets; required by the loop filter at line 195; worker always overrides).

## Decisions / non-changes
- `ds_name=None` kept (requested check): `None` = "all datasets" sentinel; removing it would break the default run-all behavior for no benefit.
- `_create_combinedpbmc_dataset.py` untouched: own kebab-case `--base-path`/`--output-dir` convention, separate dataset-specific script.
- `TODO.md` lines 410/420 left as-is: historical completed-work log entries.
- `load_input()`/`read_datasets_json()` logic otherwise untouched.

## Validation
- `python -m py_compile src/3_scrnaseq_preprocessing/1.1.1_preprocess.py src/utils/preprocess_utils.py`
- `bash -n src/3_scrnaseq_preprocessing/1.1_run_worker.sh`
- No pipeline runs, per AGENTS.md rule; full pipeline validated later on HPC with the small debugging dataset.
