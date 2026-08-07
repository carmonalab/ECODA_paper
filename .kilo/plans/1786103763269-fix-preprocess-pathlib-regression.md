# Fix Path/str regression in 3_scrnaseq_preprocessing

## Context

HPC run of `./src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh --ds_name _debug` (job 4289408, task 15) failed:

```
TypeError: unsupported operand type(s) for /: 'str' and 'str'
  File ".../1.1.1_preprocess.py", line 225, in main
    processed_file_path = output_dir / output_file_name
```

Verified root cause (confirmed against code + git history): commit `d24f73d` ("Simplify path handling in preprocess main; rename base_path to input_dir", Aug 4) removed the `resolve_path()` helper that converted `input_dir`/`output_dir` to `pathlib.Path` in `main()`. `1.1_run_worker.sh:51-56` now passes plain strings from bash (`${DATA_DIR}`, `${OUTPUT_DIR}`), so `main()` receives str and the `/` operator breaks. Not `_debug`-specific — would fail for any dataset.

Second instance of the same bug: `src/utils/preprocess_utils.py:34` (`input_path = input_dir / input_name` in `load_single_input()`) would crash next even after fixing `main()`, since the raw str flows through `load_input()` (called at `1.1.1_preprocess.py:231`).

Unaffected callers (checked): `src/2_dataset_specific_preprocessing/1.1.1_create_combinedpbmc_dataset.py` constructs its own `Path` objects (`main()` lines 217-228) before calling `load_input()` — no change needed there.

`config_path` stays a str: `read_datasets_json()` (`src/datasets_io.py:21`) uses `open(path)`, which accepts both.

## Changes

### 1. Required — fix the crash in `src/3_scrnaseq_preprocessing/1.1.1_preprocess.py` `main()`

`Path` is already imported (line 5). In `main()` (`1.1.1_preprocess.py:196`), immediately before `os.makedirs(output_dir, exist_ok=True)`:

```python
def main(config_path, input_dir, output_dir, ds_name=None, force=False):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    os.makedirs(output_dir, exist_ok=True)
```

This fixes both the line-225 crash and the `load_input()` path (Path objects flow down). Relative-path joining against `PROJECT_ROOT` (pre-d24f73d behavior) is intentionally NOT restored: the worker always passes absolute `${HPC_SCRATCH_DIR}` paths.

### 2. Optional hardening — `src/utils/preprocess_utils.py` `load_single_input()` (line 33)

Make the shared helper str-tolerant for any future caller (no behavior change for the existing combinedpbmc caller, which already passes Path):

```python
def load_single_input(input_name, input_dir, output_dir):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    input_path = input_dir / input_name
```

(Step 1 alone is sufficient to unblock the pipeline; step 2 is one line and removes the class of bug.)

## Validation

1. Local syntax check only (no pipeline run, per AGENTS.md): `python -m py_compile src/3_scrnaseq_preprocessing/1.1.1_preprocess.py src/utils/preprocess_utils.py`
2. On HPC (user): re-run `./src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh --ds_name _debug`. No `--force` needed — the failed run crashed before writing any output h5ad. Expect "Array Job finished. All tasks COMPLETED" and NAS sync.
3. If step 2 passes, full pipeline runs for other datasets are unaffected by design (Path conversion is semantics-preserving for absolute paths).

## Notes

- No docs (ARCHITECTURE.md) changes needed — this is a bug fix, not a flow change.
- Commit message suggestion if committing: `Fix str/Path regression in preprocess main (d24f73d)`
