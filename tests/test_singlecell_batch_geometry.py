"""Tiny synthetic parent only: read boundary, native LISI oracle and real UMAP."""
import copy
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import tempfile
import subprocess
import shutil
import sys
import unittest
from unittest.mock import patch

os.environ.setdefault("NUMBA_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.neighbors import NearestNeighbors

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("singlecell_pilot", ROOT / "src/notebooks/run_singlecell_batch_geometry.py")
p = importlib.util.module_from_spec(spec); spec.loader.exec_module(p)


class PilotTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(); self.root = Path(self.temp.name).resolve()
        self.n = 160; rng = np.random.default_rng(21)
        self.pcs = rng.normal(size=(self.n, 50)).astype("float32")
        ids = ["cell_" + str(i) for i in range(self.n)]
        samples = ["sample_" + str(i % 33) for i in range(self.n)]
        obs = pd.DataFrame(dict(Sample=samples, sample_id=samples,
            broad_cell_type=pd.Categorical(["type" + str(i % 4) for i in range(self.n)]),
            disease=pd.Categorical(["display_only"] * self.n), assay=pd.Categorical(["10x 3' v3"] * self.n),
            suspension_dissociation_time=pd.Categorical(["12 hour" if i % 2 else "15 minute" for i in range(self.n)]),
            **{p.LEIDEN: pd.Categorical([str(i % 3) for i in range(self.n)])}), index=ids)
        _, ix = NearestNeighbors(n_neighbors=16).fit(self.pcs).kneighbors(self.pcs)
        g = sparse.csr_matrix((np.ones(self.n * 15, dtype="float32"), (np.repeat(np.arange(self.n), 15), ix[:, 1:].ravel())), shape=(self.n, self.n))
        g = (g + g.T).tocsr(); g.data[:] = 1
        self.parent = self.root / "Breast.h5ad"
        a = ad.AnnData(X=None, obs=obs, obsm={p.PC: self.pcs}, obsp={p.NEIGHBORS + "_connectivities": g},
            uns={p.NEIGHBORS: dict(connectivities_key=p.NEIGHBORS + "_connectivities", distances_key=p.NEIGHBORS + "_distances",
                 params=dict(use_rep=p.PC, n_pcs=50, n_neighbors=15, metric="euclidean", method="umap", random_state=0))})
        a.write_h5ad(self.parent)
        with h5py.File(self.parent, "a") as h:
            h.create_dataset("X", shape=(self.n, 3), dtype="float32", fillvalue=123)
            h.create_group("raw").create_dataset("X", shape=(self.n, 3), dtype="float32", fillvalue=456)
            h["layers"].create_dataset("counts", shape=(self.n, 3), dtype="float32", fillvalue=789)
        self.alz = self.root / "Alzheimer.h5ad"
        with h5py.File(self.alz, "w") as h:
            h.create_group("obsm").create_dataset(p.PC, shape=(1305022, 50), dtype="float32", chunks=(1, 50))
            g = h.create_group("obsp").create_group(p.NEIGHBORS + "_connectivities")
            for name, shape, dtype in [("data", (30183386,), "float32"), ("indices", (30183386,), "int32"), ("indptr", (1305023,), "int64")]:
                g.create_dataset(name, shape=shape, dtype=dtype, chunks=(100,))
        config = json.loads((ROOT / "datasets.json").read_text())
        (self.root / "datasets.json").write_text(json.dumps(config))
        artifacts = []
        for ds, path in [(p.DATASET, self.parent), ("Alzheimer", self.alz)]:
            proof = dict(PATH=str(path), MD5=p.digest(path, "md5"), SIZE=str(path.stat().st_size))
            record = dict(proof, STATE="PUBLISHED", PRODUCER="stage3", RUN_ID="stage3_confounds_20261108_v3")
            owner = dict(PATH=str(path), STATE="OK", STAGE="stage3", RUN_ID=record["RUN_ID"])
            w = dict(path=str(path), size=path.stat().st_size, md5=proof["MD5"], record_fields=record, owner_fields=owner)
            for key, values in [("proof", proof), ("record", record), ("owner", owner)]:
                file = self.root / (ds + "." + key); file.write_text("".join(k + "=" + v + "\n" for k, v in values.items()))
                w[key] = dict(path=str(file), sha256=p.digest(file), size=file.stat().st_size)
            artifacts.append(dict(dataset=ds, view=p.VIEW, scratch=w))
        self.receipt = self.root / "receipt.json"; self.receipt.write_text(json.dumps(dict(artifacts=artifacts)))
        self.review = self.root / "review.json"; self.review.write_text(json.dumps(dict(status="APPROVED_SELECTED_ARTIFACT_REUSE",
            receipt_sha256=p.digest(self.receipt), original_gate_status="FAILED", failed_gate_released=False)))
        self.control = dict(schema="ecoda_singlecell_breast_unc_pilot_v1", stage3_receipt=dict(path=str(self.receipt), sha256=p.digest(self.receipt)),
            stage3_review=dict(path=str(self.review), sha256=p.digest(self.review)),
            parent_stats={a["dataset"]: {k: getattr(Path(a["scratch"]["path"]).stat(), k) for k in ("st_dev", *p.STAT)} for a in artifacts})
        self.stack = __import__("contextlib").ExitStack()
        self.stack.enter_context(patch.object(p, "RECEIPT_SHA", p.digest(self.receipt)))
        self.stack.enter_context(patch.object(p, "REVIEW_SHA", p.digest(self.review)))
        self.stack.enter_context(patch.object(p, "PILOT_CELLS", self.n))
        self.stack.enter_context(patch.dict(os.environ, ECODA_SOURCE_ROOT=str(self.root), SLURM_CPUS_PER_TASK="1"))
        read_text = Path.read_text
        def fake_mount(file, *args, **kwargs):
            if str(file) == "/proc/self/mountinfo": return "1 0 0:1 / / rw - beegfs synthetic rw\n"
            return read_text(file, *args, **kwargs)
        self.stack.enter_context(patch.object(Path, "read_text", fake_mount))

    def tearDown(self):
        self.stack.close(); self.temp.cleanup()

    def test_native_lisi_oracle_and_missing_population(self):
        from harmonypy.lisi import compute_lisi
        d, ix = NearestNeighbors(n_neighbors=91).fit(self.pcs).kneighbors(self.pcs)
        inds = []; distances = []
        for i in range(self.n):
            keep = ix[i] != i; inds.append(ix[i][keep][:90]); distances.append(d[i][keep][:90])
        ix = np.asarray(inds); d = np.asarray(distances); labels = np.arange(self.n) % 4
        actual = p.lisi_neighbors(ix, d, labels, 4)
        expected = compute_lisi(self.pcs, pd.DataFrame(dict(label=labels)), ["label"])[:, 0]
        np.testing.assert_allclose(actual, expected, rtol=1e-5, atol=1e-5)
        raw, scaled, summary = p.score_population(self.pcs, np.array(["one"] * self.n), np.array(["sample"] * self.n), "mixing", 1, (ix, d))
        np.testing.assert_allclose(raw, 1); self.assertTrue(np.isnan(scaled).all()); self.assertIsNone(summary["scaled_cellwise_mean"])
        # Only valid labels enter the evaluation population; geometry is untouched.
        labels = np.array(["a" if i % 2 else "b" for i in range(self.n)], dtype=object); labels[:10] = None
        raw, _, summary = p.score_population(self.pcs, labels, np.arange(self.n).astype(str), "preservation", 1)
        self.assertTrue(np.isnan(raw[:10]).all()); self.assertEqual(summary["evaluable_cells"], 150)

    def test_real_umap_expression_boundary_and_atomic_products(self):
        import scanpy as sc
        get = h5py.Group.__getitem__; opened = []
        def safe(group, key):
            if Path(group.file.filename) == self.parent:
                self.assertFalse(group.name == "/" and str(key).strip("/").split("/")[0] in ("X", "raw", "layers"), "Parent expression object opened")
                opened.append(str(key))
            return get(group, key)
        output = self.root / "output"
        with patch.object(h5py.Group, "__getitem__", safe), patch.object(sc.pp, "neighbors", side_effect=AssertionError("Graph refit forbidden")), patch.object(sc.pp, "pca", side_effect=AssertionError("PCA refit forbidden")):
            p.produce(self.control, output, "synthetic_pilot")
            p.validate_output(self.control, output, "synthetic_pilot")
        self.assertTrue(opened); self.assertEqual({x.name for x in output.iterdir()}, set(p.OUTPUTS) | {"files.tsv"})
        with h5py.File(output / "embedding.h5ad") as h:
            self.assertTrue("X" not in h or h["X"].shape is None)
            self.assertEqual(len(h["layers"]), 0)
        phases = json.loads((output / "phases.json").read_text()); self.assertIn("ANN90", phases)
        with self.assertRaises(ValueError): p.produce(self.control, output, "synthetic_pilot")
        with self.assertRaises(FileExistsError):
            with p.atomic_file(output / "metrics.json") as f: f.write_text("must not overwrite")

    def test_guard_before_values_and_alzheimer_header_only(self):
        broken = copy.deepcopy(self.control); broken["parent_stats"][p.DATASET]["st_size"] += 1
        with patch.object(h5py, "File", side_effect=AssertionError("H5 before native guard")):
            with self.assertRaises(ValueError): p.produce(broken, self.root / "bad", "bad")
        with patch.object(h5py.Dataset, "__getitem__", side_effect=AssertionError("Header mode read values")):
            result = p.memory_preflight(self.control)
        self.assertEqual(result["cells"], 1305022); self.assertEqual(result["coordinate_values_read"], 0)
        self.assertGreaterEqual(result["conservative_provisional_allocation_bytes"], 128 * 1024**3)
        self.assertFalse(result["full_Alzheimer_compute_approved"])

    def test_native_scalar_spool_ownership_and_failure(self):
        tree = self.root / "snapshot/tree"; utils = tree / "src/utils/bash"; utils.mkdir(parents=True)
        (tree / "src/notebooks").mkdir(); scratch = self.root / "scratch"; inputs = scratch / "_ecoda_inputs"; inputs.mkdir(parents=True)
        for name in ["singlecell_batch_geometry.sh", "ecoda_run_common.sh"]: shutil.copy2(ROOT / "src/utils/bash" / name, utils / name)
        # The real numerical worker is exercised above. Here only worker computation
        # is stubbed; run/records/checksums/owners and copied-spool bootstrap are native.
        (tree / "src/notebooks/run_singlecell_batch_geometry.py").write_text('''import argparse,hashlib,pathlib
p=argparse.ArgumentParser();p.add_argument('--control');p.add_argument('--control-sha256');p.add_argument('--output');p.add_argument('--run-id');p.add_argument('--check-control',action='store_true');p.add_argument('--validate-output',action='store_true');a=p.parse_args()
assert hashlib.sha256(pathlib.Path(a.control).read_bytes()).hexdigest()==a.control_sha256
if a.check_control or a.validate_output:raise SystemExit(0)
out=pathlib.Path(a.output);out.mkdir()
files=['embedding.h5ad','metrics.json','display.feather','pilot_pca.pdf','pilot_umap.pdf','phases.json']
for name in files:(out/name).write_bytes(b'synthetic worker payload')
with (out/'files.tsv').open('x')as f:
 for name in files:f.write(str(out/name)+'\\t'+hashlib.md5((out/name).read_bytes()).hexdigest()+'\\n')
''')
        (tree / "source.manifest").write_text("SOURCE_ROOT=" + str(tree) + "\n")
        (tree / "src/slurm_config.sh").write_text('''export PROJECT_ROOT="${TEST_TREE}" HPC_SCRATCH_DIR="${TEST_SCRATCH}"
export ECODA_RUNS_ROOT="${TEST_SCRATCH}/_ecoda_runs" ECODA_OWNERS_ROOT="${TEST_SCRATCH}/_ecoda_owners" PYTHON_BIN="${TEST_PYTHON}"
SLURM_PARTITION_BENCHMARK_CPU=shared-cpu
''')
        (utils / "ecoda_runtime.sh").write_text('''ecoda_runtime_validate_submission() { printf 'fixed fixture runtime\\n' > "${ECODA_RUN_ROOT}/manifests/runtime.identity"; }
ecoda_runtime_validate_bound_run() { test -s "${ECODA_RUN_ROOT}/manifests/runtime.identity"; }
ecoda_runtime_export_csv() { printf ECODA_RUNTIME_PROFILE=stage5; }
ecoda_runtime_reexec_worker() { test "${FAIL_WORKER:-0}" != 1; }
''')
        binpath = self.root / "bin"; binpath.mkdir()
        sbatch = binpath / "sbatch"; sbatch.write_text('''#!/bin/bash
set -euo pipefail
[[ "$*" == *'--cpus-per-task=8'* && "$*" == *'--mem=64G'* && "$*" == *'--time=02:00:00'* && "$*" == *'--partition=shared-cpu'* && "$*" != *'--constraint'* && "$*" != *'--gpus'* ]]
args=("$@");spool="${TEST_SCRATCH}/spool/${ECODA_RUN_ID}";mkdir -p "${spool}"
cp "${args[${#args[@]}-2]}" "${spool}/slurm_script"
set +e
SLURM_CPUS_PER_TASK=8 /bin/bash "${spool}/slurm_script" --worker > "${ECODA_RUN_ROOT}/logs/worker.log" 2> "${ECODA_RUN_ROOT}/logs/worker.err"
rc=$?;set -e;[[ ${rc} == 0 ]] || cat "${ECODA_RUN_ROOT}/logs/worker.err" >&2
printf '900001;fixture\\n';exit "${rc}"
'''); sbatch.chmod(0o755)
        for name in ["squeue", "sacct"]:
            file = binpath / name; file.write_text("#!/bin/bash\nexit 99\n"); file.chmod(0o755)
        receipt = inputs / "receipt.json"; review = inputs / "review.json"; receipt.write_text("{}"); review.write_text("{}")
        control = inputs / "control.json"; control.write_text(json.dumps(dict(stage3_receipt=dict(path=str(receipt)), stage3_review=dict(path=str(review)))))
        env = dict(os.environ, TEST_TREE=str(tree), TEST_SCRATCH=str(scratch), TEST_PYTHON=sys.executable,
                   ECODA_SOURCE_ROOT=str(tree), ECODA_SOURCE_MANIFEST=str(tree / "source.manifest"), ECODA_SOURCE_SNAPSHOT_REQUIRED="1", PATH=str(binpath) + os.pathsep + os.environ["PATH"])
        cmd = ["/bin/bash", str(utils / "singlecell_batch_geometry.sh")]
        def invoke(run, fail=False, validate=False):
            e = dict(env, ECODA_RUN_ID=run, ECODA_RUN_ROOT=str(scratch / "_ecoda_runs" / run), FAIL_WORKER="1" if fail else "0")
            args = ["--validate"] if validate else ["--control", str(control), "--control-sha256", hashlib.sha256(control.read_bytes()).hexdigest()]
            return subprocess.run(cmd + args, env=e, capture_output=True, text=True)
        good = invoke("pilot_good"); self.assertEqual(good.returncode, 0, good.stderr)
        run = scratch / "_ecoda_runs/pilot_good"; self.assertIn("STATE=OK", (run / "status/terminal").read_text())
        self.assertEqual((run / "manifests/scheduler_ids.tsv").read_text(), "ARRAY\t900001\n")
        self.assertEqual(invoke("pilot_good", validate=True).returncode, 0)
        self.assertNotEqual(invoke("pilot_good").returncode, 0)
        # A wrong-stage completed owner is rejected read-only and never repaired.
        owners = list((scratch / "_ecoda_owners/artifact").glob("*/owner"))
        own = next(x for x in owners if "PATH=" + str(run / "artifacts/singlecell_batch_geometry/files.tsv") + "\n" in x.read_text())
        own.write_text(own.read_text().replace("STAGE=stage5", "STAGE=stage3")); before = own.read_bytes()
        status = (run / "status/singlecell_pilot.status").read_bytes(); terminal = (run / "status/terminal").read_bytes()
        self.assertNotEqual(invoke("pilot_good", validate=True).returncode, 0); self.assertEqual(own.read_bytes(), before)
        self.assertEqual((run / "status/singlecell_pilot.status").read_bytes(), status); self.assertEqual((run / "status/terminal").read_bytes(), terminal)
        failed = invoke("pilot_failed", fail=True); self.assertNotEqual(failed.returncode, 0)
        self.assertEqual((scratch / "_ecoda_runs/pilot_failed/manifests/scheduler_ids.tsv").read_text(), "ARRAY\t900001\n")
        self.assertIn("STATE=FAIL", (scratch / "_ecoda_runs/pilot_failed/status/terminal").read_text())

    def bulk_control(self):
        """Native-format12 metadata parents; only Breast COR values are used."""
        artifacts = []; stats = {d: {} for d in p.DATASETS}
        corrected = self.root / "Breast_COR.h5ad"; shutil.copy2(self.parent, corrected)
        settings = p.row_settings(p.DATASET, "cor_harmony")
        with h5py.File(corrected, "a") as h:
            h.move("obsm/" + p.PC, "obsm/" + settings["pc"])
            h.move("obsp/" + p.NEIGHBORS + "_connectivities", "obsp/" + settings["neighbors"] + "_connectivities")
            h.move("uns/" + p.NEIGHBORS, "uns/" + settings["neighbors"])
            g = h["uns"][settings["neighbors"]]
            for key, value in [("connectivities_key", settings["neighbors"] + "_connectivities"), ("distances_key", settings["neighbors"] + "_distances")]: g[key][()] = value
            g["params/use_rep"][()] = settings["pc"]
            h.move("obs/" + p.LEIDEN, "obs/" + settings["leiden"])
        for d in p.DATASETS:
            for rep in p.REPRESENTATIONS:
                s = p.row_settings(d, rep)
                path = self.parent if (d, rep) == (p.DATASET, "unc_raw") else corrected if d == p.DATASET else self.root / (d + rep + ".h5ad")
                if not path.exists(): path.write_bytes(b"synthetic nonselected parent; no array reads")
                proof = dict(PATH=str(path), MD5=p.digest(path, "md5"), SIZE=str(path.stat().st_size))
                record = dict(proof, STATE="PUBLISHED", PRODUCER="stage3", RUN_ID="stage3_confounds_20261108_v3")
                owner = dict(PATH=str(path), STATE="OK", STAGE="stage3", RUN_ID=record["RUN_ID"])
                w = dict(path=str(path), size=path.stat().st_size, md5=proof["MD5"], record_fields=record, owner_fields=owner)
                for key, values in [("proof", proof), ("record", record), ("owner", owner)]:
                    file = self.root / (d + rep + "." + key); file.write_text("".join(k + "=" + v + "\n" for k, v in values.items()))
                    w[key] = dict(path=str(file), sha256=p.digest(file), size=file.stat().st_size)
                artifacts.append(dict(dataset=d, view=s["view"], scratch=w))
                stats[d][s["view"]] = {k: getattr(path.stat(), k) for k in ("st_dev", *p.STAT)}
        receipt = self.root / "receipt12.json"; receipt.write_text(json.dumps(dict(artifacts=artifacts)))
        review = self.root / "review12.json"; review.write_text(json.dumps(dict(status="APPROVED_SELECTED_ARTIFACT_REUSE", receipt_sha256=p.digest(receipt), original_gate_status="FAILED", failed_gate_released=False)))
        return dict(schema="ecoda_singlecell_batch_main_v2", stage3_receipt=dict(path=str(receipt), sha256=p.digest(receipt)), stage3_review=dict(path=str(review), sha256=p.digest(review)), parent_stats=stats)

    def test_bulk_fixed_scope_corrected_and_no_hpc_render(self):
        control = self.bulk_control(); output = self.root / "bulk_COR"
        counts = {d: (self.n, 33) for d in p.DATASETS}
        with patch.object(p, "RECEIPT_SHA", control["stage3_receipt"]["sha256"]), patch.object(p, "REVIEW_SHA", control["stage3_review"]["sha256"]), patch.object(p, "COHORT_SIZE", counts):
            with patch.object(p, "render", side_effect=AssertionError("No bulk HPC PDF")):
                p.produce(control, output, "synthetic_bulk", p.DATASET, "cor_harmony")
                p.validate_output(control, output, "synthetic_bulk", p.DATASET, "cor_harmony")
            self.assertEqual({x.name for x in output.iterdir()}, set(p.BULK_OUTPUTS) | {"files.tsv"})
            metrics = json.loads((output / "metrics.json").read_text())
            self.assertFalse(metrics["HPC_rendering_performed"]); self.assertNotIn("rendering", json.loads((output / "phases.json").read_text()))
            self.assertEqual(metrics["source_columns"]["author"], "broad_cell_type")
            with h5py.File(output / "embedding.h5ad") as h:
                self.assertNotIn("author_cell_type", h["obs"])  # no misleading HR-name alias in H5
                self.assertIn("broad_cell_type", h["obs"])
            with self.assertRaises(ValueError): p.produce(control, self.root / "forbidden", "x", p.DATASET, "unc_raw")
            bad = copy.deepcopy(control); bad["parent_stats"]["Stephenson"]["batch_effect_corrected"]["st_size"] += 1
            with patch.object(h5py, "File", side_effect=AssertionError("Before all12 native guards")):
                with self.assertRaises(ValueError): p.produce(bad, self.root / "bad12", "x", p.DATASET, "cor_harmony")

    def test_batch_balanced_display_determinism(self):
        ids = ["id" + str(i) for i in range(221)]
        labels = np.array(["rare"] + ["small"] * 2 + ["large"] * 218, dtype=object)
        selected, policy = p.display_indices(ids, labels, limit=20, full_threshold=100)
        self.assertEqual(len(selected), 20); self.assertEqual(len(set(selected)), 20); self.assertTrue({0, 1, 2} <= set(selected))
        order = np.random.default_rng(3).permutation(len(ids))
        reordered, _ = p.display_indices([ids[i] for i in order], labels[order], limit=20, full_threshold=100)
        self.assertEqual([ids[i] for i in selected], [ids[order[i]] for i in reordered])
        self.assertEqual(sum(x["selected"] for x in policy["strata"]), 20)
        np.testing.assert_array_equal(p.display_indices(ids[:50], labels[:50])[0], np.arange(50))

    def test_local_main_renderer_v1_pilot_reuse_single_scatter(self):
        import matplotlib.axes
        import matplotlib.figure
        rows = []; source_root = "/synthetic/native"; review_path = self.root / "geometry.review.json"; reviewed = []
        for d in p.DATASETS:
            for rep in p.REPRESENTATIONS:
                settings = p.row_settings(d, rep)
                frame = pd.DataFrame(dict(cell_id=["c" + str(i) for i in range(self.n)], Sample=["s" + str(i % (45 if d == "Alzheimer" else 3)) for i in range(self.n)],
                    author_cell_type=["cell" + str(i % 2) for i in range(self.n)], biological_label=["bio"] * self.n,
                    batch_label=["batch" + str(i % 2) for i in range(self.n)], leiden2=[str(i % 3) for i in range(self.n)],
                    PC1=np.arange(self.n, dtype=float), PC2=np.arange(self.n, dtype=float) / 2, UMAP1=np.arange(self.n, dtype=float) / 3, UMAP2=np.arange(self.n, dtype=float) / 4))
                legacy = (d, rep) == (p.DATASET, "unc_raw")
                if legacy: frame = frame.rename(columns=dict(author_cell_type="broad_cell_type", biological_label="disease", batch_label="suspension_dissociation_time", leiden2=p.LEIDEN))
                f = self.root / (d + rep + ".feather"); frame.to_feather(f)
                metric = dict(schema="ecoda_singlecell_pilot_metrics_v1" if legacy else "ecoda_singlecell_main_metrics_v2", dataset=d, run_id="synthetic", view=settings["view"],
                    representation=settings["pc"], representation_id=rep, cells=self.n, dimensions=50, n_neighbors=90, perplexity=30, source_expression_read=False,
                    source_columns=settings, batch=dict(scaled_cellwise_mean=.123456), cell_type=dict(scaled_cellwise_mean=.987654),
                    display_contract=dict(rows=self.n, ordered_cell_id_sha256=p.ids_digest(frame.cell_id), md5=p.digest(f, "md5")))
                m = self.root / (d + rep + ".json"); m.write_text(json.dumps(metric))
                refs = {"display": dict(path=str(f), md5=p.digest(f, "md5"), source_path=source_root + "/" + f.name), "metrics": dict(path=str(m), md5=p.digest(m, "md5"), source_path=source_root + "/" + m.name)}
                reviewed.extend(dict(path=r["source_path"], md5=r["md5"]) for r in refs.values())
                rows.append(dict(dataset=d, representation=rep, **refs))
        review_path.write_text(json.dumps(dict(schema="ecoda_singlecell_geometry_root_review_v1", status="APPROVED_SYNTHETIC_FIXTURE_ONLY", reviewer="TEST", reviewed_at="TEST", source_expression_read=False, run_id="synthetic", output_artifacts=reviewed)))
        for row in rows: row["review"] = dict(path=str(review_path), sha256=p.digest(review_path))
        control = dict(schema="ecoda_singlecell_main_display_control_v1", status="APPROVED_LOCAL_MAIN_DISPLAY_EXTRACTS", rows=rows)
        calls = []; sample_legends = []; scatter = matplotlib.axes.Axes.scatter; save = matplotlib.figure.Figure.savefig
        def capture(ax, x, y, **kwargs): calls.append((np.array(x), np.asarray(kwargs["c"]))); return scatter(ax, x, y, **kwargs)
        def check(fig, *args, **kwargs):
            self.assertTrue(all(len(ax.collections) == 1 for ax in fig.axes))
            self.assertTrue(all(any("50D mean iLISI mixing=0.1235" in t.get_text() and "cLISI preservation=0.9877" in t.get_text() for t in ax.texts) for ax in fig.axes))
            for ax in fig.axes:
                if ax.get_title() == "Sample": sample_legends.append({t.get_text() for t in ax.get_legend().get_texts()})
            return save(fig, *args, **kwargs)
        with patch.object(p, "COHORT_SIZE", {d: (self.n, 33) for d in p.DATASETS}), patch.object(matplotlib.axes.Axes, "scatter", capture), patch.object(matplotlib.figure.Figure, "savefig", check), patch.object(p, "ann90", side_effect=AssertionError("Local ANN forbidden")), patch.object(h5py, "File", side_effect=AssertionError("Local H5 forbidden")):
            p.render_main(control, self.root / "rendered")
        self.assertEqual(len(calls), 108); np.testing.assert_array_equal(calls[0][0], calls[4][0]); np.testing.assert_array_equal(calls[0][1], calls[4][1])
        report = json.loads((self.root / "rendered/render_manifest.json").read_text()); self.assertEqual(len(report["figures"]), 12)
        self.assertTrue(all(x["score_computations"] == 0 for x in report["figures"]))
        self.assertTrue(any(labels == {"s" + str(i) for i in range(45)} for labels in sample_legends))

    def test_display_caption_artist_only(self):
        import ast
        import inspect
        from matplotlib.figure import Figure
        # Execute only the actual caption statements, with no panels, PDFs or geometry.
        tree = ast.parse(inspect.getsource(p.render_main))
        statements = [node for node in ast.walk(tree) if
            isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == "display_policy" for t in node.targets)
            or isinstance(node, ast.Expr) and isinstance(node.value, ast.Call) and isinstance(node.value.func, ast.Attribute) and node.value.func.attr == "suptitle"]
        self.assertEqual(len(statements), 2)
        code = compile(ast.fix_missing_locations(ast.Module(body=sorted(statements, key=lambda n: n.lineno), type_ignores=[])), "caption", "exec")
        for dataset, displayed, total, policy in [("Breast_cancer", 114579, 114579, "all cells"), ("Alzheimer", 100000, 1305022, "batch-balanced; not original cell abundances")]:
            fig = Figure(); exec(code, dict(fig=fig, dataset=dataset, ids=range(displayed), summaries=[dict(cells=total)]))
            text = fig._suptitle.get_text()
            self.assertIn(f"Display: {displayed:,}/{total:,} cells per row ({policy})", text)
            self.assertIn("scores: full-population 50D cellwise means", text)

    def test_small_category_palette_artist_only(self):
        import ast
        import inspect
        import matplotlib
        import matplotlib.pyplot as plt
        from matplotlib.figure import Figure
        tree = ast.parse(inspect.getsource(p.render_main))
        statements = [n for n in ast.walk(tree) if isinstance(n, ast.Assign) and any(isinstance(t, ast.Subscript) and isinstance(t.value, ast.Name) and t.value.id == "palettes" for t in n.targets)]
        self.assertEqual(len(statements), 1)
        code = compile(ast.fix_missing_locations(ast.Module(body=statements, type_ignores=[])), "palette", "exec")
        context = dict(palettes={}, column="biological_label", levels=["breast cancer", "normal"], plt=plt, matplotlib=matplotlib, hashlib=hashlib)
        exec(code, context); palette = context["palettes"]["biological_label"]
        colors = np.asarray([palette[label] for label in context["levels"]])
        np.testing.assert_array_equal(colors, np.asarray(plt.get_cmap("tab10").colors[:2]))
        self.assertGreater(np.linalg.norm(colors[0] - colors[1]), .5)
        fig = Figure(); axes = fig.subplots(2, 1)
        artists = [ax.scatter([0, 1], [0, 1], c=colors) for ax in axes]
        np.testing.assert_array_equal(artists[0].get_facecolors(), artists[1].get_facecolors())


if __name__ == "__main__": unittest.main()
