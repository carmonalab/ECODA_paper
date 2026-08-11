# Geant4 based package import failing due to 'GLIBC_2.34' Not found

- Source: https://hpc-community.unige.ch/t/3836

- Created: 2025-02-21T16:00:54.124Z

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Stephen.Mulligan (2025-02-21T16:00:54.167Z)

## Primary informations
Username: mulligas
Cluster: baobab
## Description
I have been using a .sif container containing a custom Gean4 to python integration on baobab. I have had no issues using this until trying to run things today after maintenance was completed. Geant4 related imports now fail due to `ImportError: /lib/x86_64-linux-gnu/libc.so.6: version `GLIBC_2.34’ not found (required by /.singularity.d/libs/libGLX.so.0)`
## Steps to Reproduce
using the container `/srv/beegfs/scratch/groups/rodem/containers/minicalosim.sif` run the following in a python script `from  G4Calo import GeometryDescriptor,`
## Expected Result
GeometryDescriptor should be successfully imported and simulation should be run without issue
## Actual Result
```
Traceback (most recent call last):
  File "/usr/lib/python3.8/threading.py", line 932, in _bootstrap_inner
    self.run()
  File "/usr/lib/python3.8/threading.py", line 870, in run
    self.run()
  File "/usr/lib/python3.8/threading.py", line 870, in run
    self._target(*self._args, **self._kwargs)
  File "/home/users/m/mulligas/calo_working/tl_inner_loop/modules/experiments/reinit_testing_50_events/50_events_1_layer_tl/run_1/tools.py", line 110, in _worker
    self._target(*self._args, **self._kwargs)
  File "/home/users/m/mulligas/calo_working/tl_inner_loop/modules/experiments/reinit_testing_50_events/50_events_1_layer_tl/run_1/tools.py", line 110, in _worker
    result = task(*args, **kwargs)
  File "/home/users/m/mulligas/calo_working/tl_inner_loop/modules/experiments/reinit_testing_50_events/50_events_1_layer_tl/run_1/generator.py", line 372, in generate_single
Exception in thread Thread-3:
Traceback (most recent call last):
  File "/usr/lib/python3.8/threading.py", line 932, in _bootstrap_inner
    result = task(*args, **kwargs)
  File "/home/users/m/mulligas/calo_working/tl_inner_loop/modules/experiments/reinit_testing_50_events/50_events_1_layer_tl/run_1/generator.py", line 372, in generate_single
    self.run()
  File "/usr/lib/python3.8/threading.py", line 870, in run
    self._target(*self._args, **self._kwargs)
  File "/home/users/m/mulligas/calo_working/tl_inner_loop/modules/experiments/reinit_testing_50_events/50_events_1_layer_tl/run_1/tools.py", line 110, in _worker
    from  G4Calo import GeometryDescriptor, run_batch
ImportError: cannot import name 'GeometryDescriptor' from 'G4Calo' (/usr/local/lib/python3.8/dist-packages/G4Calo.py)
    from  G4Calo import GeometryDescriptor, run_batch
  File "/usr/local/lib/python3.8/dist-packages/G4Calo.py", line 14, in <module>
    from minicalo import GeometryDescriptor
ImportError: /lib/x86_64-linux-gnu/libc.so.6: version `GLIBC_2.34' not found (required by /.singularity.d/libs/libGLX.so.0)
    result = task(*args, **kwargs)
  File "/home/users/m/mulligas/calo_working/tl_inner_loop/modules/experiments/reinit_testing_50_events/50_events_1_layer_tl/run_1/generator.py", line 372, in generate_single
```


## Post 2 by @Stephen.Mulligan (2025-02-24T12:03:31.901Z)

fixed by manually linking geant4 libraries with singularity arg
`--env LD_LIBRARY_PATH=/opt/geant4-v11.1.2/lib:/lib/x86_64-linux-gnu:/usr/lib64`
