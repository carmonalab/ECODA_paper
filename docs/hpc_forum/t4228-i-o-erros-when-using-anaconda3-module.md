# I/O erros when using Anaconda3 module

- Source: https://hpc-community.unige.ch/t/4228

- Created: 2026-02-23T07:05:06.935Z

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Max.Briel (2026-02-23T07:05:06.997Z)

## Primary information
Username: briel
Cluster: yggdrasil
## Description
I’m having issues creating/loading conda environments, when using the provided modules: module load `Anaconda3/2024.02-1` and `Anaconda3/2020.07`
## Steps to Reproduce
```
module load Anaconda3/2024.02-1
conda --version
```
## Expected Result
Show the conda version.
## Actual Result
```
Fatal Python error: init_import_site: Failed to import the site module
Python runtime state: initialized
Traceback (most recent call last):
  File "<frozen importlib._bootstrap>", line 1176, in _find_and_load
  File "<frozen importlib._bootstrap>", line 1147, in _find_and_load_unlocked
  File "<frozen importlib._bootstrap>", line 690, in _load_unlocked
  File "<frozen importlib._bootstrap>", line 980, in exec_module
  File "<frozen site>", line 616, in <module>
  File "<frozen site>", line 602, in main
  File "<frozen site>", line 343, in addusersitepackages
  File "<frozen site>", line 226, in addsitedir
  File "<frozen site>", line 179, in addpackage
OSError: [Errno 121] Remote I/O error
```
I’m receiving these I/O related error. For different Anaconda3 version, the error is slightly different, but still I/O related.


## Post 2 by @Adrien.Albert (2026-02-23T11:01:01.109Z)

Hello @Max.Briel[@Max.Briel](https://hpc-community.unige.ch/u/max.briel)
Thank you for your report, here the related post:
[2026] Current issues on HPC Cluster[[2026] Current issues on HPC Cluster](https://hpc-community.unige.ch/t/2026-current-issues-on-hpc-cluster/4185/11) HPC ChangeLog[HPC ChangeLog](https://hpc-community.unige.ch/c/hpc-announce/hpc-changelog/9)
> Dear Users, 
Cluster: Yggdrasil 
Description
An issue is currently affecting the Home storage on the Yggdrasil cluster. 
Some users may encounter the following message: 
“I/O remote error” 

Possible disruptions accessing the Home directory
Risk of input/output errors during commands or job execution
The cluster remains available, but some operations may fail

We are working to restore normal service as soon as possible. 
An update will be posted as soon as we have more information. 
We apologiz…
