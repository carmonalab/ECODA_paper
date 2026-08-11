# Virtualenv command not found despite loading python

- Source: https://hpc-community.unige.ch/t/3281

- Created: 2024-01-31T15:01:57.606Z

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Abigail.Licata (2024-01-31T15:01:57.665Z)

## Primary informations
Username: $licata
Cluster: $BAOBAB
## Description
I load the relevant modules for Python 3.11.5 (GCCcore-13.2.0, Python 3.11.5) and when I list my loaded modules, they are both there, although I cannot create any virtual environments. I get the error: -bash: virtualenv: command not found
## Steps to Reproduce
module load GCCcore/13.2.0
module load Python/3.11.5
module list
which python
virtualenv ~/baobab_mne
## Expected Result
Normally, it creates the empty virtual environment. Recently, it isn’t recognizing the commands virtualenv or venv
## Actual Result
-bash: venv: command not found


## Post 2 by @Yann.Sagon (2024-01-31T15:27:34.005Z)

Dear @Abigail.Licata[@Abigail.Licata](https://hpc-community.unige.ch/u/abigail.licata) we discovered that `virtualenv` isn’t bundled with recent Python module anymore. You need to load `virtualenv` module in addition to Python module. We’ve installed it right now: New software installed: virtualenv version 20.24.6[New software installed: virtualenv version 20.24.6](https://hpc-community.unige.ch/t/new-software-installed-virtualenv-version-20-24-6/3282)
