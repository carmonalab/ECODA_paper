# Python pip install requires sudo privileges when building wheels

- Source: https://hpc-community.unige.ch/t/3696

- Created: 2024-10-23T07:40:40.417Z

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Denis.Mongin (2024-10-23T07:40:40.458Z)

## Primary informations
Username: mongin
Cluster: baobab
## Description
I am trying to install a python library, that requires sudo in the middle of the install. As I do not have sudo, I cannot finish the install. I am not sure how to deal with this situation (problem happens when python builds wheel, `Building wheels for collected packages: pythonmonkey, pminit`)
The why of this package is to be able to call the jsonrepair js library in python.
github.com[github.com](https://github.com/josdejong/jsonrepair)
### GitHub - josdejong/jsonrepair: Repair invalid JSON documents[GitHub - josdejong/jsonrepair: Repair invalid JSON documents](https://github.com/josdejong/jsonrepair)
Repair invalid JSON documents
## Steps to Reproduce
I am doing:
```
ml  GCC/12.3.0  OpenMPI/4.1.5 PyTorch-bundle/2.1.2-CUDA-12.1.1
~/baobab_python_env_LLM3/bin/pip install pythonmonkey
```
The output is:
```
Collecting pythonmonkey
  Using cached pythonmonkey-1.0.0.tar.gz (305 kB)
  Installing build dependencies ... done
  Getting requirements to build wheel ... done
  Preparing metadata (pyproject.toml) ... done
Collecting aiohttp<4.0.0,>=3.9.5 (from aiohttp[speedups]<4.0.0,>=3.9.5->pythonmonkey)
  Using cached aiohttp-3.10.10-cp311-cp311-manylinux_2_17_x86_64.manylinux2014_x86_64.whl.metadata (7.6 kB)
Collecting pminit>=0.4.0 (from pythonmonkey)
  Using cached pminit-1.0.0.tar.gz (7.5 kB)
  Installing build dependencies ... done
  Getting requirements to build wheel ... done
  Preparing metadata (pyproject.toml) ... done
Collecting aiohappyeyeballs>=2.3.0 (from aiohttp<4.0.0,>=3.9.5->aiohttp[speedups]<4.0.0,>=3.9.5->pythonmonkey)
  Using cached aiohappyeyeballs-2.4.3-py3-none-any.whl.metadata (6.1 kB)
Collecting aiosignal>=1.1.2 (from aiohttp<4.0.0,>=3.9.5->aiohttp[speedups]<4.0.0,>=3.9.5->pythonmonkey)
  Using cached aiosignal-1.3.1-py3-none-any.whl.metadata (4.0 kB)
Requirement already satisfied: attrs>=17.3.0 in /opt/ebsofts/Python-bundle-PyPI/2023.06-GCCcore-12.3.0/lib/python3.11/site-packages (from aiohttp<4.0.0,>=3.9.5->aiohttp[speedups]<4.0.0,>=3.9.5->pythonmonkey) (23.1.0)
Collecting frozenlist>=1.1.1 (from aiohttp<4.0.0,>=3.9.5->aiohttp[speedups]<4.0.0,>=3.9.5->pythonmonkey)
  Using cached frozenlist-1.4.1-cp311-cp311-manylinux_2_5_x86_64.manylinux1_x86_64.manylinux_2_17_x86_64.manylinux2014_x86_64.whl.metadata (12 kB)
Collecting multidict<7.0,>=4.5 (from aiohttp<4.0.0,>=3.9.5->aiohttp[speedups]<4.0.0,>=3.9.5->pythonmonkey)
  Using cached multidict-6.1.0-cp311-cp311-manylinux_2_17_x86_64.manylinux2014_x86_64.whl.metadata (5.0 kB)
Collecting yarl<2.0,>=1.12.0 (from aiohttp<4.0.0,>=3.9.5->aiohttp[speedups]<4.0.0,>=3.9.5->pythonmonkey)
  Using cached yarl-1.16.0-cp311-cp311-manylinux_2_17_x86_64.manylinux2014_x86_64.whl.metadata (63 kB)
Collecting Brotli (from aiohttp[speedups]<4.0.0,>=3.9.5->pythonmonkey)
  Using cached Brotli-1.1.0-cp311-cp311-manylinux_2_17_x86_64.manylinux2014_x86_64.whl.metadata (5.5 kB)
Collecting aiodns>=3.2.0 (from aiohttp[speedups]<4.0.0,>=3.9.5->pythonmonkey)
  Using cached aiodns-3.2.0-py3-none-any.whl.metadata (4.0 kB)
Collecting pycares>=4.0.0 (from aiodns>=3.2.0->aiohttp[speedups]<4.0.0,>=3.9.5->pythonmonkey)
  Using cached pycares-4.4.0-cp311-cp311-manylinux_2_17_x86_64.manylinux2014_x86_64.whl.metadata (4.1 kB)
Requirement already satisfied: idna>=2.0 in /opt/ebsofts/Python-bundle-PyPI/2023.06-GCCcore-12.3.0/lib/python3.11/site-packages (from yarl<2.0,>=1.12.0->aiohttp<4.0.0,>=3.9.5->aiohttp[speedups]<4.0.0,>=3.9.5->pythonmonkey) (3.4)
Collecting propcache>=0.2.0 (from yarl<2.0,>=1.12.0->aiohttp<4.0.0,>=3.9.5->aiohttp[speedups]<4.0.0,>=3.9.5->pythonmonkey)
  Using cached propcache-0.2.0-cp311-cp311-manylinux_2_17_x86_64.manylinux2014_x86_64.whl.metadata (7.7 kB)
Requirement already satisfied: cffi>=1.5.0 in /opt/ebsofts/Python-bundle-PyPI/2023.06-GCCcore-12.3.0/lib/python3.11/site-packages (from pycares>=4.0.0->aiodns>=3.2.0->aiohttp[speedups]<4.0.0,>=3.9.5->pythonmonkey) (1.15.1)
Requirement already satisfied: pycparser in /opt/ebsofts/Python-bundle-PyPI/2023.06-GCCcore-12.3.0/lib/python3.11/site-packages (from cffi>=1.5.0->pycares>=4.0.0->aiodns>=3.2.0->aiohttp[speedups]<4.0.0,>=3.9.5->pythonmonkey) (2.21)
Using cached aiohttp-3.10.10-cp311-cp311-manylinux_2_17_x86_64.manylinux2014_x86_64.whl (1.3 MB)
Using cached aiodns-3.2.0-py3-none-any.whl (5.7 kB)
Using cached aiohappyeyeballs-2.4.3-py3-none-any.whl (14 kB)
Using cached aiosignal-1.3.1-py3-none-any.whl (7.6 kB)
Using cached frozenlist-1.4.1-cp311-cp311-manylinux_2_5_x86_64.manylinux1_x86_64.manylinux_2_17_x86_64.manylinux2014_x86_64.whl (272 kB)
Using cached multidict-6.1.0-cp311-cp311-manylinux_2_17_x86_64.manylinux2014_x86_64.whl (129 kB)
Using cached yarl-1.16.0-cp311-cp311-manylinux_2_17_x86_64.manylinux2014_x86_64.whl (342 kB)
Using cached Brotli-1.1.0-cp311-cp311-manylinux_2_17_x86_64.manylinux2014_x86_64.whl (2.9 MB)
Using cached propcache-0.2.0-cp311-cp311-manylinux_2_17_x86_64.manylinux2014_x86_64.whl (236 kB)
Using cached pycares-4.4.0-cp311-cp311-manylinux_2_17_x86_64.manylinux2014_x86_64.whl (288 kB)
Building wheels for collected packages: pythonmonkey, pminit
  Building wheel for pythonmonkey (pyproject.toml) ... -
We trust you have received the usual lecture from the local System
Administrator. It usually boils down to these three things:

    #1) Respect the privacy of others.
    #2) Think before you type.
    #3) With great power comes great responsibility.

[sudo] password for mongin: 
```


## Post 2 by @Yann.Sagon (2024-10-23T12:43:16.222Z)

Denis.Mongin:
```
ml  GCC/12.3.0  OpenMPI/4.1.5 PyTorch-bundle/2.1.2-CUDA-12.1.1
~/baobab_python_env_LLM3/bin/pip install pythonmonkey
```
Hi @Denis.Mongin[@Denis.Mongin](https://hpc-community.unige.ch/u/denis.mongin)
you should try something like
```
pip install --user pythonmonkey
```
