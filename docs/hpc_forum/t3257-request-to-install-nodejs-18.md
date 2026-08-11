# Request to install nodejs 18

- Source: https://hpc-community.unige.ch/t/3257

- Created: 2024-01-19T09:29:54.881Z

- Posts: 4

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Giuseppe.Chindemi (2024-01-19T09:29:54.944Z)

## Primary informations
Username: chindemi
Cluster: Yggdrasil
## Description
Hello,
Recent versions of jupyterlab require nodejs >= 18.0.0 to install some common extensions:
```
ValueError: Please install nodejs >=18.0.0 before continuing. nodejs may be installed using conda or directly from the nodejs website.
```
Could you please install it?
For reference, I am currently loading these modules:
```
GCCcore/10.2.0 Tkinter/3.8.6 Python/3.8.6 cuDNN/8.6.0.163-CUDA-11.8.0 git-lfs/3.1.2 FFmpeg/4.3.1 nodejs/12.19.0
```
Ideally, I would prefer to simply swap nodejs12 for nodejs18, but if the compiler is too old I could update the whole environment… Python 3.8 is losing support anyway.
Thank you very much!


## Post 2 by @Adrien.Albert (2024-01-19T14:06:50.590Z)

Hi @Giuseppe.Chindemi[@Giuseppe.Chindemi](https://hpc-community.unige.ch/u/giuseppe.chindemi)
Build In progress :wink:


## Post 3 by @Adrien.Albert (2024-01-19T15:15:49.254Z)

New software installed: nodejs version 20.9.0[New software installed: nodejs version 20.9.0](https://hpc-community.unige.ch/t/new-software-installed-nodejs-version-20-9-0/3259) HPC ChangeLog[HPC ChangeLog](https://hpc-community.unige.ch/c/hpc-announce/hpc-changelog/9)
> Dear users, we have installed a new software: nodejs 20.9.0: 

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  nodejs: nodejs/20.9.0
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------…


## Post 4 by @Giuseppe.Chindemi (2024-01-22T10:58:00.396Z)

Perfect, thank you Adrien!
