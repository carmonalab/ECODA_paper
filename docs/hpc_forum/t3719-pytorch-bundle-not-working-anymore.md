# Pytorch-bundle not working anymore

- Source: https://hpc-community.unige.ch/t/3719

- Created: 2024-11-12T13:37:36.454Z

- Posts: 3

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Denis.Mongin (2024-11-12T13:37:36.490Z)

## Primary informations
Username: mongin
Cluster: baobab
## Description
Since the upgrade, I cannot load PyTorch-bundle/2.1.2-CUDA-12.1.1 , which was working just beofre
## Steps to Reproduce
Before the last updrade, I was doing:
```
ml  GCC/12.3.0  OpenMPI/4.1.5 PyTorch-bundle/2.1.2-CUDA-12.1.1
```
As decsribed by @Adrien.Albert[@Adrien.Albert](https://hpc-community.unige.ch/u/adrien.albert) New software installed: PyTorch-bundle version 2.1.2-CUDA-12.1.1[New software installed: PyTorch-bundle version 2.1.2-CUDA-12.1.1](https://hpc-community.unige.ch/t/new-software-installed-pytorch-bundle-version-2-1-2-cuda-12-1-1/3609).
But since last upgrade, `OpenMPI/4.1.5` is not listed anymore:
```
 module spider openMPI

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  OpenMPI:
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Description:
      The Open MPI Project is an open source MPI-2 implementation.

     Versions:
        OpenMPI/1.10.3
        OpenMPI/2.1.1
        OpenMPI/2.1.2
        OpenMPI/3.1.1
        OpenMPI/3.1.3
        OpenMPI/3.1.4
        OpenMPI/4.0.3
        OpenMPI/4.0.5
        OpenMPI/4.1.1
        OpenMPI/4.1.4
        OpenMPI/4.1.6

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  For detailed information about a specific "OpenMPI" package (including how to load the modules) use the module's full name.
  Note that names that have a trailing (E) are extensions provided by other modules.
  For example:

     $ module spider OpenMPI/4.1.6
-----------------------------------------------------------------------------------
```
If I do
```
ml  GCC/13.2.0  OpenMPI/4.1.6 PyTorch-bundle/2.1.2-CUDA-12.1.1
```
I instead have
```
Lmod has detected the following error:  The following module(s) are unknown: "OpenMPI/4.1.5" "PyTorch-bundle/2.1.2-CUDA-12.1.1"

Please check the spelling or version number. Also try "module spider ..."
It is also possible your cache file is out-of-date; it may help to try:
  $ module --ignore_cache load "OpenMPI/4.1.5" "PyTorch-bundle/2.1.2-CUDA-12.1.1"

Also make sure that all modulefiles written in TCL start with the string #%Module
```
Would it be possible to fix so that the bundle can be used with OpenMPI/4.1.6 ?


## Post 2 by @Yann.Sagon (2024-11-13T08:30:15.683Z)

Dear @Denis.Mongin[@Denis.Mongin](https://hpc-community.unige.ch/u/denis.mongin)
Indeed the module
```
OpenMPI/4.1.5
```
has disappeared from our installation. We often have to rebuild some libs/soft and we probably had a problem once in the process, thanks for letting us know.
We have rebuilt the module and
```
PyTorch bundle
```
works again.
For the completeness of my answer: this has nothing to do with the maintenance of our cluster, the soft are installed/upgraded outside the maintenance period.
Best


## Post 3 by @Denis.Mongin (2024-11-13T08:34:09.886Z)

Perfect, thank you so much Yann.
