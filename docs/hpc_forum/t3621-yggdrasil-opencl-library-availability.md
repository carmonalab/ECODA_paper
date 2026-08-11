# Yggdrasil openCL library availability

- Source: https://hpc-community.unige.ch/t/3621

- Created: 2024-09-02T09:33:32.215Z

- Tags: yggdrasil

- Posts: 3

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Gregory.Jevardat (2024-09-02T09:33:32.283Z)

Hello
I cannot find where are located openCL libs. Module spider cannot find them.
I Need them to build a project relying on them ! ( TornadoVM project)
Where can I find them?
Thanks for your help
Gregory


## Post 2 by @Adrien.Albert (2024-09-24T09:46:22.628Z)

Hi @Gregory.Jevardat[@Gregory.Jevardat](https://hpc-community.unige.ch/u/gregory.jevardat)
Are you talking about this : PyOpenCL[PyOpenCL](https://mathema.tician.de/software/pyopencl/)


## Post 3 by @Yann.Sagon (2024-11-19T07:47:39.297Z)

Or this one?
```
(baobab)-[sagon@login1 ~]$ ml spider pocl/4.0-CUDA-12.1.1

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  pocl: pocl/4.0-CUDA-12.1.1
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Description:
      PoCL is a portable open source (MIT-licensed) implementation of the OpenCL standard (1.2 with some 2.0 features supported).

    You will need to load all module(s) on any one of the lines below before the "pocl/4.0-CUDA-12.1.1" module is available to load.

      GCC/12.3.0

    Help:
      Description
      ===========
      PoCL is a portable open source (MIT-licensed) implementation
      of the OpenCL standard (1.2 with some 2.0 features supported).

      More information
      ================
       - Homepage: http://portablecl.org
```
Then, once loaded, you can see where the libs are located:
```
(baobab)-[sagon@login1 ~]$ pkg-config --libs pocl
-L/opt/ebsofts/pocl/4.0-GCC-12.3.0-CUDA-12.1.1/lib64 -lpocl
```
