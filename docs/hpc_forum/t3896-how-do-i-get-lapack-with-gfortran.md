# How do I get LAPACK with gfortran?

- Source: https://hpc-community.unige.ch/t/3896

- Created: 2025-04-01T04:30:43.922Z

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Christophe.Berthod (2025-04-01T04:30:43.976Z)

Good morning,
foss should provide LAPACK support:
$ module whatis foss
foss/2024a          : Description: GNU Compiler Collection (GCC) based compiler toolchain, including
OpenMPI for MPI support, OpenBLAS (BLAS and LAPACK support), FFTW and ScaLAPACK.
foss/2024a          : Homepage: https://easybuild.readthedocs.io/en/master/Common-toolchains.html#foss-toolchain[https://easybuild.readthedocs.io/en/master/Common-toolchains.html#foss-toolchain](https://easybuild.readthedocs.io/en/master/Common-toolchains.html#foss-toolchain)
foss/2024a          : URL: https://easybuild.readthedocs.io/en/master/Common-toolchains.html#foss-toolchain[https://easybuild.readthedocs.io/en/master/Common-toolchains.html#foss-toolchain](https://easybuild.readthedocs.io/en/master/Common-toolchains.html#foss-toolchain)
However I get:
$ module load foss
$ gfortran test.f90 -llapack
/opt/ebsofts/binutils/2.42-GCCcore-13.3.0/bin/ld: cannot find -llapack: No such file or directory
collect2: error: ld returned 1 exit status


## Post 2 by @Yann.Sagon (2025-04-01T13:07:41.093Z)

Dear @Christophe.Berthod[@Christophe.Berthod](https://hpc-community.unige.ch/u/christophe.berthod)
Christophe.Berthod:
> foss should provide LAPACK support
Yes it does through OpenBLAS. Thus you need to compile against openblas.
Something like `gfortran test.f90 -lopenblas`
