# Help specifying job to run a compiled Fortran program that uses MPI

- Source: https://hpc-community.unige.ch/t/3924

- Created: 2025-04-16T18:40:08.472Z

- Tags: yggdrasil

- Posts: 4

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Genevieve.Savard (2025-04-16T18:40:08.516Z)

Hello,
I am trying to run this software: GitHub - akuhara/SEIS_FILO: SEISmological transdimensional inversion tools for Flat and Isotropic Layered structure in the Ocean[GitHub - akuhara/SEIS_FILO: SEISmological transdimensional inversion tools for Flat and Isotropic Layered structure in the Ocean](https://github.com/akuhara/SEIS_FILO/tree/main)
It is written in Fortran and I have successfully compiled it (I think), using the following modules:
```
Currently Loaded Modules:
  1) GCCcore/13.3.0   4) GCC/13.3.0        7) XZ/5.4.5             10) hwloc/2.10.0     13) UCX/1.16.0        16) PRRTE/3.0.5    19) FFTW/3.3.10
  2) zlib/1.3.1       5) OpenBLAS/0.3.27   8) libxml2/2.12.7       11) OpenSSL/3        14) libfabric/1.21.0  17) UCC/1.3.0
  3) binutils/2.42    6) numactl/2.0.18    9) libpciaccess/0.18.1  12) libevent/2.1.12  15) PMIx/5.0.2        18) OpenMPI/5.0.3
```
Now I am trying to run the examples. For example, one can test the installation by going in subfolder sample/joint_inv and running the script with:
```
mpirun -np 20 ../../bin/joint_inv joint_inv.in
```
I specified by salloc job with
```
salloc  -n1 --ntasks=20 --partition=shared-cpu --time=4:00:00  --mem=16G
```
However I am getting runtime errors “SIGILL: Illegal instruction”:
For example, trying with 4 workers:
```
mpirun -np 4 ../../bin/joint_inv joint_inv.in

Program received signal SIGILL: Illegal instruction.

Backtrace for this error:

Program received signal SIGILL: Illegal instruction.

Backtrace for this error:

Program received signal SIGILL: Illegal instruction.

Backtrace for this error:

Program received signal SIGILL: Illegal instruction.

Backtrace for this error:
#0  0x14860ac3e72f in ???
#1  0x4025c7 in ???
#2  0x40249c in ???
#3  0x14860ac295cf in ???
#4  0x14860ac2967f in ???
#5  0x4024e4 in ???
#6  0xffffffffffffffff in ???
#0  0x146db0e3e72f in ???
#1  0x4025c7 in ???
#2  0x40249c in ???
#3  0x146db0e295cf in ???
#4  0x146db0e2967f in ???
#5  0x4024e4 in ???
#6  0xffffffffffffffff in ???
#0  0x1471c3e3e72f in ???
#1  0x4025c7 in ???
#2  0x40249c in ???
#3  0x1471c3e295cf in ???
#4  0x1471c3e2967f in ???
#5  0x4024e4 in ???
#6  0xffffffffffffffff in ???
#0  0x147093a3e72f in ???
#1  0x4025c7 in ???
#2  0x40249c in ???
#3  0x147093a295cf in ???
#4  0x147093a2967f in ???
#5  0x4024e4 in ???
#6  0xffffffffffffffff in ???
--------------------------------------------------------------------------
prterun noticed that process rank 2 with PID 1617008 on node cpu155 exited on
signal 4 (Illegal instruction).
--------------------------------------------------------------------------
```
I am not sure what I did wrong. It runs on my personal ubuntu laptop without issue.
Thanks in advance for any tips/help.
Kind regards
Genevieve


## Post 2 by @Gael.Rossignol (2025-04-17T13:37:11.560Z)

Dear Genevieve,
On which cluster and node did you compile software?
For example on baobab we can have cpu 10 years older than latest installed and compilations on new nodes cannot work on older nodes.
We need to check you compile on a node compatible with all other nodes.
Best regards,


## Post 3 by @Genevieve.Savard (2025-04-17T16:52:38.336Z)

Hi Gael, thanks for your reply, it is on Yggdrasil, shared-cpu partition.


## Post 4 by @Yann.Sagon (2025-04-23T15:31:48.582Z)

Dear @Genevieve.Savard[@Genevieve.Savard](https://hpc-community.unige.ch/u/genevieve.savard)
I’ve installed seisfilo with a lots of trial and error on the cluster but unfortunately when trying the same example as you did, the issue is similar.
I’ve open an issue on the github of the author here Segmentation fault · Issue #32 · akuhara/SEIS_FILO · GitHub[Segmentation fault · Issue #32 · akuhara/SEIS_FILO · GitHub](https://github.com/akuhara/SEIS_FILO/issues/32)
I don’t have an other idea right now.
Are you sure the sample code is working as expected? You said it is working on your laptop: did you compile it yourself? Did you compile it using mpifort? which version?
Best
Yann
