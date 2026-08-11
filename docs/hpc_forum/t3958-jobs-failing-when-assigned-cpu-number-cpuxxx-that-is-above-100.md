# Jobs failing when assigned cpu number (cpuXXX) that is above 100

- Source: https://hpc-community.unige.ch/t/3958

- Created: 2025-05-21T13:37:48.215Z

- Tags: baobab, yggdrasil

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Olivier.Kirchhoffer (2025-05-21T13:37:48.268Z)

## Primary informations
Username: Olivier.Kirchhoffer
Cluster: Yggdrasil
## Description
When I try to launch a job (calling for Gaussian) on a terminal window, I realized that it very often crashes and I have to repeat the launch, until I obtain an attributed cpu with a number below 100 (e.g. cpu097). It seems that whenever the job is launched with a cpu above 100 (e.g. cpu151) it crashes. The .out file in the case of a crash systematically reports the following:
Error: illegal instruction, illegal opcode
rax 0000000000da6150, rbx 00007ffed6211070, rcx 00007ffed6211070
rdx 00007ffed6211040, rsp 00007ffed6210f98, rbp 00007ffed6210fa0
rsi 00007ffed6211048, rdi 00007ffed6210fd8, r8  0000000000000060
r9  0000000000d8c6f0, r10 000000000000002a, r11 0000000000000031
r12 0000000001030930, r13 000014cf4de40830, r14 000014cf4de10cb0
r15 000014cf4de076b0
/lib64/libc.so.6(+0x3e730) [0x14d05ba3e730]
/opt/ebsofts/Gaussian/g16/l101.exe() [0xda61fa]
srun: error: cpu156: task 0: Exited with exit code 1
For now I circumvented this problem by re-launching the job until I obtained a cpu node with a number below 100, but I’m wondering if there might by a typo somewhere in my code or other code that would make it expect a two-number cpu ID or a ‘0’ as first cpu number for it to function?
## Steps to Reproduce
Re-run the following .job file:
---
#!/bin/bash
#SBATCH --job-name=A01
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=4:00:00
#SBATCH --mem=10000mb
#SBATCH --partition=shared-cpu
module load gaussian/g16
PREFIX=A01
JOBINP=${PREFIX}.com
JOBOUT=${PREFIX}.log
srun g16  < ${JOBINP} > ${JOBOUT}
---
with an appropriate .com file in the same folder:
---
%chk=A01.chk
%nprocshared=16
%mem=5GB
#b3lyp/6-31G(d,p) Opt freq SCRF=(Solvent=Methanol)
A01: Geom Optimization
0 1
C     -4.708002   -0.431985   -0.270209
C     -4.333158   -1.733221    0.010764
C     -2.988958   -2.089250    0.048546
C     -2.026146   -1.093887   -0.198955
C     -2.376035    0.219826   -0.543554
C     -3.764224    0.543950   -0.622742
N     -0.644079   -1.166018   -0.167282
C     -0.110948    0.046540   -0.510422
C     -1.138508    0.934752   -0.774126
C      1.246458    0.370048   -0.552097
N      1.600134    1.628411   -0.876301
C      0.613239    2.513221   -1.164292
C     -0.758256    2.241761   -1.127353
O     -4.178583    1.798199   -1.041815
C     -5.414459    1.775211   -1.771229
C      2.196720   -0.594369    0.057288
C      1.733514   -2.044658   -0.134203
C      2.115156   -2.391148   -1.598095
C      3.330810   -1.495932   -1.854441
C      3.601437   -0.900052   -0.464713
C      3.982143   -2.036293    0.531900
C      4.481736   -1.469261    1.874095
C      5.046787   -3.026988    0.044233
O      2.703979   -2.747508    0.721217
C      0.221955   -2.319384    0.170333
O     -0.007614   -2.641070    1.538868
H      2.223626   -0.329380    1.120484
H     -5.776076   -0.197962   -0.252391
H     -5.103609   -2.482006    0.192439
H     -2.699581   -3.109172    0.283860
H      0.967556    3.507784   -1.420354
H     -1.486662    3.013535   -1.342244
H     -5.326161    2.448471   -2.626321
H     -5.666925    0.790864   -2.185487
H     -6.234067    2.120616   -1.130088
H      2.389288   -3.449828   -1.693913
H      1.323200   -2.171791   -2.319477
H      4.180413   -2.049091   -2.268709
H      3.080414   -0.711214   -2.576229
H      4.299640   -0.060819   -0.489150
H      5.412563   -0.909206    1.733636
H      3.762499   -0.806207    2.360809
H      4.672970   -2.275934    2.591206
H      5.977372   -2.512049   -0.217821
H      5.276482   -3.758929    0.828608
H      4.714200   -3.620275   -0.811674
H     -0.140751   -3.155830   -0.437009
H      0.694798   -2.209580    2.053837
---
## Expected Result
It should run and calculate rotational energies and optimize them for a molecule according to its carthesian coordinates.
## Actual Result
It will crash whenever a cpu with a number above 100 is attributed for the job.


## Post 2 by @Adrien.Albert (2025-05-23T09:32:30.755Z)

[HPC][baobab] Gaussian error[[HPC][baobab] Gaussian error](https://hpc-community.unige.ch/t/hpc-baobab-gaussian-error/3062/3) HPC issues[HPC issues](https://hpc-community.unige.ch/c/hpc-support/hpc-issues/14)
> @Pedro.FebrerMartinez[@Pedro.FebrerMartinez](https://hpc-community.unige.ch/u/pedro.febrermartinez) Can you try on an AMD compute node (V8) the latest version we just installed? New software installed: Gaussian version 16.C.02-AVX2[New software installed: Gaussian version 16.C.02-AVX2](https://hpc-community.unige.ch/t/new-software-installed-gaussian-version-16-c-02-avx2/3824)
Could you try again with :
```
(yggdrasil)-[alberta@login1 ~]$ ml spider Gaussian

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  Gaussian: Gaussian/16.C.02-AVX2
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Description:
      Gaussian provides state-of-the-art capabilities for electronic structure modeling. Gaussian 09 is licensed for a wide variety of computer systems. All versions of Gaussian 09 contain every
      scientific/modeling feature, and none imposes any artificial limitations on calculations other than your computing resources and patience. This is the official gaussian AVX2 build. 

     Other possible modules matches:
        gaussian

    This module can be loaded directly: module load Gaussian/16.C.02-AVX2

    Help:
      Description
      ===========
      Gaussian provides state-of-the-art capabilities for electronic structure
      modeling. Gaussian 09 is licensed for a wide variety of computer
      systems. All versions of Gaussian 09 contain every scientific/modeling
      feature, and none imposes any artificial limitations on calculations
      other than your computing resources and patience.
      
      This is the official gaussian AVX2 build.
      
      
      More information
      ================
       - Homepage: https://www.gaussian.com/
```
---
To find other possible module matches execute:
```
  $ module -r spider '.*Gaussian.*'
```
