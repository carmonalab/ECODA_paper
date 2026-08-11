# FAILED, ExitCode 2 R example

- Source: https://hpc-community.unige.ch/t/3567

- Created: 2024-07-24T13:13:28.844Z

- Posts: 3

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Joel.Terschuur (2024-07-24T13:13:28.889Z)

## Primary informations
Username: terschuu
Cluster: baobab
## Description
I just started using baobab and I intend to run R simulations in parallel with it. To get started I am trying to replicate the example here: High Performance Computing - Baobab Hello World – Hugo Zzo Theme[High Performance Computing - Baobab Hello World – Hugo Zzo Theme](https://blog-dal.netlify.app/posts/baobab_hello_world/). I am trying to run it in the debug-cpu partition. But I get an email saying it failed with ExitCode1.
## Steps to Reproduce
As in the example I have the following R script
```
#Load libraries
library(foreach)

#Simulate population
mysd = 15
mymean = 50
pop = rnorm(10e5, mean = mymean, sd = mysd)

#Define nbr of iterations
B = 10e3
samplesize = 100
myresults = foreach(b = icount(B), .combine = rbind)%do%{
  mysample = sample(pop, size = samplesize)
  mean(mysample)
}

#Theoretical xbar variance
mysd^2 / samplesize

#Observed xbar variance
var(myresults[,1])

#Save results
save(myresults, "xbar_simulation_results.rda")
```
and the following .sh file (I have called it tutorial1.sh)
```
#!/bin/bash

#SBATCH --job-name=simu_R
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0:15:0
#SBATCH --partition=debug-cpu
#SBATCH --mail-type=ALL
#SBATCH --mail-user=joel.terschuur@unige.ch

module load GCC/8.2.0-2.31.1 OpenMPI/3.1.3 R/3.6.0

INFILE=simu.R

OUTFILE=report_simu.Rout

srun R CMD BATCH $INFILE $OUTFILE
```
Then I run the command sbatch tutorial1.sh in the putty console, I am using winSCP to do the file transfer, the files are located at /home/users/t/terschuu/, I have tried to put the whole directory url in infile and outfile and on the sbatch command, but also got error. The last try was Slurm Job_id=11531419.


## Post 2 by @Joel.Terschuur (2024-07-24T13:25:47.561Z)

I added
library(iterators)
and changed
save(myresults, “xbar_simulation_results.rda”)
to
save(myresults, file = “xbar_simulation_results.rda”)
and it seems to work


## Post 3 by @Adrien.Albert (2024-07-24T14:29:24.050Z)

Hi @Joel.Terschuur[@Joel.Terschuur](https://hpc-community.unige.ch/u/joel.terschuur)
Joel.Terschuur:
> and it seems to work
Good !
If you need help again, please provide the output in: report_simu.Rout (or any ouput you configure)
By the way you are using an old version of R:
Here all version available on clusters
```
(baobab)-[alberta@login1 R]$ ml spider R

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  R:
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Description:
      R is a free software environment for statistical computing and graphics.

     Versions:
        R/3.2.3
        R/3.5.1
        R/3.6.0
        R/3.6.2
        R/4.0.0
        R/4.0.4
        R/4.1.0
        R/4.1.2
        R/4.2.0
        R/4.2.1
        R/4.3.2
        R/4.3.3
     Other possible modules matches:
        APR  APR-util  AmberTools  Archive-Zip  Armadillo  Arrow  BioPerl  Bismark  Blender  Bracken  Brotli  Brunsli  CUDAcore  Cartopy  CellProfiler  CellRanger  CellRanger-ARC  CellRanger-ATAC  CmdStanR  ...

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  To find other possible module matches execute:

      $ module -r spider '.*R.*'

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  For detailed information about a specific "R" package (including how to load the modules) use the module's full name.
  Note that names that have a trailing (E) are extensions provided by other modules.
  For example:

     $ module spider R/4.3.3
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
```
