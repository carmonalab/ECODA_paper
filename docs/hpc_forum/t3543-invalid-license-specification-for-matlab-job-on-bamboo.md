# Invalid license specification for matlab job on bamboo

- Source: https://hpc-community.unige.ch/t/3543

- Created: 2024-07-16T13:28:42.542Z

- Tags: bamboo

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Nicolas.Cuny (2024-07-16T13:28:42.579Z)

Username: cuny
Cluster: bamboo
What I tried:
I’m trying to run matlab jobs on the new bamboo cluster with command
sbatch ./run.sh
with run.sh executable file of this type:
#!/bin/bash
#SBATCH --partition=shared-cpu
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --licenses=matlab@matlablm.unige.ch
unset DISPLAY
module load MATLAB/2022a
matlab -nodesktop -nosplash -nodisplay -r  “some matlab script”
What I get:
I get the following error message: “sbatch: error: Batch job submission failed: Invalid license specification”


## Post 2 by @Adrien.Albert (2024-07-16T16:52:57.473Z)

Hi @Nicolas.Cuny[@Nicolas.Cuny](https://hpc-community.unige.ch/u/nicolas.cuny)
The license has been added to Bamboo:
```
(bamboo)-[root@login1 ~]$ scontrol show lic
LicenseName=atlas_beegfs@slurmdb
    Total=42949672 Used=0 Free=42949672 Reserved=0 Remote=yes
    LastConsumed=0 LastDeficit=0 LastUpdate=2024-07-16T18:37:46
LicenseName=compiler@matlablm.unige.ch
    Total=1 Used=0 Free=1 Reserved=0 Remote=yes
    LastConsumed=0 LastDeficit=0 LastUpdate=2024-07-16T18:51:33
LicenseName=control_toolbox@matlablm.unige.ch
    Total=0 Used=0 Free=0 Reserved=0 Remote=yes
    LastConsumed=0 LastDeficit=0 LastUpdate=2024-07-16T18:51:39
LicenseName=curve_fitting_toolbox@matlablm.unige.ch
    Total=1 Used=0 Free=1 Reserved=0 Remote=yes
    LastConsumed=0 LastDeficit=0 LastUpdate=2024-07-16T18:51:13
LicenseName=distrib_computing_toolbox@matlablm.unige.ch
    Total=9 Used=0 Free=9 Reserved=0 Remote=yes
    LastConsumed=0 LastDeficit=0 LastUpdate=2024-07-16T18:49:25
LicenseName=image_toolbox@matlablm.unige.ch
    Total=8 Used=0 Free=8 Reserved=0 Remote=yes
    LastConsumed=0 LastDeficit=0 LastUpdate=2024-07-16T18:49:37
LicenseName=matlab@matlablm.unige.ch
    Total=54 Used=0 Free=54 Reserved=0 Remote=yes
    LastConsumed=0 LastDeficit=0 LastUpdate=2024-07-16T18:46:47
LicenseName=optimization-toolbox@matlablm.unige.ch
    Total=5 Used=0 Free=5 Reserved=0 Remote=yes
    LastConsumed=0 LastDeficit=0 LastUpdate=2024-07-16T18:49:43
LicenseName=reinforcement_learn_toolbox@matlablm.unige.ch
    Total=0 Used=0 Free=0 Reserved=0 Remote=yes
    LastConsumed=0 LastDeficit=0 LastUpdate=2024-07-16T18:51:45
LicenseName=signal_toolbox@matlablm.unige.ch
    Total=4 Used=0 Free=4 Reserved=0 Remote=yes
    LastConsumed=0 LastDeficit=0 LastUpdate=2024-07-16T18:51:29
LicenseName=statistics_toolbox@matlablm.unige.ch
    Total=15 Used=0 Free=15 Reserved=0 Remote=yes
    LastConsumed=0 LastDeficit=0 LastUpdate=2024-07-16T18:49:05
LicenseName=wavelet-toolbox@matlablm.unige.ch
    Total=4 Used=0 Free=4 Reserved=0 Remote=yes
    LastConsumed=0 LastDeficit=0 LastUpdate=2024-07-16T18:45:46
```
