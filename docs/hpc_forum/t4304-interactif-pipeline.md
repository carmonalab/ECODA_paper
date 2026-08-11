# Interactif pipeline

- Source: https://hpc-community.unige.ch/t/4304

- Created: 2026-06-09T11:06:38.674Z

- Posts: 2

- Category: 13

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Camille.Christe (2026-06-09T11:06:38.732Z)

Dear HPC team,
For a master project, we have build modules with the workflow language Nextflow.
Now I would like to link these modules with a script that ask at the end of each module if the quality of the data are good or some samples have to be remove than the module rerun before going to the next one.
Depending of the number of sample, each module can take hours before the script will ask to review the results. Each nextflow module take care of memory, cpu and time to be allocated at each sub-steps. It runs well via sbatch for the moment.
I would like to know how I can launch the new pipeline that integrate the different modules with script that allow interactivity? Is it the public-interactive-cpu that I have to use? But it has 8hours time and 10GB memory limits.
Thank you very much.
Camille


## Post 2 by @Yann.Sagon (2026-06-10T12:22:34.066Z)

Dear @Camille.Christe[@Camille.Christe](https://hpc-community.unige.ch/u/camille.christe)
The purpose of the partition public-interactive-cpu is to be sure to have a small amount of resources always available: for example if you develop your code, you don’t want to wait until 3AM to have the resource allocated.
If I understand well, In your case, your workflow isn’t interactive as the interaction is done by a script? If yes, you can launch this script as a normal sbatch script?
Maybe this would help you? Slurm Workload Manager - sbatch[Slurm Workload Manager - sbatch](https://slurm.schedmd.com/sbatch.html#OPT_dependency)
If possible, can you share a simple diagram of what you want to achieve, it will be easier to help.
