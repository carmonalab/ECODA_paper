# Installing terra and mgcv for biomod2 in R on Baobab

- Source: https://hpc-community.unige.ch/t/4186

- Created: 2026-01-06T17:17:00.335Z

- Tags: baobab

- Posts: 3

- Category: 8

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Sandrine.Ranger (2026-01-06T17:17:00.388Z)

Hi,
I’m using biomod2 package v4.2.6.2 in R on Baobab, and my script works except for GAM models. I need mgcv package, but with R 4.3.2, mgcv won’t install because R is too old (I need R 4.4.0 or more).
I tried: module load GCC/12.3.0
module load GCCcore/12.3.0 module load UDUNITS/2.2.28 module load OpenMPI/4.1.5 module load GDAL/3.7.1 module load PROJ/9.2.0
module load R/4.3.2
and module load GCC/13.3.0 module load R/4.4.2
With R 4.4.2: terra package (required for biomod2) won’t install.
Does anyone know a recent R version and module setup on Baobab that allows installing terra package?
Thank you for your help :slight_smile:


## Post 2 by @pablo.strasser1 (2026-01-06T21:03:00.163Z)

Hi, I don’t use R. However, a generic method to install software not available on the cluster is by using apptainer (previously singularity). There are official Docker images for R available. A singularity image can be build from a docker image. If needed, the docker image can be customized using a Dockerfile.
In summary apptainer allows to create your own environment installing all the software you want. By default your home directory will be mounted.


## Post 3 by @Yann.Sagon (2026-01-07T10:23:33.610Z)

Dear @Sandrine.Ranger[@Sandrine.Ranger](https://hpc-community.unige.ch/u/sandrine.ranger)
The problem is probably due to the fact that you are changing R versions: the packages you installed through R are in the `Rpackages` directory in your home directory.
If you install from R/4.3.2 and then switch to R/4.4.2, there will be conflicts and you will need to recompile the packages. The cleanest solution is to delete all the contents of `Rpackages` and start again.
I see that you have a lot of R packages in Rpackages: many of them are available in the `R-bundle-CRAN/2024.11` module. I suggest you load it, as you will have much less to compile.
If compilation is too slow, I suggest you compile on a computing node.
doc.eresearch.unige.ch[doc.eresearch.unige.ch](https://doc.eresearch.unige.ch/hpc/slurm#interactive_jobs)
### hpc:slurm [eResearch Doc][hpc:slurm [eResearch Doc]](https://doc.eresearch.unige.ch/hpc/slurm#interactive_jobs)
