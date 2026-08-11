# Install.packages failed in Rstudio

- Source: https://hpc-community.unige.ch/t/3398

- Created: 2024-03-27T15:51:55.262Z

- Posts: 5

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Mariana.Bouvier-Rabe (2024-03-27T15:51:55.322Z)

Hello,
I tried to install packages in Rstudio (on Baobab), but Rstudio gives me an error message below.
How could I solve the problem ?
Thank you,
Mariana
## Results
```
install.packages("tidync")
library(tidync)
```
Error: package or namespace load failed for ‘tidync’ in dyn.load(file, DLLpath = DLLpath, …):
unable to load shared object ‘/home/users/b/R/x86_64-pc-linux-gnu-library/4.2/RNetCDF/libs/RNetCDF.so’:
libnetcdf.so.19: cannot open shared object file: No such file or directory


## Post 2 by @Gael.Rossignol (2024-04-03T08:19:45.305Z)

Dear Mariana,
Did you try to load NetCDF software with command :
`ml GCC/12.2.0 OpenMPI/4.1.4 netCDF/4.9.0`
Before to install this software?
I seems you miss the library to compile the software.
Best regards,


## Post 3 by @Mariana.Bouvier-Rabe (2024-04-04T08:13:38.109Z)

Hello,
Thank you for your answer. However, I loaded NetCDF (attached) and it still doesn’t work, same message.
image
Best regards,
Mariana Bouvier Rabe


## Post 4 by @Gael.Rossignol (2024-04-04T09:11:41.878Z)

Hello,
I think you have to follow this documentation in order to install R packages :
https://doc.eresearch.unige.ch/hpc/applications_and_libraries#r_packages[https://doc.eresearch.unige.ch/hpc/applications_and_libraries#r_packages](https://doc.eresearch.unige.ch/hpc/applications_and_libraries#r_packages)
Best regards,


## Post 5 by @Adrien.Albert (2024-04-16T13:44:56.054Z)

Dear @Mariana.Bouvier-Rabe[@Mariana.Bouvier-Rabe](https://hpc-community.unige.ch/u/mariana.bouvier-rabe)
The interactive session for R was implemented by @Julien.Prados[@Julien.Prados](https://hpc-community.unige.ch/u/julien.prados) and he also maintains an R singularity image for his own research group.
He has kindly built a new image by inserting NetCDF. :pray:t3:
You can use the image store here :
`/acanas/m-BioinfoSupport/singularity/ngs_latest.sif`
I tested it, and it’s working for me.
