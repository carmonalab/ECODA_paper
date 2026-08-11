# Problem to install packages "rgdal"

- Source: https://hpc-community.unige.ch/t/3875

- Created: 2025-03-20T09:16:18.745Z

- Posts: 5

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Cecilia.Barouillet (2025-03-20T09:16:18.780Z)

## Primary informations
Username: $barouill
Cluster: $baobab
## Description
I am using the R application for my work and I was able to install all the packages that I needed except for the “terra” which is essential for my work.
In order to install the “terra” package, I need to install the “rgdal” package; however, I was not able to install the latest. I tried to install it using several versions of R but the issue is still persisting.
## Steps to Reproduce
Lines of codes used in the PuTTY terminal:
```
module load GCC/13.2.0 R/4.4.1 
R
install.packages("rgdal") # selected CRAN mirror 63
```
### error message :
```
Warning message package "rgdal" is not available for this version of R. 
#when trying to install it from source
url<-"https://cran.r-project.org/src/contrib/Archive/rgdal/rgdal_1.6-7.tar.gz"
install.packages(url, type="source", repos=NULL)
installation of packages ‘/tmp/RtmphEagt0/downloaded_packages/rgdal_1.6-7.tar.gz’ had non-zero exit status.
```
when trying another version of R
```
GCC/8.3.0  OpenMPI/3.1.4 R/3.6.2
```
reiterated the lines above and I had the same issues.


## Post 2 by @Simon.Hug (2025-04-22T14:23:14.638Z)

Dear Cecilia Barrouillet
As I encountered similar problems, here a few hints:
i) The rgdal package has been removed from the cran repository, because the author retired and no one is continuing the development of this package (that’s why you can not install it directly from a repository).
ii) While your approach (i.e., installing from source an archived version of the package) is the correct one, you were a bit too pessimistic regarding how far back you had to go for R versions. The last archived version of rgdal is of 2023, and if my memory is corrected it required at least R 4.0.
iii) In the past (it does no longer work), I was able to use rgdal on HPC by loading the following modules:
module load GCC/9.3.0  OpenMPI/4.0.3 R/4.0.0
module load rgdal/1.4-8-R-4.0.0
While these module still load, R is unable to load rgdal.
iv) According to
module spider rgdal
using
module load GCC/10.2.0  OpenMPI/4.0.5
module load rgdal/1.5-23-R-4.0.4
might be your best shot (as I have currently other HPC problems, I can not try this out).
best wishes


## Post 3 by @Yann.Sagon (2025-04-24T08:08:30.359Z)

Dear @Simon.Hug[@Simon.Hug](https://hpc-community.unige.ch/u/simon.hug) thanks for the post!
Simon.Hug:
> as I have currently other HPC problems, I can not try this out
Feel free to ask for help in a dedicated post.
Best


## Post 4 by @Simon.Hug (2025-04-24T12:59:50.869Z)

Luckily enough you were able to solve my other issues (other post). However, what I suggested as solution for the rgdal package does not work (probably due to what you mentioned in your other post: old R version). Just to keep a trace on this:
While
module spider rgdal
suggests using
module load GCC/10.2.0 OpenMPI/4.0.5
module load rgdal/1.5-23-R-4.0.4
when using R and trying to install dependency packages for rgdal, I get the same error message:
Warning: unable to access index for repository https://stat.ethz.ch/CRAN/src/contrib:[https://stat.ethz.ch/CRAN/src/contrib:](https://stat.ethz.ch/CRAN/src/contrib:)
internet routines cannot be loaded
Warning messages:
1: In download.file(url, destfile = f, quiet = TRUE) :
unable to load shared object ‘/opt/ebsofts/R/4.0.4-foss-2020b/lib64/R/modules//internet.so’:
libssl.so.10: cannot open shared object file: No such file or directory
2: package ‘sp’ is not available for this version of R
best wishes


## Post 5 by @Simon.Hug (2025-04-25T07:17:50.328Z)

Dear all
So there is a way to still use the (retired) rgdal package in R:
module load  GCC/11.3.0  OpenMPI/4.1.4
module load rgdal/1.6-6
module load  R/4.2.1
Then you can install from source (rgdal_1.6-6.tar.gz) the rgdal package in R.
As this is not documented in
module spider rgdal
It might be useful to refer to it (i.e., how to load the rgdal module with the corresponding/compatible R module).
best wishes
