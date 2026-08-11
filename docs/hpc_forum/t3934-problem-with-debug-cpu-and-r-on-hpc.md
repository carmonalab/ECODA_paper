# Problem with debug-cpu and R on HPC

- Source: https://hpc-community.unige.ch/t/3934

- Created: 2025-04-22T12:56:26.462Z

- Tags: baobab, bamboo, yggdrasil

- Posts: 4

- Category: 13

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Simon.Hug (2025-04-22T12:56:26.510Z)

Dear all
I encounter some strange behavior when using the debug-cpu nodes on bamboo, baobab, yggdrasil when using R. On these nodes I load the following modules (to run JAGS):
```
module load GCC/9.3.0  OpenMPI/4.0.3 R/4.0.0
module load V8/3.4.0-R-4.0.0
```
A first strange behavior is that when running R, some packages that are installed (amongst them runjags and rjags: when submitting jobs with the exact same setup to public-cpu, R finds these packages) appear not to be installed.
When I then try to install the missing packages (e.g. runjags)
```
install.packages("runjags")
```
and selection the Swiss repository I get the following warnings (please note that I get the same last warning (“not available (for R version 4.0.0)”) for each and every package I wish to install:
```
Warning: unable to access index for repository https://stat.ethz.ch/CRAN/src/contrib:
  internet routines cannot be loaded
Warning messages:
1: In download.file(url, destfile = f, quiet = TRUE) :
  unable to load shared object '/opt/ebsofts/R/4.0.0-foss-2020a/lib64/R/modules//internet.so':
  libssl.so.10: cannot open shared object file: No such file or directory
2: package ‘runjags’ is not available (for R version 4.0.0)
```
I have proceeded similarly in the past, and did not encounter any of these problems. Thanks for any help and best wishes, simon
P.S. This problem appears also related to my query about the V8 module. If on debug-cpu I load the new V8 module as instructed I’m unable to start R:
```
(baobab)-[hugs@cpu001 ~]$ module load GCC/13.2.0
(baobab)-[hugs@cpu001 ~]$ module load V8/6.0.1-R-4.4.1
(baobab)-[hugs@cpu001 ~]$ R --no-restore
Illegal instruction (core dumped)
```


## Post 2 by @Yann.Sagon (2025-04-24T09:19:50.024Z)

Dear @Simon.Hug[@Simon.Hug](https://hpc-community.unige.ch/u/simon.hug)
Simon.Hug:
> This problem appears also related to my query about the V8 module. If on debug-cpu I load the new V8 module as instructed I’m unable to start R
I’ve checked as well and I confirm it doesn’t work on cpu001 on Baobab. The reason is that openblas is compiled for “haswell” architecture and this server is a virtual machine that doesn’t supports all the needed stuff. I’ve tried on another host and it works. We’ll see if we can fix the issue for cpu001. You can use the debug partition yggdrasil and bamboo, it should work.
Simon.Hug:
> `unable to load shared object '/opt/ebsofts/R/4.0.0-foss-2020a/lib64/R/modules//internet.so':`
This indicates ssl provided with R is too old. I’ll recompile it,I’ll let you know if it works.
Is it possible to use a more recent version of R? I’m not sure I’ll be able to recompile it.


## Post 3 by @Simon.Hug (2025-04-24T12:50:26.273Z)

Simon.Hug:
> `R --no-restore`
Thanks, the issue with the new V8 module on the debug-cpu nodes has disappeared (even on baobab).
Regarding the other issue: I have been able to resolve the issue by using modules that are compatible with R/4.1.1 (from “module spider JAGS” it was unclear to me whether this might also work with more recent versions of R).
Thanks for your help and best wishes


## Post 4 by @Yann.Sagon (2025-04-25T09:10:42.091Z)

Simon.Hug:
> even on baobab
My colleague @Gael.Rossignol[@Gael.Rossignol](https://hpc-community.unige.ch/u/Gael.Rossignol) fixed that, good to hear it is working!
I have as well compiled a newer rjags version for R 4.3.2, I’ll publish the info later.
