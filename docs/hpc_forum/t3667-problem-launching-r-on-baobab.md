# Problem launching R on baobab

- Source: https://hpc-community.unige.ch/t/3667

- Created: 2024-10-02T10:37:33.380Z

- Posts: 9

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Maria.Suveges (2024-10-02T10:37:33.419Z)

Username: sueveges
Cluster: baobab
Module: R
Job: 12811948
Starting to work with R on baobab, I tried the hello.R / launchR.sh scripts. Removing the version numbers from the modules loaded (so that I get the latest version of everything) the script I tried to run is the following:
```
#!/bin/bash
#SBATCH --job-name=testR
#SBATCH --cpus-per-task=1
#SBATCH --time=0:15:0
#SBATCH --partition=debug-cpu

module load GCC OpenMPI R

INFILE=hello.R
OUTFILE=hello.Rout

srun R CMD BATCH $INFILE $OUTFILE
```
I get two outfiles, the attached `slurm-12811948.out` and `hello.Rout` files. R stops with an error (`hello.Rout`), and slurm terminates with an exit code 1. Same result if I try with a modified R script, just trying to print out exp(1:3) instead of “hello world”.
What am I doing wrong?
Thank you in advance!
Best,
Maria
admin edit: code block syntax[code block syntax](https://commonmark.org/help/)


## Post 2 by @Yann.Sagon (2024-10-04T13:18:49.285Z)

Dear @Maria.Suveges[@Maria.Suveges](https://hpc-community.unige.ch/u/maria.suveges)
I don’t see the two attached files? I tried the example file you probably picked up from here[here](https://gitlab.unige.ch/hpc/softs/-/tree/master/r/R) and it worked.
Best
Yann


## Post 3 by @Maria.Suveges (2024-10-07T08:12:23.380Z)

Hi @Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon)
Sorry, indeed I forgot them. I tried again, with the same result. As I don’t find the way to upload other than image files, here are the files:
hello.Rout:
```
R version 4.2.1 (2022-06-23) -- "Funny-Looking Kid"
Copyright (C) 2022 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu (64-bit)

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

  Natural language support but running in an English locale

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

Error: 3:1: unexpected symbol
2: 
3: local
   ^
Execution halted
```
Is it possible that the R source file has some unwanted characters? I work under Windows 10, and used Notepad++, set to Unix (LF) and UTF-8 encoding. I don’t think I have changed anything in that file, but I looked at it on my laptop and uploaded from there to my homedir. Here is the version on my laptop:
hello.R:
`print("Hello World!")`
Finally, just for completeness, slurm-12866214.out:
```
Inactive Modules:
  1) nettle/3.5.1

Due to MODULEPATH changes, the following have been reloaded:
  1) Szip/2.1.1     2) bzip2/1.0.8     3) libpng/1.6.37

The following have been reloaded with a version change:
  1) FFTW/3.3.8 => FFTW/3.3.10
  2) GCC/8.3.0 => GCC/11.3.0
  3) GCCcore/8.3.0 => GCCcore/11.3.0
  4) GLib/2.62.0 => GLib/2.72.1
  5) GMP/6.1.2 => GMP/6.2.1
  6) GSL/2.6 => GSL/2.7
  7) Ghostscript/9.50 => Ghostscript/9.56.1
  8) HDF5/1.10.5 => HDF5/1.12.2
  9) ICU/64.2 => ICU/71.1
 10) ImageMagick/7.0.9-5 => ImageMagick/7.1.0-37
 11) JasPer/2.0.14 => JasPer/2.0.33
 12) LLVM/9.0.0 => LLVM/14.0.3
 13) LibTIFF/4.0.10 => LibTIFF/4.3.0
 14) LittleCMS/2.9 => LittleCMS/2.13.1
 15) Mesa/19.1.7 => Mesa/22.0.3
 16) NASM/2.14.02 => NASM/2.15.05
 17) NLopt/2.6.1 => NLopt/2.7.1
 18) OpenBLAS/0.3.7 => OpenBLAS/0.3.20
 19) OpenMPI/3.1.4 => OpenMPI/4.1.4
 20) PCRE/8.43 => PCRE/8.45
 21) R/3.6.2 => R/4.2.1
 22) SQLite/3.29.0 => SQLite/3.38.3
 23) ScaLAPACK/2.0.2 => ScaLAPACK/2.2.0-fb
 24) Tcl/8.6.9 => Tcl/8.6.12
 25) Tk/8.6.9 => Tk/8.6.12
 26) UDUNITS/2.2.26 => UDUNITS/2.2.28
 27) X11/20190717 => X11/20220504
 28) XZ/5.2.4 => XZ/5.2.5
 29) binutils/2.32 => binutils/2.38
 30) cURL/7.66.0 => cURL/7.83.0
 31) cairo/1.16.0 => cairo/1.17.4
 32) expat/2.2.7 => expat/2.4.8
 33) fontconfig/2.13.1 => fontconfig/2.14.0
 34) freetype/2.10.1 => freetype/2.12.1
 35) gettext/0.20.1 => gettext/0.21
 36) hwloc/1.11.12 => hwloc/2.7.1
 37) libGLU/9.0.1 => libGLU/9.0.2
 38) libdrm/2.4.99 => libdrm/2.4.110
 39) libffi/3.2.1 => libffi/3.4.2
 40) libjpeg-turbo/2.0.3 => libjpeg-turbo/2.1.3
 41) libpciaccess/0.14 => libpciaccess/0.16
 42) libreadline/8.0 => libreadline/8.1.2
 43) libsndfile/1.0.28 => libsndfile/1.1.0
 44) libunwind/1.3.1 => libunwind/1.6.2
 45) libxml2/2.9.9 => libxml2/2.9.13
 46) ncurses/6.1 => ncurses/6.3
 47) numactl/2.0.12 => numactl/2.0.14
 48) pixman/0.38.0 => pixman/0.40.0
 49) util-linux/2.34 => util-linux/2.38
 50) zlib/1.2.11 => zlib/1.2.12

srun: error: cpu001: task 0: Exited with exit code 1
```
Thanks in advance!
Maria


## Post 4 by @Yann.Sagon (2024-10-07T08:47:11.559Z)

Maria.Suveges:
> Finally, just for completeness, slurm-12866214.out:
You have a lot of modules version change. You probably had already other modules loaded before submitting your job. You can issue`module purge` to cleanup your environment prior to submitting the job.
I’ve checked your account: you have R/3.6.2 and its dependency loaded when you login.  This is in done in your `.bashrc`, you should remove the automatic load and logout/login again.
Can you give another try?
Best


## Post 5 by @Maria.Suveges (2024-10-07T09:48:33.402Z)

Thanks! I did so, now all those info lines about reloaded modules went away, all that remains in the slurm output is the former last line:
`srun: error: cpu001: task 0: Exited with exit code 1`
The .Rout file stayed the same, R stopped with an unexpected symbol message.
Cheers, M


## Post 6 by @Yann.Sagon (2024-10-07T09:59:07.980Z)

Hello, I got it!
the issue is the content of your file `~/.Rprofile`. It wasn’t very clear in our doc[doc](https://doc.eresearch.unige.ch/hpc/applications_and_libraries#r_packages) that the first line was just a way to show the content of the file and mustn’t appear in the file itself, I’ve updated the doc.
In short: remove the first line in your `~/.Rprofile` and it should work!


## Post 7 by @Maria.Suveges (2024-10-07T10:12:27.847Z)

Thanks, it works now! Best, Maria


## Post 8 by @Maria.Suveges (2024-10-09T11:22:19.232Z)

Hi Yann, sorry, just a comment more: I went to the doc again, and saw this line after the now corrected .Rprofile content:
`The first line is purely informative. The output may break things such as package installation. Feel free to comment out this line or remove it.`
But I guess now this too should be removed, right?


## Post 9 by @Yann.Sagon (2024-10-09T13:05:11.943Z)

Maria.Suveges:
> But I guess now this too should be removed, right?
Indeed, thanks for the notification.
