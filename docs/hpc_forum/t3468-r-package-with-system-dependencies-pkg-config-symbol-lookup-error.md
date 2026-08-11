# R package with system dependencies: pkg-config: symbol lookup error

- Source: https://hpc-community.unige.ch/t/3468

- Created: 2024-06-03T12:30:53.784Z

- Posts: 6

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Matthieu.Stigler (2024-06-03T12:30:53.830Z)

I just started using the cluster, trying use R with the standard package tidyverse.
Upon using install.packages(“tidyverse”), there is an error due to the fact that the R-dependency `textshaping` has a system dependency.
> Configuration failed to find the harfbuzz freetype2 fribidi library. Try installing:
> - deb: libharfbuzz-dev libfribidi-dev (Debian, Ubuntu, etc)
> - rpm: harfbuzz-devel fribidi-devel (Fedora, EPEL)
> - csw: libharfbuzz_dev libfribidi_dev (Solaris)
> - brew: harfbuzz fribidi (OSX)
So I tried to load these dependencies with:
> module load GCCcore/12.2.0 HarfBuzz/5.3.1
> module load GCCcore/13.2.0 FriBidi/1.0.13
But I still get a (slightly different) error:
> pkg-config: symbol lookup error: pkg-config: undefined symbol: pkgconf_pkg_dir_list_build
> pkg-config: symbol lookup error: pkg-config: undefined symbol: pkgconf_pkg_dir_list_build
> Using PKG_CFLAGS=
> Using PKG_LIBS=-lfreetype -lharfbuzz -lfribidi -lpng
> If harfbuzz freetype2 fribidi is already installed, check that ‘pkg-config’ is in your
> PATH and PKG_CONFIG_PATH contains a harfbuzz freetype2 fribidi.pc file. If pkg-config
> is unavailable you can set INCLUDE_DIR and LIB_DIR manually via:
> R CMD INSTALL --configure-vars=‘INCLUDE_DIR=… LIB_DIR=…’
Could you please advise me on what is the `pkg-config: symbol lookup error: pkg-config: undefined symbol: pkgconf_pkg_dir_list_build` error, and what should be done to have
`install.packages("textshaping")`
to run?
Thanks!
### Minimal reprex:
```
module load GCC/11.3.0  OpenMPI/4.1.4 R/4.2.1
module load GCCcore/12.2.0 HarfBuzz/5.3.1
module load GCCcore/13.2.0 FriBidi/1.0.13
R
install.packages("textshaping")
```


## Post 2 by @Yann.Sagon (2024-06-06T14:13:05.391Z)

Dear @Matthieu.Stigler[@Matthieu.Stigler](https://hpc-community.unige.ch/u/matthieu.stigler)
you can’t load the software that are requesting different version of GCC or OpenMPI at the same time: you need to stick the the same common base for all the software you load.
For your case, this should work:
```
(baobab)-[sagon@login2 ~]$ ml purge
(baobab)-[sagon@login2 ~]$ ml GCC/11.3.0  OpenMPI/4.1.4 R/4.2.1 HarfBuzz/4.2.1 FriBidi/1.0.12
```
Best


## Post 3 by @Matthieu.Stigler (2024-06-06T15:33:31.350Z)

Thanks a lot @Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon) for looking into it!
Unfortunately, I still see the same error message, despite using your new line, see below output and a minimal reproducible example.
I wonder whether I am missing the `-dev` versions of HarfBuzz or FriBidi, noting that the error message mentions:
> Configuration failed to find the harfbuzz freetype2 fribidi library. Try installing:
> - deb: libharfbuzz-dev libfribidi-dev (Debian, Ubuntu, etc)
> - rpm: harfbuzz-devel fribidi-devel (Fedora, EPEL)
> - csw: libharfbuzz_dev libfribidi_dev (Solaris)
Is it possible that baobab’s `FriBidi/1.0.12` corresponds to the non-devel version of `harfbuzz-devel`? How can I check that, in particular whether the file `hb-ft.h` is included?
Thanks!!
### minimum reproducible code
```
module load GCC/11.3.0 OpenMPI/4.1.4 R/4.2.1 HarfBuzz/4.2.1 FriBidi/1.0.12
R
install.packages("textshaping") 
```
### error message
> pkg-config: symbol lookup error: pkg-config: undefined symbol: pkgconf_pkg_dir_list_build
> pkg-config: symbol lookup error: pkg-config: undefined symbol: pkgconf_pkg_dir_list_build
> Using PKG_CFLAGS=
> Using PKG_LIBS=-lfreetype -lharfbuzz -lfribidi -lpng
> --------------------------- [ANTICONF] --------------------------------
> Configuration failed to find the harfbuzz freetype2 fribidi library. Try installing:
> - deb: libharfbuzz-dev libfribidi-dev (Debian, Ubuntu, etc)
> - rpm: harfbuzz-devel fribidi-devel (Fedora, EPEL)
> - csw: libharfbuzz_dev libfribidi_dev (Solaris)
> - brew: harfbuzz fribidi (OSX)
> - If harfbuzz freetype2 fribidi is already installed, check that ‘pkg-config’ is in your
> - PATH and PKG_CONFIG_PATH contains a harfbuzz freetype2 fribidi.pc file	. If pkg-config
> - is unavailable you can set INCLUDE_DIR and LIB_DIR manually via:
> - R CMD INSTALL --configure-vars=‘INCLUDE_DIR=… LIB_DIR=…’
> - -------------------------- [ERROR MESSAGE] ---------------------------
> - :1:10: fatal error: hb-ft.h: No such file or directory
> - compilation terminated.
---


## Post 4 by @Yann.Sagon (2024-06-07T07:06:38.730Z)

Hello,
Matthieu.Stigler:
> How can I check that, in particular whether the file `hb-ft.h` is included?
You can check like that:
```
(baobab)-[testsagon1@login2 ~]$ ls $EBROOTHARFBUZZ/include/harfbuzz/hb-ft.h
/opt/ebsofts/HarfBuzz/4.2.1-GCCcore-11.3.0/include/harfbuzz/hb-ft.h
```
It seems the error is due to
```
pkg-config: symbol lookup error: pkg-config: undefined symbol: pkgconf_pkg_dir_list_build
```
You are using `pkg-config` from the OS which is version 1.4.2. I tried with `pkg-config` from EasyBuild and I could install the package successfully (this time I tested it!)
```
module load GCC/11.3.0 OpenMPI/4.1.4 R/4.2.1 HarfBuzz/4.2.1 FriBidi/1.0.12 pkg-config/0.29.2
```
Best
Yann


## Post 5 by @Matthieu.Stigler (2024-06-07T09:29:58.746Z)

Excellent, that did the trick, thanks a lopt @Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon) !!
As an aside, I’m slightly confused by the package numbers… are we here really downgrading pkg-config from  1.4.2.  to 0.29, or are these different version numbers?
Thanks!


## Post 6 by @Yann.Sagon (2024-07-17T14:42:55.914Z)

Matthieu.Stigler:
> I’m slightly confused by the package numbers… are we here really downgrading pkg-config from 1.4.2. to 0.29, or are these different version numbers?
I have no idea, I didn’t checked.
