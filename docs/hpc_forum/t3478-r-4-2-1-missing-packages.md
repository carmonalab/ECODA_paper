# R/4.2.1: missing packages?

- Source: https://hpc-community.unige.ch/t/3478

- Created: 2024-06-07T10:42:21.374Z

- Posts: 11

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Matthieu.Stigler (2024-06-07T10:42:21.417Z)

It seems that the latest R/4.2.1 version is missing some packages, although they are mentioned in the announce post New software installed R 4.2.1[New software installed R 4.2.1](https://hpc-community.unige.ch/t/new-software-installed-r-4-2-1/2332) ?
Take for example `tidyverse` and `sf`, both are on the list but not installed? It seems it is an omission, as `tidyverse` was there for almost all versions (check `find /opt/ebsofts/R -type d -name "*tidyverse*" -print`).
Also, when running `module spider R/4.2.1`, zero `Included extensions` are listed, unlike previous versions (and unlike in  the pos[pos](https://hpc-community.unige.ch/t/new-software-installed-r-4-2-1/2332)t?)
### Minimum reproducible code
```
module load GCC/11.3.0  OpenMPI/4.1.4 R/4.2.1
R
pkgs_here <- installed.packages(lib.loc="/opt/ebsofts/R/4.2.1-foss-2022a/lib64/R/library")[,"Package"]
any(pkgs_here %in% c("tidyverse", "sf"))
```
will return `FALSE` whereas it should be `TRUE`


## Post 2 by @Matthieu.Stigler (2024-06-17T06:34:01.285Z)

Just following-up on this, besides missing packages (`tidyverse`, `sf`), it seems there are also packages that were not well installed?
See the empty folder for package `terra` (`ls /opt/ebsofts/R/4.2.1-foss-2022a/lib64/R/library/terra`), which is not empty for other versions (check `ls /opt/ebsofts/R/4.2.0-foss-2021b/lib64/R/library/terra` for example).
The problem seems to be the same on yggdrasil, with the same issues (missing packages, weird terra folder), presumably installed with the same problematic script?


## Post 3 by @Adrien.Albert (2024-06-17T11:32:54.817Z)

Hi Matthieu,
I just installed the latest version of R, could you try it and let me know if it’s working well ?
New software installed: R version 4.3.3[New software installed: R version 4.3.3](https://hpc-community.unige.ch/t/new-software-installed-r-version-4-3-3/3494/1) HPC ChangeLog[HPC ChangeLog](https://hpc-community.unige.ch/c/hpc-announce/hpc-changelog/9)
> Dear users, we have installed a new software: R 4.3.3: 

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  R: R/4.3.3
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Descriptio…


## Post 4 by @Matthieu.Stigler (2024-06-17T14:30:36.116Z)

Oh nice, great to hear about the new installation!
Unfortunately, it looks like the problem was “solved” by offering much less extensions than before? I count up to 1000 packages installed with `4.2.0`, but only ~110 with `4.3.3`?
Ideally, do you think you could maintain consistency with the same packages pre-compiled across R versions? If not, would it be possible to include at least the package `tidyverse` with R 4.3.3? It has become the most popular R package, so is likely so be useful to many users here!?
Thanks!
For a rough count of packages:
```
find /opt/ebsofts/R/4.3.3-*/lib64/R/library/* -maxdepth 0 -type d | wc -l
find /opt/ebsofts/R/4.2.0-*/lib64/R/library/* -maxdepth 0 -type d | wc -l
```


## Post 5 by @Matthieu.Stigler (2024-06-17T14:33:51.598Z)

Also, note that R `4.3.3` is not the latest, as in the meanwhile `4.4.0` and `4.4.1` came out (ok, the latter very recently, on June 14) :wink:


## Post 6 by @Adrien.Albert (2024-06-20T10:43:24.348Z)

Dear @Matthieu.Stigler[@Matthieu.Stigler](https://hpc-community.unige.ch/u/matthieu.stigler),
We are using EasyBuild[EasyBuild](https://easybuild.io/) to build software on our HPC cluster. Each build recipe is submitted through a strict validation process to ensure reproducibility across different systems.
We follow their release as closely as possible to avoid any unexpected behavior.
I built `R/4.3.3` using the latest release of EasyBuild, published last month and available here[here](https://github.com/easybuilders/easybuild-easyconfigs/blob/develop/easybuild/easyconfigs/r/R/R-4.3.3-gfbf-2023b.eb). The configuration for `R/4.4.4` is not yet available.
Building software can sometimes be a challenging task (as easybuilder once told me, “a pain in the ass”).
---
SoI tried to install `tidyvers` with `R/4.3.3` and it’s not working due to version depenencies:
```
install.packages("tidyverse")
```
logs:
```
ERROR: dependency ‘Matrix’ is not available for package ‘mgcv’
* removing ‘/home/users/a/alberta/Rpackages/mgcv’
ERROR: dependencies ‘MASS’, ‘mgcv’ are not available for package ‘ggplot2’
* removing ‘/home/users/a/alberta/Rpackages/ggplot2’
ERROR: dependency ‘ggplot2’ is not available for package ‘tidyverse’
* removing ‘/home/users/a/alberta/Rpackages/tidyverse’
```
---
I notice another R-bundle-cran module exists  (which include `tidyvers`) but only compatible with `R/4.3.2`. I am bulding both, R and the bundle.
I keep you informed


## Post 7 by @Matthieu.Stigler (2024-06-20T14:08:40.927Z)

The problem is that the most recent version of `Matrix` depends on R 4.4 (see CRAN: Package Matrix[CRAN: Package Matrix](https://cran.r-project.org/web/packages/Matrix/index.html)).
So to install manually, you need to run:
```
install.packages("lattice") # might be necessary
install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.6-5.tar.gz", repo=NULL)
```
I’m not familiar with EasyBuild but it seems strange to me that they omitted `Matrix`, which is a package used by many other packages, 50 depend on it (check `length(tools::dependsOnPkgs("Matrix", dependencies = "strong"))`) and it is therefore the 19th (out of 20000) most important package in that metric.


## Post 8 by @Adrien.Albert (2024-06-21T07:44:24.006Z)

Hi @Matthieu.Stigler[@Matthieu.Stigler](https://hpc-community.unige.ch/u/matthieu.stigler),
The `Matrix`, `tidyverse`, and other R packages are considered extensions and are, by definition, optional. It makes sense not to include these packages in the core `R` build for several reasons.
- By keeping the core and extensions separate, the core R installation remains lighter and more modular.
- If an error occurs in building an extension package, it does not cause the entire R build process to fail. This separation ensures that issues with individual packages do not affect the stability of the core R installation.
- Users have the flexibility to install (or load other module) only the packages they need, which can be tailored to their specific use cases.
A more telling example:
It took me 30 minutes to compile R/4.3.2 without extensions.
Yesterday I started compiling R-bundle-cran, which includes almost 1100 extensions, it’s still a work in progress… :slight_smile:


## Post 9 by @Adrien.Albert (2024-06-21T12:55:36.562Z)

Here the R module build with the extension:
New software installed: R version 4.3.2[New software installed: R version 4.3.2](https://hpc-community.unige.ch/t/new-software-installed-r-version-4-3-2/3504) HPC ChangeLog[HPC ChangeLog](https://hpc-community.unige.ch/c/hpc-announce/hpc-changelog/9)
> Dear users, we have installed a new software: R 4.3.2: 

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  R: R/4.3.2
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Descriptio…
New software installed: R-bundle-CRAN 2023.12[New software installed: R-bundle-CRAN 2023.12](https://hpc-community.unige.ch/t/new-software-installed-r-bundle-cran-2023-12/3505) HPC ChangeLog[HPC ChangeLog](https://hpc-community.unige.ch/c/hpc-announce/hpc-changelog/9)
> Dear users, 
we have installed a new software: R-bundle-CRAN 2023.12: 

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  R-bundle-CRAN: R-bundle-CRAN/2023.12
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------…


## Post 10 by @Matthieu.Stigler (2024-06-21T14:32:31.653Z)

oh very nice, thanks!! That seems a great solution to have both R and R-bundle separately.
I see that R-bundle also comes with R (ie I don’t need to `module load R`), though a slightly lower version 4.3.2? Not sure if you can edit the output of `module spider`, but might be a useful info to add?
Merci encore!


## Post 11 by @Adrien.Albert (2024-06-21T15:04:44.692Z)

Hi @Matthieu.Stigler[@Matthieu.Stigler](https://hpc-community.unige.ch/u/matthieu.stigler)
You do not need to load R before as R as been marked as dependency. If you load another R version it will automatically adapt all module version used.
```
(baobab)-[alberta@login1 ~]$  ml GCC/12.3.0  OpenMPI/4.1.5 R-bundle-CRAN/2023.12
(baobab)-[alberta@login1 ~]$ ml

Currently Loaded Modules:
  1) EasyBuild/4.9.2    17) OpenMPI/4.1.5       33) gzip/1.12        49) Tcl/8.6.13              65) GMP/6.2.1             81) Gdk-Pixbuf/2.42.10    97) PROJ/9.2.0                 113) Xerces-C++/3.2.4
  2) GCCcore/12.3.0     18) OpenBLAS/0.3.23     34) lz4/1.9.4        50) SQLite/3.42.0           66) NLopt/2.7.1           82) Pango/1.50.14         98) libgeotiff/1.7.1           114) Imath/3.1.7
  3) zlib/1.2.13        19) FlexiBLAS/3.3.1     35) zstd/1.5.5       51) NASM/2.16.01            67) libogg/1.3.5          83) libepoxy/1.5.10       99) cffi/1.15.1                115) OpenEXR/3.1.7
  4) binutils/2.40      20) FFTW/3.3.10         36) libdrm/2.4.115   52) libjpeg-turbo/2.1.5.1   68) FLAC/1.4.2            84) Wayland/1.22.0       100) cryptography/41.0.1        116) Highway/1.0.4
  5) GCC/12.3.0         21) FFTW.MPI/3.3.10     37) libglvnd/1.6.0   53) jbigkit/2.1             69) libvorbis/1.3.7       85) GTK3/3.24.37         101) virtualenv/20.23.1         117) Brunsli/0.1
  6) numactl/2.0.16     22) ScaLAPACK/2.2.0-fb  38) libunwind/1.6.2  54) libdeflate/1.18         70) libopus/1.4           86) Ghostscript/10.01.2  102) Python-bundle-PyPI/2023.06 118) Qhull/2020.2
  7) XZ/5.4.2           23) bzip2/1.0.8         39) LLVM/16.0.6      55) LibTIFF/4.5.0           71) LAME/3.100            87) JasPer/4.0.0         103) pybind11/2.11.1            119) LERC/4.0.0
  8) libxml2/2.11.4     24) expat/2.5.0         40) Mesa/23.1.4      56) Java/11 -> Java/11.0.2  72) libsndfile/1.2.2      88) LittleCMS/2.15       104) SciPy-bundle/2023.07       120) OpenJPEG/2.5.0
  9) libpciaccess/0.17  25) libpng/1.6.39       41) libGLU/9.0.3     57) PCRE/8.45               73) Szip/2.1.1            89) ImageMagick/7.1.1-15 105) libtirpc/1.3.3             121) SWIG/4.1.1
 10) hwloc/2.9.1        26) Brotli/1.0.9        42) pixman/0.42.2    58) libgit2/1.7.1           74) HDF5/1.14.0           90) GLPK/5.0             106) HDF/4.2.16-2               122) GDAL/3.7.1
 11) OpenSSL/1.1        27) freetype/2.13.0     43) libffi/3.4.4     59) cURL/8.0.1              75) UDUNITS/2.2.28        91) nodejs/18.17.1       107) Boost/1.82.0               123) MPFR/4.2.0
 12) libevent/2.1.12    28) ncurses/6.4         44) gettext/0.21.1   60) Tk/8.6.13               76) GSL/2.7               92) Python/3.11.3        108) arpack-ng/3.9.0            124) PostgreSQL/16.1
 13) UCX/1.14.1         29) util-linux/2.39     45) PCRE2/10.42      61) ICU/73.2                77) ATK/2.38.0            93) netCDF/4.9.2         109) Armadillo/12.6.2           125) R-bundle-CRAN/2023.12
 14) libfabric/1.18.0   30) fontconfig/2.14.2   46) GLib/2.77.1      62) HarfBuzz/5.3.1          78) DBus/1.15.4           94) GEOS/3.12.0          110) CFITSIO/4.3.0
 15) PMIx/4.2.4         31) xorg-macros/1.20.0  47) cairo/1.17.8     63) FriBidi/1.0.12          79) at-spi2-core/2.49.91  95) libarchive/3.6.2     111) giflib/5.2.1
 16) UCC/1.2.0          32) X11/20230603        48) libreadline/8.2  64) R/4.3.2                 80) at-spi2-atk/2.38.0    96) nlohmann_json/3.11.2 112) json-c/0.16
```
I wonder if this bundle wouldn’t be compatible with higher versions of R. :thinking:
However, it was also compiled with a defined version of GCC and OpenMPI, so it would have to be recompiled with the correct versions. :melting_face:
