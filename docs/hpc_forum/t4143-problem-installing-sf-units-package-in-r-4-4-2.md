# Problem installing ‘sf’ & 'units' package in R/4.4.2

- Source: https://hpc-community.unige.ch/t/4143

- Created: 2025-11-10T13:24:42.948Z

- Tags: bamboo

- Posts: 9

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Fabien.Cottier (2025-11-10T13:24:43.032Z)

## Primary informations
Username: $cottierf
Cluster: $bamboo
## Description
Hello. I am trying to install the packages sf, and its dependencies “units” on R/4.4.2
When I trying to run the command, I am running into the error, below. Spefically, when attempting to load the package after built, R throws a `*** caught illegal operation ***`
`address 0x15004256bc0d, cause illegal operand` error.  No sure what is going wrong…
## Steps to Reproduce
> module load GCC/13.3.0 OpenMPI/5.0.3 GEOS/3.12.2 GDAL/3.10.0 PROJ/9.4.1 CMake/3.29.3 UDUNITS/2.2.26 R/4.4.2
> R
> install.packages(“sf”)
## Expected Result
The R sf packages includes the critical dependency R “units” should install correctly. However, it fails (see error message below and sessionInfo() outputs at the bottom).
## Actual Results
Error message received, when installing R dependency package “units” (depends on system library udunits2):
> trying URL ‘https://stat.ethz.ch/CRAN/src/contrib/units_1.0-0.tar.gz’[https://stat.ethz.ch/CRAN/src/contrib/units_1.0-0.tar.gz’](https://stat.ethz.ch/CRAN/src/contrib/units_1.0-0.tar.gz%E2%80%99)
> Content type ‘application/x-gzip’ length 367179 bytes (358 KB)
> downloaded 358 KB
> trying URL ‘https://stat.ethz.ch/CRAN/src/contrib/sf_1.0-22.tar.gz’[https://stat.ethz.ch/CRAN/src/contrib/sf_1.0-22.tar.gz’](https://stat.ethz.ch/CRAN/src/contrib/sf_1.0-22.tar.gz%E2%80%99)
> Content type ‘application/x-gzip’ length 9874803 bytes (9.4 MB)
> downloaded 9.4 MB
> - installing source package ‘units’ …
> - ** package ‘units’ successfully unpacked and MD5 sums checked
> - ** using staged installation
> - configure: units: 1.0-0
> - checking for default C++ standard… 201703
> - checking whether the C++ compiler works… yes
> - checking for C++ compiler default output file name… a.out
> - checking for suffix of executables…
> - checking whether we are cross compiling… no
> - checking for suffix of object files… o
> - checking whether the compiler supports GNU C++… yes
> - checking whether g++ -std=gnu++17 accepts -g… yes
> - checking for g++ -std=gnu++17 option to enable C++11 features… none needed
> - checking for stdio.h… yes
> - checking for stdlib.h… yes
> - checking for string.h… yes
> - checking for inttypes.h… yes
> - checking for stdint.h… yes
> - checking for strings.h… yes
> - checking for sys/stat.h… yes
> - checking for sys/types.h… yes
> - checking for unistd.h… yes
> - checking for _Bool… no
> - checking for stdbool.h that conforms to C99 or later… yes
> - checking for error_at_line… yes
> - checking for gcc… gcc
> - checking whether the compiler supports GNU C… yes
> - checking whether gcc accepts -g… yes
> - checking for gcc option to enable C11 features… none needed
> - checking for brew… no
> - checking for XML_ParserCreate in -lexpat… yes
> - checking for udunits2.h… yes
> - checking for ut_read_xml in -ludunits2… yes
> - configure: creating ./config.status
> - config.status: creating src/Makevars
> - ** libs
> - using C++ compiler: ‘g++ (GCC) 13.3.0’
> - g++ -std=gnu++17 -I"/opt/ebsofts/R/4.4.2-gfbf-2024a/lib64/R/include" -DNDEBUG -DUDUNITS2_DIR=0  -I’/home/users/c/cottierf/R/x86_64-pc-linux-gnu-library/4.4/Rcpp/include’ -I/opt/ebsofts/FriBidi/1.0.15-GCCcore-13.3.0/include -I/opt/ebsofts/HarfBuzz/9.0.0-GCCcore-13.3.0/include -I/opt/ebsofts/Tk/8.6.14-GCCcore-13.3.0/include -I/opt/ebsofts/cURL/8.7.1-GCCcore-13.3.0/include -I/opt/ebsofts/OpenSSL/3/include -I/opt/ebsofts/libgit2/1.8.1-GCCcore-13.3.0/include -I/opt/ebsofts/Java/17.0.2/include -I/opt/ebsofts/LibTIFF/4.6.0-GCCcore-13.3.0/include -I/opt/ebsofts/libjpeg-turbo/3.0.1-GCCcore-13.3.0/include -I/opt/ebsofts/libpng/1.6.43-GCCcore-13.3.0/include -I/opt/ebsofts/PCRE2/10.43-GCCcore-13.3.0/include -I/opt/ebsofts/SQLite/3.45.3-GCCcore-13.3.0/include -I/opt/ebsofts/zlib/1.3.1-GCCcore-13.3.0/include -I/opt/ebsofts/XZ/5.4.5-GCCcore-13.3.0/include -I/opt/ebsofts/bzip2/1.0.8-GCCcore-13.3.0/include -I/opt/ebsofts/ncurses/6.5-GCCcore-13.3.0/include -I/opt/ebsofts/libreadline/8.2-GCCcore-13.3.0/include -I/opt/ebsofts/cairo/1.18.0-GCCcore-13.3.0/include -I/opt/ebsofts/libGLU/9.0.3-GCCcore-13.3.0/include -I/opt/ebsofts/Mesa/24.1.3-GCCcore-13.3.0/include -I/opt/ebsofts/X11/20240607-GCCcore-13.3.0/include -I/opt/ebsofts/pkgconf/2.2.0-GCCcore-13.3.0/include -I/opt/ebsofts/FFTW/3.3.10-GCC-13.3.0/include -I/opt/ebsofts/FlexiBLAS/3.4.4-GCC-13.3.0/include -I/opt/ebsofts/FlexiBLAS/3.4.4-GCC-13.3.0/include/flexiblas    -fpic  -O2 -ftree-vectorize -march=native -fno-math-errno   -c RcppExports.cpp -o RcppExports.o
> - g++ -std=gnu++17 -I"/opt/ebsofts/R/4.4.2-gfbf-2024a/lib64/R/include" -DNDEBUG -DUDUNITS2_DIR=0  -I’/home/users/c/cottierf/R/x86_64-pc-linux-gnu-library/4.4/Rcpp/include’ -I/opt/ebsofts/FriBidi/1.0.15-GCCcore-13.3.0/include -I/opt/ebsofts/HarfBuzz/9.0.0-GCCcore-13.3.0/include -I/opt/ebsofts/Tk/8.6.14-GCCcore-13.3.0/include -I/opt/ebsofts/cURL/8.7.1-GCCcore-13.3.0/include -I/opt/ebsofts/OpenSSL/3/include -I/opt/ebsofts/libgit2/1.8.1-GCCcore-13.3.0/include -I/opt/ebsofts/Java/17.0.2/include -I/opt/ebsofts/LibTIFF/4.6.0-GCCcore-13.3.0/include -I/opt/ebsofts/libjpeg-turbo/3.0.1-GCCcore-13.3.0/include -I/opt/ebsofts/libpng/1.6.43-GCCcore-13.3.0/include -I/opt/ebsofts/PCRE2/10.43-GCCcore-13.3.0/include -I/opt/ebsofts/SQLite/3.45.3-GCCcore-13.3.0/include -I/opt/ebsofts/zlib/1.3.1-GCCcore-13.3.0/include -I/opt/ebsofts/XZ/5.4.5-GCCcore-13.3.0/include -I/opt/ebsofts/bzip2/1.0.8-GCCcore-13.3.0/include -I/opt/ebsofts/ncurses/6.5-GCCcore-13.3.0/include -I/opt/ebsofts/libreadline/8.2-GCCcore-13.3.0/include -I/opt/ebsofts/cairo/1.18.0-GCCcore-13.3.0/include -I/opt/ebsofts/libGLU/9.0.3-GCCcore-13.3.0/include -I/opt/ebsofts/Mesa/24.1.3-GCCcore-13.3.0/include -I/opt/ebsofts/X11/20240607-GCCcore-13.3.0/include -I/opt/ebsofts/pkgconf/2.2.0-GCCcore-13.3.0/include -I/opt/ebsofts/FFTW/3.3.10-GCC-13.3.0/include -I/opt/ebsofts/FlexiBLAS/3.4.4-GCC-13.3.0/include -I/opt/ebsofts/FlexiBLAS/3.4.4-GCC-13.3.0/include/flexiblas    -fpic  -O2 -ftree-vectorize -march=native -fno-math-errno   -c tokenizer.cpp -o tokenizer.o
> - g++ -std=gnu++17 -I"/opt/ebsofts/R/4.4.2-gfbf-2024a/lib64/R/include" -DNDEBUG -DUDUNITS2_DIR=0  -I’/home/users/c/cottierf/R/x86_64-pc-linux-gnu-library/4.4/Rcpp/include’ -I/opt/ebsofts/FriBidi/1.0.15-GCCcore-13.3.0/include -I/opt/ebsofts/HarfBuzz/9.0.0-GCCcore-13.3.0/include -I/opt/ebsofts/Tk/8.6.14-GCCcore-13.3.0/include -I/opt/ebsofts/cURL/8.7.1-GCCcore-13.3.0/include -I/opt/ebsofts/OpenSSL/3/include -I/opt/ebsofts/libgit2/1.8.1-GCCcore-13.3.0/include -I/opt/ebsofts/Java/17.0.2/include -I/opt/ebsofts/LibTIFF/4.6.0-GCCcore-13.3.0/include -I/opt/ebsofts/libjpeg-turbo/3.0.1-GCCcore-13.3.0/include -I/opt/ebsofts/libpng/1.6.43-GCCcore-13.3.0/include -I/opt/ebsofts/PCRE2/10.43-GCCcore-13.3.0/include -I/opt/ebsofts/SQLite/3.45.3-GCCcore-13.3.0/include -I/opt/ebsofts/zlib/1.3.1-GCCcore-13.3.0/include -I/opt/ebsofts/XZ/5.4.5-GCCcore-13.3.0/include -I/opt/ebsofts/bzip2/1.0.8-GCCcore-13.3.0/include -I/opt/ebsofts/ncurses/6.5-GCCcore-13.3.0/include -I/opt/ebsofts/libreadline/8.2-GCCcore-13.3.0/include -I/opt/ebsofts/cairo/1.18.0-GCCcore-13.3.0/include -I/opt/ebsofts/libGLU/9.0.3-GCCcore-13.3.0/include -I/opt/ebsofts/Mesa/24.1.3-GCCcore-13.3.0/include -I/opt/ebsofts/X11/20240607-GCCcore-13.3.0/include -I/opt/ebsofts/pkgconf/2.2.0-GCCcore-13.3.0/include -I/opt/ebsofts/FFTW/3.3.10-GCC-13.3.0/include -I/opt/ebsofts/FlexiBLAS/3.4.4-GCC-13.3.0/include -I/opt/ebsofts/FlexiBLAS/3.4.4-GCC-13.3.0/include/flexiblas    -fpic  -O2 -ftree-vectorize -march=native -fno-math-errno   -c udunits.cpp -o udunits.o
> - g++ -std=gnu++17 -shared -L/opt/ebsofts/R/4.4.2-gfbf-2024a/lib64/R/lib -L/opt/ebsofts/FriBidi/1.0.15-GCCcore-13.3.0/lib64 -L/opt/ebsofts/FriBidi/1.0.15-GCCcore-13.3.0/lib -L/opt/ebsofts/HarfBuzz/9.0.0-GCCcore-13.3.0/lib64 -L/opt/ebsofts/HarfBuzz/9.0.0-GCCcore-13.3.0/lib -L/opt/ebsofts/Tk/8.6.14-GCCcore-13.3.0/lib64 -L/opt/ebsofts/Tk/8.6.14-GCCcore-13.3.0/lib -L/opt/ebsofts/cURL/8.7.1-GCCcore-13.3.0/lib64 -L/opt/ebsofts/cURL/8.7.1-GCCcore-13.3.0/lib -L/opt/ebsofts/OpenSSL/3/lib64 -L/opt/ebsofts/OpenSSL/3/lib -L/opt/ebsofts/libgit2/1.8.1-GCCcore-13.3.0/lib64 -L/opt/ebsofts/libgit2/1.8.1-GCCcore-13.3.0/lib -L/opt/ebsofts/Java/17.0.2/lib64 -L/opt/ebsofts/Java/17.0.2/lib -L/opt/ebsofts/LibTIFF/4.6.0-GCCcore-13.3.0/lib64 -L/opt/ebsofts/LibTIFF/4.6.0-GCCcore-13.3.0/lib -L/opt/ebsofts/libjpeg-turbo/3.0.1-GCCcore-13.3.0/lib64 -L/opt/ebsofts/libjpeg-turbo/3.0.1-GCCcore-13.3.0/lib -L/opt/ebsofts/libpng/1.6.43-GCCcore-13.3.0/lib64 -L/opt/ebsofts/libpng/1.6.43-GCCcore-13.3.0/lib -L/opt/ebsofts/PCRE2/10.43-GCCcore-13.3.0/lib64 -L/opt/ebsofts/PCRE2/10.43-GCCcore-13.3.0/lib -L/opt/ebsofts/SQLite/3.45.3-GCCcore-13.3.0/lib64 -L/opt/ebsofts/SQLite/3.45.3-GCCcore-13.3.0/lib -L/opt/ebsofts/zlib/1.3.1-GCCcore-13.3.0/lib64 -L/opt/ebsofts/zlib/1.3.1-GCCcore-13.3.0/lib -L/opt/ebsofts/XZ/5.4.5-GCCcore-13.3.0/lib64 -L/opt/ebsofts/XZ/5.4.5-GCCcore-13.3.0/lib -L/opt/ebsofts/bzip2/1.0.8-GCCcore-13.3.0/lib64 -L/opt/ebsofts/bzip2/1.0.8-GCCcore-13.3.0/lib -L/opt/ebsofts/ncurses/6.5-GCCcore-13.3.0/lib64 -L/opt/ebsofts/ncurses/6.5-GCCcore-13.3.0/lib -L/opt/ebsofts/libreadline/8.2-GCCcore-13.3.0/lib64 -L/opt/ebsofts/libreadline/8.2-GCCcore-13.3.0/lib -L/opt/ebsofts/cairo/1.18.0-GCCcore-13.3.0/lib64 -L/opt/ebsofts/cairo/1.18.0-GCCcore-13.3.0/lib -L/opt/ebsofts/libGLU/9.0.3-GCCcore-13.3.0/lib64 -L/opt/ebsofts/libGLU/9.0.3-GCCcore-13.3.0/lib -L/opt/ebsofts/Mesa/24.1.3-GCCcore-13.3.0/lib64 -L/opt/ebsofts/Mesa/24.1.3-GCCcore-13.3.0/lib -L/opt/ebsofts/X11/20240607-GCCcore-13.3.0/lib64 -L/opt/ebsofts/X11/20240607-GCCcore-13.3.0/lib -L/opt/ebsofts/pkgconf/2.2.0-GCCcore-13.3.0/lib64 -L/opt/ebsofts/pkgconf/2.2.0-GCCcore-13.3.0/lib -L/opt/ebsofts/FFTW/3.3.10-GCC-13.3.0/lib64 -L/opt/ebsofts/FFTW/3.3.10-GCC-13.3.0/lib -L/opt/ebsofts/FlexiBLAS/3.4.4-GCC-13.3.0/lib64 -L/opt/ebsofts/FlexiBLAS/3.4.4-GCC-13.3.0/lib -L/opt/ebsofts/GCCcore/13.3.0/lib64 -L/opt/ebsofts/GCCcore/13.3.0/lib -o units.so RcppExports.o tokenizer.o udunits.o -lexpat -lexpat -ludunits2 -L/opt/ebsofts/R/4.4.2-gfbf-2024a/lib64/R/lib -lR
> - installing to /home/users/c/cottierf/R/x86_64-pc-linux-gnu-library/4.4/00LOCK-units/00new/units/libs
> - ** R
> - ** demo
> - ** inst
> - ** byte-compile and prepare package for lazy loading
> - ** help
> - *** installing help indices
> - ** building package indices
> - ** installing vignettes
> - ** testing if installed package can be loaded from temporary location
> *** caught illegal operation ***
> address 0x15004256bc0d, cause ‘illegal operand’
> Traceback:
> 1: ud_init(path)
> 2: load_units_xml()
> 3: fun(libname, pkgname)
> 4: doTryCatch(return(expr), name, parentenv, handler)
> 5: tryCatchOne(expr, names, parentenv, handlers[[1L]])
> 6: tryCatchList(expr, classes, parentenv, handlers)
> 7: tryCatch(fun(libname, pkgname), error = identity)
> 8: runHook(“.onLoad”, env, package.lib, package)
> 9: loadNamespace(package, lib.loc)
> 10: doTryCatch(return(expr), name, parentenv, handler)
> 11: tryCatchOne(expr, names, parentenv, handlers[[1L]])
> 12: tryCatchList(expr, classes, parentenv, handlers)
> 13: tryCatch({    attr(package, “LibPath”) ← which.lib.loc    ns ← loadNamespace(package, lib.loc)    env ← attachNamespace(ns, pos = pos, deps, exclude, include.only)}, error = function(e) {    P ← if (!is.null(cc ← conditionCall(e)))         paste(" in", deparse(cc)[1L])    else “”    msg ← gettextf(“package or namespace load failed for %s%s:\n %s”,         sQuote(package), P, conditionMessage(e))    if (logical.return && !quietly)         message(paste(“Error:”, msg), domain = NA)    else stop(msg, call. = FALSE, domain = NA)})
> 14: library(pkg_name, lib.loc = lib, character.only = TRUE, logical.return = TRUE)
> 15: withCallingHandlers(expr, packageStartupMessage = function(c) tryInvokeRestart(“muffleMessage”))
> 16: suppressPackageStartupMessages(library(pkg_name, lib.loc = lib,     character.only = TRUE, logical.return = TRUE))
> 17: doTryCatch(return(expr), name, parentenv, handler)
> 18: tryCatchOne(expr, names, parentenv, handlers[[1L]])
> 19: tryCatchList(expr, classes, parentenv, handlers)
> 20: tryCatch(expr, error = function(e) {    call ← conditionCall(e)    if (!is.null(call)) {        if (identical(call[[1L]], quote(doTryCatch)))             call ← sys.call(-4L)        dcall ← deparse(call, nlines = 1L)        prefix ← paste(“Error in”, dcall, ": ")        LONG ← 75L        sm ← strsplit(conditionMessage(e), “\n”)[[1L]]        w ← 14L + nchar(dcall, type = “w”) + nchar(sm[1L], type = “w”)        if (is.na(w))             w ← 14L + nchar(dcall, type = “b”) + nchar(sm[1L],                 type = “b”)        if (w > LONG)             prefix ← paste0(prefix, "\n  ")    }    else prefix ← "Error : "    msg ← paste0(prefix, conditionMessage(e), “\n”)    .Internal(seterrmessage(msg[1L]))    if (!silent && isTRUE(getOption(“show.error.messages”))) {        cat(msg, file = outFile)        .Internal(printDeferredWarnings())    }    invisible(structure(msg, class = “try-error”, condition = e))})
> 21: try(suppressPackageStartupMessages(library(pkg_name, lib.loc = lib,     character.only = TRUE, logical.return = TRUE)))
> 22: tools:::.test_load_package(“units”, “/home/users/c/cottierf/R/x86_64-pc-linux-gnu-library/4.4/00LOCK-units/00new”)
> An irrecoverable exception occurred. R is aborting now …
> ERROR: loading failed
> - removing ‘/home/users/c/cottierf/R/x86_64-pc-linux-gnu-library/4.4/units’
> - ERROR: dependency ‘units’ is not available for package ‘sf’
> - removing ‘/home/users/c/cottierf/R/x86_64-pc-linux-gnu-library/4.4/sf’
> The downloaded source packages are in
> ‘/tmp/RtmpmgD2OA/downloaded_packages’
> Warning messages:
> 1: In install.packages(“sf”) :
> installation of package ‘units’ had non-zero exit status
> 2: In install.packages(“sf”) :
> installation of package ‘sf’ had non-zero exit status
Some additional informations:
> R> sessionInfo()
> R version 4.4.2 (2024-10-31)
> Platform: x86_64-pc-linux-gnu
> Running under: Rocky Linux 9.6 (Blue Onyx)
> Matrix products: default
> BLAS/LAPACK: FlexiBLAS OPENBLAS;  LAPACK version 3.12.0
> locale:
> [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
> [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
> [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
> [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
> [9] LC_ADDRESS=C               LC_TELEPHONE=C
> [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C
> time zone: Europe/Zurich
> tzcode source: system (glibc)
> attached base packages:
> [1] stats     graphics  grDevices utils     datasets  methods   base
> other attached packages:
> [1] Rcpp_1.1.0
> loaded via a namespace (and not attached):
> [1] compiler_4.4.2 tools_4.4.2


## Post 2 by @Adrien.Albert (2025-11-13T10:46:40.958Z)

hi @Fabien.Cottier[@Fabien.Cottier](https://hpc-community.unige.ch/u/fabien.cottier)
Indeed I get an error with installing units. However I needed to load `CMake` in addition of R, to compile depedencies:
```
ml GCC/13.3.0 R CMake
```
Here my install R script, but the `units-1.0-1` seems incompatible or (missing something I do not know) installing the previous version works:
```
$ cat install_R_pkg.R

installed.packages()
install.packages("sf", type = "source")
install.packages("https://cran.r-project.org/src/contrib/Archive/units/units_0.8-7.tar.gz", repos = NULL, type = "source", lib="~/Rpackages")

find.package("units", verbose = TRUE)
packageVersion("units")
```
```
(baobab)-[alberta@cpu349 R]$ Rscript install_R_pkg.R
[... Build blabla...]

installing to /home/users/a/alberta/Rpackages/00LOCK-units/00new/units/libs
** R
** demo
** inst
** byte-compile and prepare package for lazy loading
** help
*** installing help indices
** building package indices
** installing vignettes
** testing if installed package can be loaded from temporary location
** checking absolute paths in shared objects and dynamic libraries
** testing if installed package can be loaded from final location
** testing if installed package keeps a record of temporary installation path
* DONE (units)
[1] "/home/users/a/alberta/Rpackages/units"
[1] ‘0.8.7’
```
TADA !


## Post 3 by @Fabien.Cottier (2025-11-14T13:28:28.982Z)

Ciao Adrien,
Merci bien de ton message.J’ai effectivement pu réussir à installer le package “units”, selon tes indications (à la fois la version 0.8-7 et 1.0.0). Il semble que l’erreur soit associée avec le module UDUNITS/2.2.28. Si le module n’est pas “load”, alors “units” est installé correctement.
Je rencontre par contre toujours une erreur avec l’installation du package “sf” à proprement parlé, il semble que le programme ne puisse pas detecté les libraries GDAL, même si le module est chargé (l’erreur semble similaire à celle raportée ici: Problem install 'sf' package in R/4.2.1[Problem install 'sf' package in R/4.2.1](https://hpc-community.unige.ch/t/problem-install-sf-package-in-r-4-2-1/2865) )
Encore merci pour ton assistance,
```
module load foss/2024a  GEOS/3.12.2 GDAL/3.10.0 PROJ/9.4.1  CMake R/4.4.2

R

install.packages(“sf”)
```
Erreur
> installing source package ‘sf’ …
> ** package ‘sf’ successfully unpacked and MD5 sums checked
> ** using staged installation
> configure: CC: gcc
> configure: CXX: g++ -std=gnu++17
> checking for gdal-config… /opt/ebsofts/GDAL/3.10.0-foss-2024a/bin/gdal-config
> checking gdal-config usability… yes
> configure: GDAL: 3.10.0
> checking GDAL version >= 2.0.1… yes
> checking for gcc… gcc
> checking whether the C compiler works… yes
> checking for C compiler default output file name… a.out
> checking for suffix of executables…
> checking whether we are cross compiling… no
> checking for suffix of object files… o
> checking whether the compiler supports GNU C… yes
> checking whether gcc accepts -g… yes
> checking for gcc option to enable C11 features… none needed
> checking for stdio.h… yes
> checking for stdlib.h… yes
> checking for string.h… yes
> checking for inttypes.h… yes
> checking for stdint.h… yes
> checking for strings.h… yes
> checking for sys/stat.h… yes
> checking for sys/types.h… yes
> checking for unistd.h… yes
> checking for gdal.h… yes
> checking GDAL: linking with --libs only… no
> checking GDAL: linking with --libs and --dep-libs… no
> /opt/ebsofts/binutils/2.42-GCCcore-13.3.0/bin/ld: warning: libOpenCL.so.1, needed by /opt/ebsofts/GDAL/3.10.0-foss-2024a/lib/libgdal.so, not found (try using -rpath or -rpath-link)
> /opt/ebsofts/binutils/2.42-GCCcore-13.3.0/bin/ld: /opt/ebsofts/GDAL/3.10.0-foss-2024a/lib/libgdal.so: undefined reference to `clReleaseMemObject@OPENCL_1.0' /opt/ebsofts/binutils/2.42-GCCcore-13.3.0/bin/ld: /opt/ebsofts/GDAL/3.10.0-foss-2024a/lib/libgdal.so: undefined reference to `clGetKernelWorkGroupInfo@OPENCL_1.0’
> [more similar messages]
> collect2: error: ld returned 1 exit status
> /opt/ebsofts/binutils/2.42-GCCcore-13.3.0/bin/ld: cannot find -lOpenCL: No such file or directory
> collect2: error: ld returned 1 exit status
> configure: Install failure: compilation and/or linkage problems.
> configure: error: GDALAllRegister not found in libgdal.
> *** Installing this package from source requires the prior


## Post 4 by @Fabien.Cottier (2025-11-14T15:55:13.232Z)

Here are some supplementary information. The error is caused by GDAL 3.10 missing some critical libraries. In addition, GCC 13.3.0 does not have the library libOpenCL.so
Indeed, running the following command throws an error
> `gdalinfo --version`
> Caught signal 4 (Illegal instruction: illegal operand)
> ==== backtrace (tid:2407155) ====
> 0 0x000000000003ebf0 __GI___sigaction()  :0
> 1 0x00000000001954cd osgeo::proj::util::LocalName::LocalName()  :0
> 2 0x0000000000195568 osgeo::proj::util::NameSpace::createGLOBAL()  :0
> 3 0x00000000000d6df7 __static_initialization_and_destruction_0()  static.cpp:0
> 4 0x000000000000551e call_init()  /usr/src/debug/glibc-2.34-168.el9_6.23.x86_64/elf/dl-init.c:70
> 5 0x000000000000560c _dl_init()  /usr/src/debug/glibc-2.34-168.el9_6.23.x86_64/elf/dl-init.c:117
> 6 0x000000000001da9a _dl_start_user()  :0
On the other hand, I am able to install sf if I load GCC/13.2.0  OpenMPI/4.1.6 and GDAL/3.10.0 instead and then Rv4.4.1 (instead of R4.4.2).
Running `gdalinfo --version`with GDAL 3.9,00 give the following output ":  `GDAL 3.9.0, released 2024/05/07`
I am not sure what is going wrong, but it seems that GCC 13.3.0 and GDAL 3.10 may be corrupted.


## Post 5 by @Yann.Sagon (2025-11-17T13:33:31.375Z)

Dear @Fabien.Cottier[@Fabien.Cottier](https://hpc-community.unige.ch/u/fabien.cottier)
I found something interesting. We are building Baobab software on an apptainer container. And it appears that when we start apptainer with gpu support, the needed libs are automagically available and gdal used them. Later when you run the soft, the library isn’t available anymore.
I’ll rebuild gdal, but I need first to check how to ignore the library OpenCL.


## Post 6 by @Fabien.Cottier (2025-11-17T14:00:54.615Z)

HI Yann,
That’s interesting. I am not enough well versed in Linux sytems to be of much assistance, here. However, while I was trying to find a solution on Friday, chatGPT suggested that I looked for locations where the library libOpenCL.so was installed.
I did not find any within the GCC libraries using the linux shell “find” tool, but openCL was installed on multiple CUDA libraries, e.g. /opt/ebsofts/CUDA/12.4.0/targets/x86_64-linux/lib/libOpenCL.so.
Not sure if this helps.


## Post 7 by @Yann.Sagon (2025-11-18T10:32:53.314Z)

Dear @Fabien.Cottier[@Fabien.Cottier](https://hpc-community.unige.ch/u/fabien.cottier)
I’ve recompiled GDAL, GEOS and PROJ. It should work better now, please give a try.
Best regards
Yann


## Post 8 by @Fabien.Cottier (2025-11-19T08:22:59.266Z)

Hi @yann.Sagon[@yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon),
Thanks. I can confirm that everything is now working. I was able to recompile the sf and other related packages under R 4.4.2 without any issue. It seemst to be working well, as I can judge.
Just of out of curiosity: how did you solve the issue?
Fabien


## Post 9 by @Yann.Sagon (2025-11-19T09:12:12.810Z)

Fabien.Cottier:
> how did you solve the issue
I’ve recompiled the software on a “legacy” compute node without gpu, thus `libOpenCL.so` not present and not linked.
