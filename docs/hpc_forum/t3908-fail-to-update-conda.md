# Fail to update conda

- Source: https://hpc-community.unige.ch/t/3908

- Created: 2025-04-03T13:50:57.743Z

- Posts: 7

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Camille.Christe (2025-04-03T13:50:57.788Z)

I try to install a program following the instruction given in the github page : GitHub - cpouchon/ORTHOSKIM: ORTHOSKIM allows in silico capture of targeted sequences in genomic or transcriptomic libraries.[GitHub - cpouchon/ORTHOSKIM: ORTHOSKIM allows in silico capture of targeted sequences in genomic or transcriptomic libraries.](https://github.com/cpouchon/ORTHOSKIM)
I follow the instructions but when asking for
```
conda install -c bioconda spades exonerate diamond blast mafft trimal numpy joblib scipy -y
```
it fails with the message :
Solving environment: failed with initial frozen solve. Retrying with flexible solve.
Solving environment: failed with repodata from current_repodata.json, will retry with next repodata source.
ResolvePackageNotFound:
- python=3.1
Apparently python 3.1 is not compatible with conda 4.10.1
I tried to update conda but it fails as well.
I tried as well after loading the last version of Anaconda3 (Anaconda3/2024.02-1)
I tried two times the same process, removing the virtual environment and recreating it again with the same error messages.
My colleagues are using the same program, no problem to install it. When they ask for the conda version it is the 24.11 by default.
I am not 100% sure that it is the problem but I hope it is.
Best regards,
camille


## Post 2 by @Adrien.Albert (2025-04-09T14:33:54.477Z)

Hi @Camille.Christe[@Camille.Christe](https://hpc-community.unige.ch/u/camille.christe)
Since Anaconda requires the installation/creation of a lot of files, I’m giving you the procedure to create a container image of it.
- Load the module cotainr:
```
(bamboo)-[alberta@login1 ~]$  module load cotainr
```
- Create a file with the following definition:
```
(bamboo)-[alberta@login1 ~]$ cat bioenv.yml 
name: bioenv
channels:
  - bioconda
  - conda-forge
  - defaults
dependencies:
  - _libgcc_mutex=0.1=conda_forge
  - _openmp_mutex=4.5=2_gnu
  - blast=2.16.0=h66d330f_5
  - bzip2=1.0.8=h4bc722e_7
  - c-ares=1.34.5=hb9d3cd8_0
  - ca-certificates=2025.1.31=hbcca054_0
  - curl=8.13.0=h332b0f4_0
  - diamond=2.1.11=h5ca1c30_1
  - entrez-direct=22.4=he881be0_0
  - exonerate=2.4.0=h09da616_8
  - gawk=5.3.1=hcd3d067_0
  - glib=2.84.1=h07242d1_0
  - glib-tools=2.84.1=h4833e2c_0
  - gmp=6.3.0=hac33072_2
  - joblib=1.4.2=pyhd8ed1ab_1
  - keyutils=1.6.1=h166bdaf_0
  - krb5=1.21.3=h659f571_0
  - ld_impl_linux-64=2.43=h712a8e2_4
  - libasprintf=0.23.1=h8e693c7_0
  - libblas=3.9.0=31_h59b9bed_openblas
  - libcblas=3.9.0=31_he106b2a_openblas
  - libcurl=8.13.0=h332b0f4_0
  - libedit=3.1.20250104=pl5321h7949ede_0
  - libev=4.33=hd590300_2
  - libexpat=2.7.0=h5888daf_0
  - libffi=3.4.6=h2dba641_1
  - libgcc=14.2.0=h767d61c_2
  - libgcc-ng=14.2.0=h69a702a_2
  - libgettextpo=0.23.1=h5888daf_0
  - libgfortran=14.2.0=h69a702a_2
  - libgfortran5=14.2.0=hf1ad2bd_2
  - libglib=2.84.1=h2ff4ddf_0
  - libgomp=14.2.0=h767d61c_2
  - libiconv=1.18=h4ce23a2_1
  - libidn2=2.3.8=ha4ef2c3_0
  - liblapack=3.9.0=31_h7ac8fdf_openblas
  - liblzma=5.8.1=hb9d3cd8_0
  - libmpdec=4.0.0=h4bc722e_0
  - libnghttp2=1.64.0=h161d5f1_0
  - libopenblas=0.3.29=pthreads_h94d23a6_0
  - libsqlite=3.49.1=hee588c1_2
  - libssh2=1.11.1=hf672d98_0
  - libstdcxx=14.2.0=h8f9b012_2
  - libstdcxx-ng=14.2.0=h4852527_2
  - libunistring=0.9.10=h7f98852_0
  - libuuid=2.38.1=h0b41bf4_0
  - libxcrypt=4.4.36=hd590300_1
  - libzlib=1.3.1=hb9d3cd8_2
  - mafft=7.525=h031d066_1
  - mpfr=4.2.1=h90cbb55_3
  - ncbi-vdb=3.2.1=h9948957_0
  - ncurses=6.5=h2d0b736_3
  - numpy=2.2.4=py313h103f029_0
  - openssl=3.5.0=h7b32b05_0
  - packaging=24.2=pyhd8ed1ab_2
  - pcre=8.45=h9c3ff4c_0
  - pcre2=10.44=hba22ea6_2
  - perl=5.32.1=7_hd590300_perl5
  - perl-archive-tar=2.40=pl5321hdfd78af_0
  - perl-carp=1.38=pl5321hdfd78af_4
  - perl-common-sense=3.75=pl5321hdfd78af_0
  - perl-compress-raw-bzip2=2.201=pl5321h87f3376_1
  - perl-compress-raw-zlib=2.105=pl5321h87f3376_0
  - perl-encode=3.19=pl5321hec16e2b_1
  - perl-exporter=5.72=pl5321hdfd78af_2
  - perl-exporter-tiny=1.002002=pl5321hdfd78af_0
  - perl-extutils-makemaker=7.70=pl5321hd8ed1ab_0
  - perl-io-compress=2.201=pl5321h503566f_5
  - perl-io-zlib=1.14=pl5321hdfd78af_0
  - perl-json=4.10=pl5321hdfd78af_1
  - perl-json-xs=4.03=pl5321h9948957_4
  - perl-list-moreutils=0.430=pl5321hdfd78af_0
  - perl-list-moreutils-xs=0.430=pl5321h7b50bb2_5
  - perl-parent=0.236=pl5321hdfd78af_2
  - perl-pathtools=3.75=pl5321hec16e2b_3
  - perl-scalar-list-utils=1.62=pl5321hec16e2b_1
  - perl-types-serialiser=1.01=pl5321hdfd78af_0
  - pip=25.0.1=pyh145f28c_0
  - python=3.13.2=h4724d56_1_cp313t
  - python_abi=3.13=5_cp313t
  - readline=8.2=h8c095d6_2
  - rpsbproc=0.5.0=hd6d6fdc_3
  - scipy=1.15.2=py313h7f7b39c_0
  - setuptools=78.1.0=pyhff2d567_0
  - spades=4.1.0=haf24da9_0
  - tk=8.6.13=noxft_h4845f30_101
  - trimal=1.5.0=h9948957_2
  - tzdata=2025b=h78e105d_0
  - wget=1.21.4=hda4d442_0
  - zlib=1.3.1=hb9d3cd8_2
  - zstd=1.5.7=hb8e6e7a_2
prefix: /home/users/a/alberta/.conda/envs/bioenv
```
---
:light_bulb: To generate this file :
- Install the conda environment
```
$ conda create -n bioenv -c bioconda -c conda-forge spades exonerate diamond blast mafft trimal numpy joblib scipy -y
```
- Generate the file (you may need to edit/delete the prefix field)
```
conda activate bioenv
conda env export > bioenv.yml
```
– “Installing the conda env just to generate the file is useless” :speaking_head: you will tell me…
– Yes, but this file can be generated once (by the maintainer, for example).
---
- build the image (I use ubuntu but you can use rockylinux[rockylinux](https://hub.docker.com/_/rockylinux) )
```
$ bin/cotainr build bioenv.sif --base-image=docker://ubuntu:latest --accept-licenses --conda-env=bioenv.yml
```
And then you can use the image like this:
```
(bamboo)-[alberta@login1 ~]$ apptainer exec bioenv.sif python3 -c "import numpy; print(numpy.__version__)"
2.2.4
```
### :green_circle: What the benefit using this method ?
- Instead of having a billion files, you only have one (our storage will thank you)
- It’s exportable (you can share it with your amazing colleague)
- When we update the cluster, it should not impact the image, as all dependencies are inside it
- Because I took the time to make this tutorial :heart_hands:t3:
### :cross_mark: What’s the disadvantage?
- This method does not make fondue for two :fondue:


## Post 3 by @Camille.Christe (2025-04-10T14:08:54.638Z)

Thank you very much for you help and the time spent ! It is working now!
That’s very sad that it does not make foundue for two


## Post 4 by @Adrien.Albert (2025-04-11T07:54:15.440Z)

Is it working with conda or container ? or both ? :slight_smile:


## Post 5 by @Adrien.Albert (2025-04-11T09:30:14.426Z)

Hi @Camille.Christe[@Camille.Christe](https://hpc-community.unige.ch/u/camille.christe)
I have update our documentation about this subject:
doc.eresearch.unige.ch[doc.eresearch.unige.ch](https://doc.eresearch.unige.ch/hpc/applications_and_libraries#how_to_create_a_conda_environment_in_a_container)
### hpc:applications_and_libraries [eResearch Doc][hpc:applications_and_libraries [eResearch Doc]](https://doc.eresearch.unige.ch/hpc/applications_and_libraries#how_to_create_a_conda_environment_in_a_container)
Using cotainr which have been installed on cluster:
New software installed: cotainr version 2025.3.0[New software installed: cotainr version 2025.3.0](https://hpc-community.unige.ch/t/new-software-installed-cotainr-version-2025-3-0/3917) HPC ChangeLog[HPC ChangeLog](https://hpc-community.unige.ch/c/hpc-announce/hpc-changelog/9)
> Dear users, we have installed a new software: cotainr 2025.3.0: 

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  cotainr: cotainr/2025.3.0
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------…


## Post 6 by @Camille.Christe (2025-04-15T06:23:46.242Z)

With the container (apptainer …), I have not tried with conda yet


## Post 7 by @Adrien.Albert (2025-04-15T12:48:38.802Z)

That sound perfect, thank you for feedback. :folded_hands:
PS: if everything works fine, there’s no need to test with conda (it’s even better).
