# Compile with OpenGL

- Source: https://hpc-community.unige.ch/t/3568

- Created: 2024-07-25T11:57:42.506Z

- Posts: 6

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yuzheng.Kang (2024-07-25T11:57:42.587Z)

## Primary informations
Username: Kang
Cluster: Yggdrasil
## Description
I am wondering if there is an option or an alternative way to compile a program with OpenGL?
```
cmake-3.30/Modules/FindPackageHandleStandardArgs.cmake:233 (message):
  Could NOT find OpenGL (missing: OPENGL_opengl_LIBRARY OPENGL_glx_LIBRARY
  OPENGL_INCLUDE_DIR)
```


## Post 2 by @Adrien.Albert (2024-07-25T12:50:40.424Z)

Hi @Yuzheng.Kang[@Yuzheng.Kang](https://hpc-community.unige.ch/u/yuzheng.kang)
I do not have the context, could give me more information about what you want to do ?
Do you have a sbatch ? a command to execute ?


## Post 3 by @Yuzheng.Kang (2024-07-25T20:45:06.924Z)

Hi @Adrien.Albert[@Adrien.Albert](https://hpc-community.unige.ch/u/Adrien.Albert)
Thanks for your response.
I tried to compile a software called colmap https://colmap.github.io/install.html[https://colmap.github.io/install.html](https://colmap.github.io/install.html) on the cluster. To do that I submit a .sh file which as you said request computing node using #SBATCH.
Basically what I did is running these code. For the cmake part, I used module load to get most of the packages I need, but I missed OpenGL module.
```
git clone https://github.com/colmap/colmap.git
cd colmap
mkdir build
cd build
cmake .. -GNinja
ninja
sudo ninja install
```
Then there was an error message about `cmake-3.30/Modules/FindPackageHandleStandardArgs.cmake:233 (message):   Could NOT find OpenGL (missing: OPENGL_opengl_LIBRARY OPENGL_glx_LIBRARY   OPENGL_INCLUDE_DIR)`


## Post 4 by @Adrien.Albert (2024-07-25T22:03:02.918Z)

Hi @Yuzheng.Kang[@Yuzheng.Kang](https://hpc-community.unige.ch/u/yuzheng.kang),
Please note that the use of ‘sudo’ is not permitted on the cluster.
I will proceed with installing COLMAP 3.8. The build receipe are available on EasyBuild, the tool we use to compile software on our HPC cluster.


## Post 5 by @Yuzheng.Kang (2024-07-26T07:27:20.865Z)

Hi @Adrien.Albert[@Adrien.Albert](https://hpc-community.unige.ch/u/Adrien.Albert)
Thanks for the help!


## Post 6 by @Adrien.Albert (2024-08-06T08:34:21.736Z)

Hi @Yuzheng.Kang[@Yuzheng.Kang](https://hpc-community.unige.ch/u/yuzheng.kang)
I forgot to answer:
New software installed: COLMAP version 3.8[New software installed: COLMAP version 3.8](https://hpc-community.unige.ch/t/new-software-installed-colmap-version-3-8/3571) HPC ChangeLog[HPC ChangeLog](https://hpc-community.unige.ch/c/hpc-announce/hpc-changelog/9)
> Dear users, we have installed a new software: COLMAP 3.8: 

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  COLMAP: COLMAP/3.8
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
   …
