# CMake was not found on the PATH

- Source: https://hpc-community.unige.ch/t/3689

- Created: 2024-10-17T13:35:59.279Z

- Posts: 3

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Samuel.Orso (2024-10-17T13:35:59.337Z)

## Primary informations
Username: orsos
Cluster: baobab
## Description
I would like to install an R package called `nloptr` that requires cmake. During installation, I encounter the following error:
```
------------------ CMAKE NOT FOUND --------------------

CMake was not found on the PATH. Please install CMake:

 - sudo yum install cmake          (Fedora/CentOS; inside a terminal)
 - sudo apt install cmake          (Debian/Ubuntu; inside a terminal).
 - sudo pacman -S cmake            (Arch Linux; inside a terminal).
 - brew install --cask cmake       (MacOS; inside a terminal with Homebrew)
 - sudo port install cmake         (MacOS; inside a terminal with MacPorts)

Alternatively install CMake from: <https://cmake.org/>
```
## Steps to Reproduce
- Launch R/4.3.3: `ml GCC/13.2.0 R/4.3.3`
- In R, install package `nloptr`: `install.packages("nloptr")`
## Expected Result
Package installation without an error.
## Actual Result
`ERROR: configuration failed for package ‘nloptr’`, reason invoked: `CMake was not found on the PATH`.


## Post 2 by @Adrien.Albert (2024-10-17T17:01:11.596Z)

Hi @Samuel.Orso[@Samuel.Orso](https://hpc-community.unige.ch/u/samuel.orso)
This documentation should interest you:
https://doc.eresearch.unige.ch/hpc/applications_and_libraries#module_-_lmod[https://doc.eresearch.unige.ch/hpc/applications_and_libraries#module_-_lmod](https://doc.eresearch.unige.ch/hpc/applications_and_libraries#module_-_lmod)
Solution for lazy people
```
   .-.
  .-. |
.-. | |
| | | | ___
| | | |/ _ \     ___
|       /o\_'||.' o `.
|       \=/ /|| `---'
 \      ___//  \
  \____/   (`__')
       ///////\\\\\\\
        `--.____.--'  . . . Be seeing you! (forgive me it's late)
```
```
(baobab)-[alberta@login1 ~]$ ml spider cmake

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  CMake:
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Description:
      CMake, the cross-platform, open-source build system. CMake is a family of tools designed to build, test and package software.

     Versions:
        CMake/3.0.0
        CMake/3.4.1
        CMake/3.4.3
        CMake/3.5.2
        CMake/3.6.1
        CMake/3.6.2
        CMake/3.7.2
        CMake/3.9.1
        CMake/3.9.5
        CMake/3.9.6
        CMake/3.10.2
        CMake/3.11.4
        CMake/3.12.1
        CMake/3.13.3
        CMake/3.15.3
        CMake/3.16.4
        CMake/3.18.4
        CMake/3.20.1
        CMake/3.21.1
        CMake/3.22.1
        CMake/3.23.1
        CMake/3.24.3
        CMake/3.25.1
        CMake/3.25.2
        CMake/3.26.3
        CMake/3.27.6

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  For detailed information about a specific "CMake" package (including how to load the modules) use the module's full name.
  Note that names that have a trailing (E) are extensions provided by other modules.
  For example:

     $ module spider CMake/3.27.6
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
```
:joy:.


## Post 3 by @Samuel.Orso (2024-10-17T17:17:27.540Z)

Arf, of course :sweat_smile: Thanks!
