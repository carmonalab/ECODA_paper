# Need to install GDAL for python (admin privilege required)

- Source: https://hpc-community.unige.ch/t/4286

- Created: 2026-04-27T13:58:30.228Z

- Tags: baobab

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Philippe.Glass (2026-04-27T13:58:30.323Z)

## Primary informations
Username: glassp
Cluster: cpu298.baobab
Subject: software | need to install GDAL for python (admin privilege required)
jobid: 7845327
## Description
Dear, I am trying to install gdal for python as I need to perform some operations on TIFF raster files that contain geodata. This is my first time using the HPC environment, which is why I need to install all the libraries required to run the machine learning programs, some of which are applied to TIFF images. Do we need to submit a request for each installation that requires admin privileges?
## Steps to Reproduce
I am trying  to launch the 3 installation commands via the VSCode application from the terminal window:
sudo dnf install -y epel-release
sudo dnf install -y gdal
sudo dnf install -y python3-gdal
(execution path: /home/users/g/glassp/ai4sweng-KIO6-model-unit-tester/ModelUnitTester/test)
## Expected Result
installation of GDAL (GDAL present in pip list)
## Actual Result
The execution is denied after each password request.
“We trust you have received the usual lecture from the local System
Administrator. It usually boils down to these three things:
```
#1) Respect the privacy of others.

#2) Think before you type.

#3) With great power comes great responsibility.
```
[sudo] password for glassp:
Sorry, try again.“
Although I consider the privilege issue to be resolved, I am not yet certain that this will lead to the complete installation of the GDAL library for python.
Thank  you for your support.
Best regards,


## Post 2 by @Adrien.Albert (2026-04-27T14:17:00.894Z)

Dear @Philippe.Glass[@Philippe.Glass](https://hpc-community.unige.ch/u/philippe.glass)
GDAL is already installed on cluster and available via module:
```
(baobab)-[alberta@login1 ~]$ ml spider gdal

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  GDAL:
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Description:
      GDAL is a translator library for raster geospatial data formats that is released under an X/MIT style Open Source license by the Open Source Geospatial Foundation. As a library, it presents a single
      abstract data model to the calling application for all supported formats. It also comes with a variety of useful commandline utilities for data translation and processing.

     Versions:
        GDAL/2.0.2
        GDAL/2.1.0
        GDAL/2.2.3-Python-3.6.4
        GDAL/2.2.3-Python-3.6.6
        GDAL/3.0.0-Python-3.7.2
        GDAL/3.0.2-Python-3.7.4
        GDAL/3.0.4-Python-3.8.2
        GDAL/3.2.1
        GDAL/3.3.0-Python-3.9.5
        GDAL/3.3.0
        GDAL/3.3.2
        GDAL/3.5.0
        GDAL/3.6.2
        GDAL/3.7.1
        GDAL/3.9.0
        GDAL/3.10.0
```
Let me know if something is not working correctly.
