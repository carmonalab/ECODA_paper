# How to use DeepLabCut?

- Source: https://hpc-community.unige.ch/t/3385

- Created: 2024-03-21T12:51:45.161Z

- Posts: 2

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Maxime.Revel (2024-03-21T12:51:45.220Z)

Hello,
I’m trying to use DeepLabCut on the Cluster. I have a little bit of experience with it on my personal machine as I already trained my model and analyzed a few videos, but I would like to take advantage of the GPUs of the cluster to speed up the analysis.
I imported the DeepLabCut/2.2.0.6 module (It’s already a win by my standards)
I’m not entirely sure how to proceed now, and I would appreciate a lot if anyone could share their experience using this module.
Am I supposed to run it using python or using the GUI?
Do I need to create a conda environment?
Thank a lot for your help,
Maxime


## Post 2 by @Yann.Sagon (2024-04-03T06:52:52.363Z)

Hi, I’m a beginner too with DeepLabCut as I never used it, I only installed it on the cluster but this is very different stuff. How I would proceed:
Load the module:
```
(baobab)-[sagon@login2 ~]$ ml GCC/10.3.0  OpenMPI/4.1.1 DeepLabCut/2.2.0.6-CUDA-11.3.1
```
Then I want to “see” what this module is doing:
```
(baobab)-[sagon@login2 ~]$ ml show DeepLabCut
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
   /opt/ebmodules/all/MPI/GCC/10.3.0/OpenMPI/4.1.1/DeepLabCut/2.2.0.6-CUDA-11.3.1.lua:
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
help([[
Description
===========
Markerless tracking of user-defined features with deep learning

More information
================
 - Homepage: http://www.mousemotorlab.org/deeplabcut

Included extensions
===================
DeepLabCut-2.2.0.6, filterpy-1.4.5, msgpack-numpy-0.4.8, tensorpack-0.11,
tf_slim-1.1.0
]])
whatis("Description: Markerless tracking of user-defined features with deep learning")
whatis("Homepage: http://www.mousemotorlab.org/deeplabcut")
whatis("URL: http://www.mousemotorlab.org/deeplabcut")
whatis("Extensions: DeepLabCut-2.2.0.6, filterpy-1.4.5, msgpack-numpy-0.4.8, tensorpack-0.11, tf_slim-1.1.0")
conflict("DeepLabCut")
prepend_path("CMAKE_PREFIX_PATH","/opt/ebsofts/DeepLabCut/2.2.0.6-foss-2021a-CUDA-11.3.1")
prepend_path("LIBRARY_PATH","/opt/ebsofts/DeepLabCut/2.2.0.6-foss-2021a-CUDA-11.3.1/lib")
prepend_path("PATH","/opt/ebsofts/DeepLabCut/2.2.0.6-foss-2021a-CUDA-11.3.1/bin")
setenv("EBROOTDEEPLABCUT","/opt/ebsofts/DeepLabCut/2.2.0.6-foss-2021a-CUDA-11.3.1")
setenv("EBVERSIONDEEPLABCUT","2.2.0.6")
setenv("EBDEVELDEEPLABCUT","/opt/ebsofts/DeepLabCut/2.2.0.6-foss-2021a-CUDA-11.3.1/easybuild/MPI-GCC-10.3.0-OpenMPI-4.1.1-DeepLabCut-2.2.0.6-CUDA-11.3.1-easybuild-devel")
prepend_path("PYTHONPATH","/opt/ebsofts/DeepLabCut/2.2.0.6-foss-2021a-CUDA-11.3.1/lib/python3.9/site-packages")
setenv("EBEXTSLISTDEEPLABCUT","filterpy-1.4.5,tensorpack-0.11,tf_slim-1.1.0,msgpack-numpy-0.4.8,DeepLabCut-2.2.0.6")
```
Let see what binaries this module provides:
```
(baobab)-[sagon@login2 ~]$ ls $EBROOTDEEPLABCUT/bin
dlc  download.sh
```
Unfortunately `dlc` binary doesn’t seems to work?! Anyway, you can import deeplacut from Python:
```
(baobab)-[sagon@gpu021 2.2.0.6-foss-2021a-CUDA-11.3.1]$ python
Python 3.9.5 (default, Aug 31 2021, 16:56:58)
[GCC 10.3.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import deeplabcut
DLC loaded in light mode; you cannot use any GUI (labeling, relabeling and standalone GUI)
>>>
```
Or more modern using `ipython`:
```
(baobab)-[sagon@gpu021 2.2.0.6-foss-2021a-CUDA-11.3.1]$ ipython
Python 3.9.5 (default, Aug 31 2021, 16:56:58)
Type 'copyright', 'credits' or 'license' for more information
IPython 7.25.0 -- An enhanced Interactive Python. Type '?' for help.

In [1]: import deeplabcut
DLC loaded in light mode; you cannot use any GUI (labeling, relabeling and standalone GUI)

In [2]: deeplabcut.__version__
Out[2]: '2.2.0.6'

In [3]:
```
Is this ok like that?
Best
Yann
