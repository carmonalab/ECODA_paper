# Spyder on Yggrasil

- Source: https://hpc-community.unige.ch/t/3357

- Created: 2024-03-07T09:30:39.969Z

- Posts: 12

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Damien.Benis (2024-03-07T09:30:40.028Z)

Dear admins and users,
I have huge trouble using MNE-Python on Yggdrasil. In baobab, a simple code allows me to use spyder after MNE-Python pip installation of a MNE environnement and Spyder installed using anaconda-navigator
```
module load Anaconda3
conda activate mne
spyder
```
In Yggdrasil I run into errors, you cannot install spyder in this environnement or Failed to initialize GLX. The node attribution have the flag X11
Is there somthing I am doing wrong?
Best
Damien


## Post 2 by @Damien.Benis (2024-03-07T09:44:48.753Z)

Update : l’output terminal aprés une clean install de Spyder sur un environnement dédié
```
X server does not support XInput 2
QStandardPaths: XDG_RUNTIME_DIR not set, defaulting to '/tmp/runtime-benis'
qglx_findConfig: Failed to finding matching FBConfig for QSurfaceFormat(version 2.0, options QFlags<QSurfaceFormat::FormatOption>(), depthBufferSize -1, redBufferSize 1, greenBufferSize 1, blueBufferSize 1, alphaBufferSize -1, stencilBufferSize -1, samples -1, swapBehavior QSurfaceFormat::SingleBuffer, swapInterval 1, colorSpace QSurfaceFormat::DefaultColorSpace, profile  QSurfaceFormat::NoProfile)
qglx_findConfig: Failed to finding matching FBConfig for QSurfaceFormat(version 2.0, options QFlags<QSurfaceFormat::FormatOption>(), depthBufferSize -1, redBufferSize 1, greenBufferSize 1, blueBufferSize 1, alphaBufferSize -1, stencilBufferSize -1, samples -1, swapBehavior QSurfaceFormat::SingleBuffer, swapInterval 1, colorSpace QSurfaceFormat::DefaultColorSpace, profile  QSurfaceFormat::NoProfile)
Could not initialize GLX
```
admin edit: code syntax


## Post 3 by @Damien.Benis (2024-03-07T09:51:31.096Z)

I tried
```
(spyder) (yggdrasil)-[benis@login1 ~]$ module load GCC/9.3.0  OpenMPI/4.0.3
(spyder) (yggdrasil)-[benis@login1 ~]$ module load Spyder/4.1.5
(spyder) (yggdrasil)-[benis@login1 ~]$ spyder
```
Results :
```
Traceback (most recent call last):
  File "/home/benis/.conda/envs/spyder/bin/spyder", line 7, in <module>
    from spyder.app.start import main
  File "/opt/ebsofts/Spyder/4.1.5-foss-2020a-Python-3.8.2/lib/python3.8/site-packages/spyder/app/start.py", line 21, in <module>
    import zmq
  File "/opt/ebsofts/IPython/7.15.0-foss-2020a-Python-3.8.2/lib/python3.8/site-packages/zmq/__init__.py", line 47, in <module>
    from zmq import backend
  File "/opt/ebsofts/IPython/7.15.0-foss-2020a-Python-3.8.2/lib/python3.8/site-packages/zmq/backend/__init__.py", line 40, in <module>
    reraise(*exc_info)
  File "/opt/ebsofts/IPython/7.15.0-foss-2020a-Python-3.8.2/lib/python3.8/site-packages/zmq/utils/sixcerpt.py", line 34, in reraise
    raise value
  File "/opt/ebsofts/IPython/7.15.0-foss-2020a-Python-3.8.2/lib/python3.8/site-packages/zmq/backend/__init__.py", line 27, in <module>
    _ns = select_backend(first)
          ^^^^^^^^^^^^^^^^^^^^^
  File "/opt/ebsofts/IPython/7.15.0-foss-2020a-Python-3.8.2/lib/python3.8/site-packages/zmq/backend/select.py", line 28, in select_backend
    mod = __import__(name, fromlist=public_api)
          ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/opt/ebsofts/IPython/7.15.0-foss-2020a-Python-3.8.2/lib/python3.8/site-packages/zmq/backend/cython/__init__.py", line 6, in <module>
    from . import (constants, error, message, context,
ImportError: cannot import name 'constants' from partially initialized module 'zmq.backend.cython' (most likely due to a circular import) (/opt/ebsofts/IPython/7.15.0-foss-2020a-Python-3.8.2/lib/python3.8/site-packages/zmq/backend/cython/__init__.py)
(spyder) (yggdrasil)-[benis@login1 ~]$
```


## Post 4 by @Gael.Rossignol (2024-04-04T09:16:30.158Z)

Dear Damien,
Did you solve the issue, I have done some research but this is not trivial. Let me know if you need more help.
Best regards,


## Post 5 by @Damien.Benis (2024-06-13T14:25:25.967Z)

Update: I solved the problem. I really dont get it. It should not change anything but it does :slight_smile:
> ssh -Y cpu100
> module load Anaconda3
> module load GCC/9.3.0  OpenMPI/4.0.3
> module load Spyder/4.1.5-Python-3.8.2
> module unload Spyder/4.1.5-Python-3.8.2
> conda activate mne
> spyder
In case someone runs into the same problem.
Best
Damien


## Post 6 by @Damien.Benis (2024-06-14T09:11:43.351Z)

Update no, it loads a spyder from the base, no access to MNE. A potential explanation : The installed `conda` version should be `23.10.0` or newer… Perhaps a conda update --all by the administrator could solve the problem?
https://mne.tools/stable/install/manual_install.html[https://mne.tools/stable/install/manual_install.html](https://mne.tools/stable/install/manual_install.html)


## Post 7 by @Damien.Benis (2024-06-14T12:24:25.496Z)

As an alternative, I tried  Jupyter Lab, but the tab keep crashing, and I do not know why.


## Post 8 by @Yann.Sagon (2024-06-17T09:57:20.145Z)

Here it is:
New software installed: Anaconda3 version 2024.02-1[New software installed: Anaconda3 version 2024.02-1](https://hpc-community.unige.ch/t/new-software-installed-anaconda3-version-2024-02-1/3493) HPC ChangeLog[HPC ChangeLog](https://hpc-community.unige.ch/c/hpc-announce/hpc-changelog/9)
> Dear users, we have installed a new software: Anaconda3 2024.02-1: 

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  Anaconda3: Anaconda3/2024.02-1
--------------------------------------------------------------------------------------------------------------------------------------------------------------------…


## Post 9 by @Damien.Benis (2024-06-18T14:07:44.836Z)

I loaded the module, now the same error appears on Baobab, i cannot use Spyder in Baobab either
> X server does not support XInput 2
> QStandardPaths: XDG_RUNTIME_DIR not set, defaulting to ‘/scratch/runtime-benis’
> qglx_findConfig: Failed to finding matching FBConfig for QSurfaceFormat(version 2.0, options QFlagsQSurfaceFormat::FormatOption(), depthBufferSize -1, redBufferSize 1, greenBufferSize 1, blueBufferSize 1, alphaBufferSize -1, stencilBufferSize -1, samples -1, swapBehavior QSurfaceFormat::SingleBuffer, swapInterval 1, colorSpace QSurfaceFormat::DefaultColorSpace, profile  QSurfaceFormat::NoProfile)
> qglx_findConfig: Failed to finding matching FBConfig for QSurfaceFormat(version 2.0, options QFlagsQSurfaceFormat::FormatOption(), depthBufferSize -1, redBufferSize 1, greenBufferSize 1, blueBufferSize 1, alphaBufferSize -1, stencilBufferSize -1, samples -1, swapBehavior QSurfaceFormat::SingleBuffer, swapInterval 1, colorSpace QSurfaceFormat::DefaultColorSpace, profile  QSurfaceFormat::NoProfile)
> Could not initialize GLX
> Aborted (core dumped)


## Post 10 by @Yann.Sagon (2024-06-18T15:21:12.683Z)

This is what I tried:
- connect to openondemand[openondemand](https://doc.eresearch.unige.ch/hpc/how_to_use_openondemand), launch a desktop session on a compute node
- open a terminal and type `module load Anaconda3`
- launch spyder on the terminal `spyder`
It opens the graphical windows without issue.
Best


## Post 11 by @Damien.Benis (2024-06-19T06:24:43.069Z)

Ok, I never used openondemand, I will try, thank you for the reply.


## Post 12 by @Damien.Benis (2024-07-24T09:23:16.979Z)

Hello,
I found a solution,
> conda deactivate
> conda activate mne
I needed to use the command deactivate before activating my environnment. Using only activate created conflict.
Best
