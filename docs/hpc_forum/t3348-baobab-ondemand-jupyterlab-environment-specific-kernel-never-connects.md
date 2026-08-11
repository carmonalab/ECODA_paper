# Baobab OnDemand: jupyterLab environment-specific kernel never connects

- Source: https://hpc-community.unige.ch/t/3348

- Created: 2024-03-01T10:55:49.059Z

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Erica.Lastufka (2024-03-01T10:55:49.162Z)

## Primary informations
Username: lastufka
Cluster: Baobab
## Description
When trying to launch a notebook with a custom ipykernel in Baobab OnDemand JupyterLab, the kernel never connects.
## Steps to Reproduce
- launch JupyterLab via Baobab OnDemand
- try to create new notebook with custom kernel, or open existing notebook and select that kernel
## Expected Result
Kernel connects
## Actual Result
Kernel fails to connect. When trying the default ipykernel instead, basic modules like numpy are not found. Using module load [either numpy versions and the required modules] from the notebook results in an error. JupyterLab is unusable for running Python.


## Post 2 by @Gael.Rossignol (2024-04-04T09:46:43.964Z)

Dear Erica,
Yould you please check the tutorial :
[tutorial] load modules from JuypterLab launched on OpenOnDemand[[tutorial] load modules from JuypterLab launched on OpenOnDemand](https://hpc-community.unige.ch/t/tutorial-load-modules-from-juypterlab-launched-on-openondemand/3241) HPC Technical[HPC Technical](https://hpc-community.unige.ch/c/hpc-technical/5)
> When you start a JupyterLab session using OpenOnDemand, you won’t be able to specify the modules you want to load as you would do it from the command line. To get around this problem, we have installed a JupyterLab plugin called jupyterlmod[jupyterlmod](https://github.com/cmd-ntrf/jupyter-lmod) that allows you to list the modules that are already loaded and to load and unload the available modules. See the screenshot below. 
 [image][[image]](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/50f0a2dc499e7b8a5b4fffe99d41397d83b764e0.png) 
This plugin doesn’t do a module spider (search), but a module avail (software available, depending on which compiler…
It could solve your issue.
Best regards,
