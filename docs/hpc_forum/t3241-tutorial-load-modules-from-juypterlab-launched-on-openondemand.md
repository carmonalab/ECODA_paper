# [tutorial] load modules from JuypterLab launched on OpenOnDemand

- Source: https://hpc-community.unige.ch/t/3241

- Created: 2024-01-08T14:13:27.065Z

- Posts: 1

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yann.Sagon (2024-01-08T14:13:27.169Z)

When you start a JupyterLab session using OpenOnDemand, you won’t be able to specify the modules you want to load as you would do it from the command line. To get around this problem, we have installed a JupyterLab plugin called jupyterlmod[jupyterlmod](https://github.com/cmd-ntrf/jupyter-lmod) that allows you to list the modules that are already loaded and to load and unload the available modules. See the screenshot below.
image
image445×763 16 KB
[image445×763 16 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/50f0a2dc499e7b8a5b4fffe99d41397d83b764e0.png)
This plugin doesn’t do a `module spider` (search), but a `module avail` (software available, depending on which compiler/openmpi is loaded).
Example: You would like to load the `SciPy-bundle` module.
First, check what are the `SciPy-bundle` dependencies of the same GCCcore that is loaded when starting JupyterLab. In our case it is `GCCcore/11.3.0`.
List the available `SciPy-bundle` versions (on a console):
```
(baobab)-[root@admin1 ~]$ ml spider SciPy-bundle

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  SciPy-bundle:
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Description:
      Bundle of Python packages for scientific software.

     Versions:
        SciPy-bundle/2019.03
        SciPy-bundle/2019.10-Python-2.7.16
        SciPy-bundle/2019.10-Python-3.7.4
        SciPy-bundle/2020.03-Python-3.8.2
        SciPy-bundle/2020.11
        SciPy-bundle/2021.05
        SciPy-bundle/2021.10
        SciPy-bundle/2022.05
        SciPy-bundle/2022.11
        SciPy-bundle/2023.02
        SciPy-bundle/2023.07
```
We now check which version corresponds to `GCCcore/11.3.0` (by accident, we found it on the first try :magic_wand: )
```
(baobab)-[root@admin1 ~]$ ml spider SciPy-bundle/2022.05

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  SciPy-bundle: SciPy-bundle/2022.05
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Description:
      Bundle of Python packages for scientific software

    You will need to load all module(s) on any one of the lines below before the "SciPy-bundle/2022.05" module is available to load.

      GCC/11.3.0  OpenMPI/4.1.4
```
Return to JupyterLab and load `GCC/11.3.0` and `OpenMPI/4.1.4`. To do so, click on the step 3 after searching the module on step 1. Once done, you can see all the available modules, such as SciPy-bundle and load them.
Enjoy!
