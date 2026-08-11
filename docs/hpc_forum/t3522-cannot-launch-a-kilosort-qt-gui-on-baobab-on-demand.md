# Cannot launch a Kilosort Qt GUI on baobab on demand

- Source: https://hpc-community.unige.ch/t/3522

- Created: 2024-07-02T10:07:34.012Z

- Posts: 3

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Njiva.Andrianarivelo (2024-07-02T10:07:34.058Z)

Hi Guys,
My name is Andry, postdoc in the Bellone’s Lab.
I’m trying to use `Kilosort`, an electrophysiological spike sorter that use GPU on baobab on demand.
`https://github.com/MouseLand/Kilosort`
Kilosort 4 uses `pyqt6` and `pyqtgraph` for its GUI.
I am using the `PRIVATE-BELONE` node and `gpu046` through a windows desktop.
My problem is that when I launch the program, the GUI doesnt load and show this error:
```
(base) (baobab)-[andrianj@gpu046 ~]$ conda activate ks4
(ks4) (baobab)-[andrianj@gpu046 ~]$ python -m kilosort
qt.qpa.plugin: From 6.5.0, xcb-cursor0 or libxcb-cursor0 is needed to load the Qt xcb platform plugin.
qt.qpa.plugin: Could not load the Qt platform plugin "xcb" in "" even though it was found.
This application failed to start because no Qt platform plugin could be initialized. Reinstalling the application may fix this problem.

Available platform plugins are: vkkhrdisplay, wayland, minimalegl, vnc, offscreen, eglfs, minimal, xcb, wayland-egl, linuxfb.

Aborted (core dumped)
```
I installed `Kilosort 4`  on my local machine and didn’t have any trouble also I don’t have any problem with other software using GUI on baobab on demand (for instance `Phy` (GitHub - cortex-lab/phy: phy: interactive visualization and manual spike sorting of large-scale ephys data[GitHub - cortex-lab/phy: phy: interactive visualization and manual spike sorting of large-scale ephys data](https://github.com/cortex-lab/phy/)) and the display interface support is mainly why I use baobab on demand.
So far I tried:
uninstalling/intalling `pyqt6` and `pyqtgraph`,
downgrading to `pyqt6=6.5.0/ 6.4.0`,
even using `pyqt5` instead,
trying `unset DISPLAY`,
> export XDG_RUNTIME_DIR=/tmp/runtime-andrianj
```
export QT_QPA_PLATFORM=wayland
```
(and each of vkkhrdisplay, minimalegl, vnc, offscreen, eglfs, minimal, xcb, wayland-egl, linuxfb)
Can you please help me ? This step is very important for analyzing my data…
Best,
Andry


## Post 2 by @Adrien.Albert (2024-07-19T18:15:20.256Z)

Hi @Njiva.Andrianarivelo[@Njiva.Andrianarivelo](https://hpc-community.unige.ch/u/njiva.andrianarivelo)
I tried on OpenOnDemand ans it’s working:
image
image1909×906 151 KB
[image1909×906 151 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/fcee577507a97b9d7555d76e6d58590400e74fbc.png)
Baobab OpenOnDemand[Baobab OpenOnDemand](https://openondemand.baobab.hpc.unige.ch/pun/sys/dashboard/) select advanced desktop pick your need and submit job.
Open a terminal and load the module kilosort:
New software installed: kilosort version 4.0.13[New software installed: kilosort version 4.0.13](https://hpc-community.unige.ch/t/new-software-installed-kilosort-version-4-0-13/3551/1) HPC ChangeLog[HPC ChangeLog](https://hpc-community.unige.ch/c/hpc-announce/hpc-changelog/9)
> Dear users, we have installed a new software: kilosort 4.0.13: 

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  kilosort: kilosort/4.0.13
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------…
and execute `python -m kilosort` to open the soft


## Post 3 by @Njiva.Andrianarivelo (2024-09-24T13:06:54.829Z)

Adrien.Albert:
> select advanced desktop pick
This is wonderful thanks !
