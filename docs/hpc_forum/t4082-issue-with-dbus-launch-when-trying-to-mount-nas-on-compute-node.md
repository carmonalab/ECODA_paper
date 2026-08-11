# Issue with dbus-launch when trying to mount NAS on compute node

- Source: https://hpc-community.unige.ch/t/4082

- Created: 2025-09-05T09:34:00.339Z

- Tags: bamboo

- Posts: 4

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Aurelie.Kamoun (2025-09-05T09:34:00.400Z)

## Primary informations
Username: kamouna
Cluster: bamboo
## Description
I would like to mount a nas on a compute node but i get an error message after running dbus-launch and i do not know how to address it.
Note : The nas is already mounted on my login node and i have no issue with it.
## Steps to Reproduce
Trying
```
salloc
module load GCCcore/12.3.0 Python
dbus-launch bash
```
(I load the above module because if i do not do this i get the following error : /opt/ebsofts/Python/3.11.3-GCCcore-12.3.0/bin/python: error while loading shared libraries: libpython3.11.so.1.0: cannot open shared object file: No such file or directory)
## Expected Result
I expect to then be able to run the following command to get my nas accessible
```
gio mount smb://nasac-m2.unige.ch/m-CDGG < ~/.credentials
```
## Actual Result
This is the error i get
```
flatpak: /opt/ebsofts/util-linux/2.39-GCCcore-12.3.0/lib/libmount.so.1: version `MOUNT_2_40' not found (required by /lib64/libgio-2.0.so.0)
```


## Post 2 by @Adrien.Albert (2025-09-05T10:54:49.321Z)

Hi @Aurelie.Kamoun[@Aurelie.Kamoun](https://hpc-community.unige.ch/u/aurelie.kamoun)
Be carefull, you need to use the gio from the system and not from a module.
Please deactivated conda environment and all module before running dbus and gio.


## Post 3 by @Aurelie.Kamoun (2025-09-08T11:42:34.075Z)

Hello
If I do
```
salloc

module purge
dbus-launch bash
```
I get
```
/opt/ebsofts/Python/3.11.3-GCCcore-12.3.0/bin/python: error while loading shared libraries: libpython3.11.so.1.0: cannot open shared object file: No such file or directory
```


## Post 4 by @Yann.Sagon (2025-09-15T09:14:14.925Z)

Dear @Aurelie.Kamoun[@Aurelie.Kamoun](https://hpc-community.unige.ch/u/aurelie.kamoun)
I’ve checked what can be cause the error message `error while loading shared libraries: libpython3.11.so.1.0...` when you do an salloc.
In the path `.local/bin` you have a binary named `register-python-argcomplete` that seems to come from the EasyBuild Python version `/opt/ebsofts/Python/3.11.3-GCCcore-12.3.0/`. This script is called by `/etc/bash_completion.d/tracer` when you create a new session for example when using `salloc`.
What you can do:
- remove the unneeded stuff in `.local/bin`
- use a version of this script that don’t needs this specific Python version
