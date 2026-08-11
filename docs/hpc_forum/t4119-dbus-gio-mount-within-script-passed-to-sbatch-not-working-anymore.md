# Dbus + gio mount within script passed to sbatch not working anymore

- Source: https://hpc-community.unige.ch/t/4119

- Created: 2025-10-10T14:45:19.950Z

- Tags: bamboo

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Aurelie.Kamoun (2025-10-10T14:45:20.026Z)

## Primary informations
Username: kamouna
Cluster: bamboo
## Description
The mounting of my NAS share does not work anymore when i want to include it in a script that i launch with sbatch on bamboo. However if i run the dbus-launch + gio mount lines one after the other on an salloc session it works (if i copy paste the 2 lines in one go it fails). I tried to add sleep 5 after dbus-launch bash and before gio mount but it did not solve the issue. The exact same script used to work before (2 weeks ago).
These are the lines in the script that i pass to sbatch
`module purge`
`module load GCCcore/12.3.0 Python`
`dbus-launch bash`
`gio mount smb://nasac-m2.unige.ch/m-CDGG < .credentials`
When i run the sbatch command calling this script, i get the following error in the logfile:
```
Error creating proxy: Could not connect: Connection refused (g-io-error-quark, 39)
gio: smb://nasac-m2.unige.ch/m-CDGG: volume doesn’t implement mount
```
## Steps to Reproduce
I tried to debug using salloc and copy pasted these 4 lines in the shell
`module purge`
`module load GCCcore/12.3.0 Python`
`dbus-launch bash`
`gio mount smb://nasac-m2.unige.ch/m-CDGG < .credentials`
→ The mount was not successful
However when i sequentially paste the 3 first lines:
`module purge`
`module load GCCcore/12.3.0 Python`
`dbus-launch bash`
Then the last one:
`gio mount smb://nasac-m2.unige.ch/m-CDGG < .credentials`
→ The mount is successful
## Expected Result
I expect that the mount to be successful after
`dbus-launch bash `
`gio mount smb://nasac-m2.unige.ch/m-CDGG < .credentials `
when they are run via sbatch
## Actual Result
I get the following error
```
Error creating proxy: Could not connect: Connection refused (g-io-error-quark, 39)
gio: smb://nasac-m2.unige.ch/m-CDGG: volume doesn’t implement mount
```


## Post 2 by @Adrien.Albert (2025-10-15T08:07:13.683Z)

Hello @Aurelie.Kamoun[@Aurelie.Kamoun](https://hpc-community.unige.ch/u/aurelie.kamoun)
You need to start with a clean session. Before running `salloc`, make sure you do not have any GIO session running. When using `sbatch`, you should add a `sleep` between `dbus-launch` and the GIO command. We suspect that GIO is not fully initialized before the GIO command executes, which is why the delay is necessary.
In your `sbatch` script:
```
dbus-launch bash 
sleep 5
gio mount smb://nasac-m2.unige.ch/m-CDGG < .credentials
```
Using slurm on a NAS share leads to slurmstepd: error: couldn't chdir to?[Using slurm on a NAS share leads to slurmstepd: error: couldn't chdir to?](https://hpc-community.unige.ch/t/using-slurm-on-a-nas-share-leads-to-slurmstepd-error-couldnt-chdir-to/3960/6) HPC issues[HPC issues](https://hpc-community.unige.ch/c/hpc-support/hpc-issues/14)
> Hi @Matthieu.Stigler[@Matthieu.Stigler](https://hpc-community.unige.ch/u/matthieu.stigler) 
Here are my tests. The script in the SMB share contains only: 
echo toto

Reproducing the current situation (Not working):
(baobab)-[alberta@login1 ~]$ cat !$
cat sbatch_test.sh
#!/bin/sh
#SBATCH --job-name test_gio
#SBATCH --cpus-per-task 1
#SBATCH --time 00:05:00
#SBATCH --partition debug-cpu

dbus-launch bash
gio  mount smb://isis.unige.ch/nasac/hpc_exchange/backup < .credentials
bash /var/run/user/401775/gvfs/smb-share\:server\=isis.unige.ch\,share\=nasac/hpc_exchange/b…
