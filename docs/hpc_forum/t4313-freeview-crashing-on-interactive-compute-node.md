# Freeview crashing on interactive compute node

- Source: https://hpc-community.unige.ch/t/4313

- Created: 2026-06-12T15:26:34.750Z

- Posts: 2

- Category: 1

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Abigail.Licata (2026-06-12T15:26:34.817Z)

Hi HPC team,
I am trying to work with Freeview from FreeSurfer in an interactive session. It crashes immediately upon launching and I have not been able to figure out why. At first, I got the error :slight_smile:
No matching fbConfigs or visuals found glx: failed to create drisw screen Aborted (core dumped)
But now it is just crashing with message Aborted (core dumped).
This is my command line history:
ssh -Y myusername@login1.baobab.hpc.unige.ch
salloc -n1 -c2 --partition=public-interactive-cpu --time=00:15:00 --x11
module load FreeSurfer/7.3.2-centos7_x86_64
which freeview —> /opt/ebsofts/FreeSurfer/7.3.2-centos7_x86_64/bin/freeview
and when I just try to open freeview, I get this error.
Thank you in advance for any help you can provide!


## Post 2 by @Yann.Sagon (2026-06-26T08:24:10.465Z)

Dear @Abigail.Licata[@Abigail.Licata](https://hpc-community.unige.ch/u/abigail.licata)
The version you are trying is quite old, it was for CentOS7 and we are now using Rocky9.
Please give a try with the latest version we have on the cluster: `FreeSurfer/7.4.1-centos8_x86_64`.
As well, you may want to launch FreeView directly from OpenOnDemand[OpenOnDemand](https://doc.eresearch.unige.ch/hpc/how_to_use_openondemand) using the virtual desktop.
