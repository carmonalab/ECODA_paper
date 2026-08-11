# Communication error on send on Baobab's /srv/beegfs/dpnc

- Source: https://hpc-community.unige.ch/t/4159

- Created: 2025-12-07T00:32:14.172Z

- Tags: baobab

- Posts: 4

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Hugo.Boutin (2025-12-07T00:32:14.225Z)

## Primary informations
Username: boutinh
Cluster: Baobab
Beegfs server: /srv/beegfs/dpnc/ on specific files
## Description
Commands trying to access files: cat, vim, md5sum, rsync or zip locally or trying to move files to another server will either freeze (`ctrl+c` not working) or return “Communication error on send” and then freeze. htop shows the command with status S. `ctrl+c` changes that status to D.
This appears only on the beegfs server mentioned above and for certain files (I tried other directories and files and the commands execute). It appears the files I have issues with were all moved there/generated in the past month.
This issue appeared mid-afternoon yesterday (Sat 6th).
## Steps to Reproduce
from `/srv/beegfs/dpnc/groups/dampe/users/hugo/libeb_flux/MC_select/`:
`md5sum Filtered_Etruth_Ereco_Weight_PSDCharge_PSDCut_PSDRawCharge_PSDChargeOld_PSDChargeBoron_STKHits_noSTKFilter_Proton.npy ~`
## Expected Result
For the checksum to print
## Actual Result
`md5sum Filtered_Ereco_PSDCharge_PSDCut_PSDChargeOld_PSDChargeBoron_STKHits_20*`
`md5sum: Filtered_Ereco_PSDCharge_PSDCut_PSDChargeOld_PSDChargeBoron_STKHits_2015.npy: Communication error on send`


## Post 2 by @Adrien.Albert (2025-12-08T10:44:18.763Z)

Hello @Hugo.Boutin[@Hugo.Boutin](https://hpc-community.unige.ch/u/hugo.boutin)
Do you still have the issue ?
On my side:
```
(baobab)-[root@login1 MC_select]$  md5sum Filtered_Etruth_Ereco_Weight_PSDCharge_PSDCut_PSDRawCharge_PSDChargeOld_PSDChargeBoron_STKHits_noSTKFilter_Proton.npy
9e1ae3ebf62ba3d3a319a714b7f2ff18  Filtered_Etruth_Ereco_Weight_PSDCharge_PSDCut_PSDRawCharge_PSDChargeOld_PSDChargeBoron_STKHits_noSTKFilter_Proton.npy
```


## Post 3 by @Hugo.Boutin (2025-12-08T12:45:57.076Z)

Hi Adrien,
Indeed it is working now, thanks! What was the issue?


## Post 4 by @Adrien.Albert (2025-12-08T13:25:40.176Z)

I have too few information to investigate more :confused:
