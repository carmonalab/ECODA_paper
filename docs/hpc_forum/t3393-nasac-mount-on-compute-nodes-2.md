# NASAC mount on compute nodes 2

- Source: https://hpc-community.unige.ch/t/3393

- Created: 2024-03-26T08:29:42.463Z

- Posts: 2

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Nada.Kojovic (2024-03-26T08:29:42.525Z)

Hi All
I am desperately trying to start processing by mounting nasac on the compute nodes but is been weeks and I am very lost cause none of the suggested solutions worked.
I have tried all the steps: hpc:storage_on_hpc [eResearch Doc][hpc:storage_on_hpc [eResearch Doc]](https://doc.eresearch.unige.ch/hpc/storage_on_hpc#access_external_storage)
Nothing happens! I would like to have a brief in person exchange so we understand where is the problem stemming from, is it my permissions or something else
MANY thanks !!
Nada


## Post 2 by @Adrien.Albert (2024-03-26T17:54:59.943Z)

The SMB mount point was not functional due to the presence of a directory under the destination. Example
I want to mount the directory bio_med_info from NasAc but the directory already exists in the mount destination `(~/.gvfs`). The command will run without an error message, but the data will not be from NasAc but from the underlying directory.
The root cause is unknown
To remedy this:
- Be sure to unmount the SMB (`gio mount -u [...]`) (to avoid deletion on NasAc)
- Delete the underlying directory (the data being on the cluster it will not be deleted on the NasAc since it has been unmounted)
- Kill all your processes (`pkill -u $USER`) unlog then relog on cluster and start the procedure the SMB from scratch.
This procedure has been tested and approuved with @Nada.Kojovic[@Nada.Kojovic](https://hpc-community.unige.ch/u/nada.kojovic)


## Post 3 by @Yann.Sagon (2024-05-02T14:02:43.361Z)

A post was split to a new topic: Unable to mount NASAC share on[Unable to mount NASAC share on](https://hpc-community.unige.ch/t/unable-to-mount-nasac-share-on/3435)
