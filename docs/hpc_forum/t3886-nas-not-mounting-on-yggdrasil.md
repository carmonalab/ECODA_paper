# NAS not mounting on yggdrasil

- Source: https://hpc-community.unige.ch/t/3886

- Created: 2025-03-26T15:53:41.858Z

- Tags: yggdrasil

- Posts: 2

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Roberto.Sierra (2025-03-26T15:53:41.906Z)

Dear HPC team,
When trying to mount the NAS on yggdrasil it seems as if it mounts but to an unknown or unretrievable path.
The NAS appeared in `~/.gvfs`, but no longer the case and in the alternative path `/run/user/($(u -id)/`, the user ID folder is not created.
```
(yggdrasil)-[sierrami@login1 ~]$ dbus-launch bash
(yggdrasil)-[sierrami@login1 ~]$ gio mount 'smb://ISIS;sierrami@nasac-m2.unige.ch/m-Andrey_KPC'
Authentication Required
Enter password for share “m-andrey_kpc” on “nasac-m2.unige.ch”:
Password: 
(yggdrasil)-[sierrami@login1 ~]$ gio mount 'smb://ISIS;sierrami@nasac-m2.unige.ch/m-Andrey_KPC'
gio: smb://ISIS;sierrami@nasac-m2.unige.ch/m-andrey_kpc/: Location is already mounted
```
Is there a new location/path?
Thank you,
Roberto


## Post 2 by @Gael.Rossignol (2025-03-31T14:35:17.722Z)

Dear Roberto,
Actually on yggdrasil it’s waiting for next maintenance to be fixed, but there are possibilities to access files using gio.
Documentation is available :
doc.eresearch.unige.ch[doc.eresearch.unige.ch](https://doc.eresearch.unige.ch/hpc/storage_on_hpc#sometimes_mount_is_not_available_but_you_can_browsecopyinterract_with_gio_commands)
### hpc:storage_on_hpc [eResearch Doc][hpc:storage_on_hpc [eResearch Doc]](https://doc.eresearch.unige.ch/hpc/storage_on_hpc#sometimes_mount_is_not_available_but_you_can_browsecopyinterract_with_gio_commands)
Best regards,
