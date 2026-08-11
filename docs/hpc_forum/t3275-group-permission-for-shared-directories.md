# Group permission for shared directories

- Source: https://hpc-community.unige.ch/t/3275

- Created: 2024-01-26T12:31:30.248Z

- Posts: 5

- Category: 10

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Matthias.Kruckow (2024-01-26T12:31:30.308Z)

There is already an automatic script in place to ensure the permission in every body’s home directory to allow only the user to have permissions, see documentation[documentation](https://doc.eresearch.unige.ch/hpc/storage_on_hpc#home_directory).
Would it be possible to get a similar script to ensure group permissions in shared directories? Hence getting the group permission to be set to rwX?
Additionally, it would be nice if the correct group ownership could be ensured in the same way.
Currently we have a bash script to do this, but it needs to be run by any user separately and regularly to capture all new files. This is always an extensive job on the file system. Hence, we’d like to get a more automated and efficient solution for that.
To give an example on yggdrasil: all group members of `GL_S_Astro_POSYDON` should have access to all files in the full file structure below `/srv/beegfs/scratch/shares/astro/posydon`.
I guess, it will be the case for most stuff in shares, that one can assign a group to it.


## Post 2 by @Gael.Rossignol (2024-01-31T12:46:46.593Z)

Dear Matthias,
Actually we don’t want to manage this kind of scripts, I have change folder permission to set sticky bit, that permits all user to create new files and directory under group permission GL_S_Astro_POSYDON.
Is that fix can help you?
Best regards,


## Post 3 by @Matthias.Kruckow (2024-02-01T13:43:31.031Z)

The stick bit won’t help us.
OK, we’ll have to look for an own solution then.


## Post 4 by @Yann.Sagon (2024-02-01T15:16:29.032Z)

Hi,
my colleague @Gael.Rossignol[@Gael.Rossignol](https://hpc-community.unige.ch/u/Gael.Rossignol) set the sticky bit and did a test with one of the GL_S_Astro_POSYDON group member: he created a file and the later was member of the GL_S_Astro_POSYDON group. Is this not want you would like to achieve?
If this would be ok, we’ll set it like that next week. Right now, the sticky bit was removed by the script which enforce the permissions, thus we need to do some modification to it.


## Post 5 by @Matthias.Kruckow (2024-02-02T08:29:51.550Z)

Yes, it is one part to ensure the correct group.
The issue we have is, that sometimes files are created in a non-standard way, causing the group and permission for the group to be wrong. This we’d need to catch and correct for.
