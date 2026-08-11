# Unable to access '/srv/beegfs/dpnc/groups/atlas': Permission denied

- Source: https://hpc-community.unige.ch/t/4140

- Created: 2025-11-07T09:03:55.435Z

- Tags: baobab

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Vilius.Cepaitis (2025-11-07T09:03:55.494Z)

Dear HPC team,
Since this morning, I’m no longer able to access the following directory:
`ls /srv/beegfs/dpnc/groups/atlas   ls: cannot open directory`
Could you please take a look into this?
Cheers,
Vilius


## Post 2 by @Gael.Rossignol (2025-11-10T12:51:48.923Z)

Dear Vilius,
Acl was not correctly set, I modify ACL and now you may access to the share.
Best regards,
