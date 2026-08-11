# [Baobab] Is scratch down?

- Source: https://hpc-community.unige.ch/t/3274

- Created: 2024-01-26T09:59:37.866Z

- Posts: 4

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @maciej.falkiewicz (2024-01-26T09:59:37.921Z)

what did you try: Copy to / list on scratch
what didn’t work: none of them
what was the expected result: successful copying / listing of files
what was the error message: timeout / frozen terminal / `Communication error on send`
path to the relevant files (logs, sbatch script, etc): /home/users/f/falkiewi/scratch


## Post 2 by @Alexander.Couzens (2024-01-26T10:09:14.058Z)

I also can’t SSH to Baobab (or access the scratch)
Could not establish connection to “login2.baobab.hpc.unige.ch”: Connecting with SSH timed out.


## Post 3 by @Gael.Rossignol (2024-01-26T10:22:38.228Z)

Dear users,
Login2 is overloaded that is explains connection problems. We are trying to identify the process but it seem this is related of the global usage of login2.
Does users using rsync command can be relaunch it from a node?
Best regards,


## Post 4 by @Gael.Rossignol (2024-01-26T10:34:29.185Z)

Dear all,
We have found the issue, one beegfs storage services was down because of too many open file at the same time.
This process has been restarted and all is now working fine. We will debrief to check the origin of the issue.
Best regards,
