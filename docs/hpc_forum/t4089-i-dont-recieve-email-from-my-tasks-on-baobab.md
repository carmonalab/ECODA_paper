# I don't recieve email from my tasks on baobab

- Source: https://hpc-community.unige.ch/t/4089

- Created: 2025-09-17T07:50:31.856Z

- Tags: baobab

- Posts: 2

- Category: 13

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Leonhard.Unterlerchner (2025-09-17T07:50:31.915Z)

Dear hpc community,
Since the 5th of september I didn’t revieve any email fomr my tasks (on baobab). I have not changed anything in my batch files, here the arguments:
`#SBATCH` `--mail-user=myEmailAdress`
`#SBATCH` `--mail-type=BEGIN`
`#SBATCH` `--mail-type=END`
`#SBATCH` `--mail-type=FAIL`
Two questions:
- Am I the only one with this issue?
- Can I do something to fix this?
Many thanks, best regards,
Leonhard.


## Post 2 by @Adrien.Albert (2025-09-17T13:22:45.974Z)

Hi @Leonhard.Unterlerchner[@Leonhard.Unterlerchner](https://hpc-community.unige.ch/u/leonhard.unterlerchner)
- Are you sure MyEmailAdress is your correct email.
- is `$ getent aliases $USER` returns your email adress (no syntax errors ?)
I have tested as user:
```
(baobab)-[alberta@login1 ~]$ srun --mail-type=BEGIN,END  hostname
cpu001.baobab
```
I receive all emails
