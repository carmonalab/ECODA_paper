# Python not available on all the nodes

- Source: https://hpc-community.unige.ch/t/3800

- Created: 2025-01-24T15:14:19.991Z

- Posts: 6

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Kinga.Wozniak (2025-01-24T15:00:21.964Z)

hello,
I have a very similar problem, but without using venv:
I `module load GCCcore/12.3.0 Python/3.11.3` for which I have pip installed snakemake together with the slurm-adaptor.
when launching jobs I get the error
```
/var/spool/slurmd/job14437715/slurm_script: line 4: /opt/ebsofts/Python/3.11.3-GCCcore-12.3.0/bin/python: No such file or directory
```
on a non-deterministic basis (depending on the time 20-100% of jobs failing)
thank you


## Post 2 by @Kinga.Wozniak (2025-01-24T15:14:19.991Z)

adding that since a couple of hours all of my jobs fail with that error message


## Post 3 by @Malte.Algren (2025-01-24T15:28:30.925Z)

I have had the same error starting this week
Cheers,
Malte


## Post 4 by @Adrien.Albert (2025-01-24T16:10:44.227Z)

On which cluster you have this error ?


## Post 5 by @Malte.Algren (2025-01-27T07:58:23.601Z)

Hi again
Its on baobab
Malte


## Post 6 by @Yann.Sagon (2025-01-28T10:08:59.312Z)

Hi, I’ve split the original post at this are two distinct problems I guess.
In your case, the software that isn’t seen is `/opt/ebsofts/Python/3.11.3-GCCcore-12.3.0/bin/python` which is hosted on a central NFS share, thus not related with your home directory.
I’ve checked if this location is reachable from all the nodes… and no! It wasn’t reachable from cpu333 which was the node you used! This is now fixed.
Best


## Post 7 by @Yann.Sagon (2025-01-28T10:10:11.982Z)
