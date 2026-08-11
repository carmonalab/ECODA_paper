# Scratch overburdened

- Source: https://hpc-community.unige.ch/t/3532

- Created: 2024-07-09T11:53:14.229Z

- Posts: 7

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Debajyoti.Sengupta (2024-07-09T11:53:14.273Z)

Dear HPC team,
I currently have a workflow that performs i/o on scratch on Baobab, but it looks like the scratch is quite slow at the moment. Even `ls` and `cd` take a while.
Can something be done regarding this?
Regards,
Deb


## Post 2 by @Adrien.Albert (2024-07-09T12:00:01.992Z)

Hi @Debajyoti.Sengupta[@Debajyoti.Sengupta](https://hpc-community.unige.ch/u/debajyoti.sengupta)
Could you please give us more details, on wich cluster do you have this issue ? which directories ?


## Post 3 by @Debajyoti.Sengupta (2024-07-09T12:02:39.333Z)

Hi Adrien,
Apologies for the lack of information.
This is on baobab. And the workflow needs access to `/srv/beegfs/scratch/groups/rodem/skycurtains/fullskyscan/`


## Post 4 by @Matthew.Leigh (2024-07-09T12:09:46.002Z)

I am finding the same on baobab. Anything on scratch, my own space or on /srv/beegfs/scratch/groups/rodem/ takes forever. As Deb said, just running ls takes a while.
I am trying to clear up a lot of space by deleting old runs, but the issue persisted before. This morning many of my jobs would fail after timing out trying to save simple log files.


## Post 5 by @Samuel.Klein (2024-07-09T13:47:25.319Z)

I am having the same issue, everything on scratch is slow


## Post 6 by @Adrien.Albert (2024-07-09T14:20:32.570Z)

Hi @Matthew.Leigh[@Matthew.Leigh](https://hpc-community.unige.ch/u/matthew.leigh)
I am looking for it but it’s quite difficult to find the root cause.


## Post 7 by @Yann.Sagon (2024-07-10T08:33:06.702Z)

Dear users,
I remind you that we have now a new cluster: Bamboo.
This cluster have a much faster storage! Home is fully on SSD.
Maybe a good time to give a try?!
