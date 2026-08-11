# Time limit on private partition

- Source: https://hpc-community.unige.ch/t/4177

- Created: 2025-12-30T12:15:14.629Z

- Tags: baobab

- Posts: 7

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Maura.Brunetti (2025-12-30T12:15:14.677Z)

Dear HPC team,
I noticed that, despite I asked 7 days on my partition (private-gapnl-cpu), after 1 day the simulation stops `due to time limit’. Could you please check its setting?
Best wishes
Maura


## Post 2 by @Yann.Sagon (2026-01-05T14:32:09.300Z)

Dear @Maura.Brunetti[@Maura.Brunetti](https://hpc-community.unige.ch/u/maura.brunetti),
are you using sbatch or salloc? Can you please share a job id which has the issue?
Best regards
Yann


## Post 3 by @Laure.Moinat (2026-01-13T12:32:39.367Z)

Dear Yann,
I also have the same problem with the private-gapnl-cpu partition, where the simulations stop after a couple of hours without any reason. When I launch the same simulation, but with the public-cpu partition, it runs for 4 days.
Do you have an explanation for this behavior?
Best,
Laure


## Post 4 by @Yann.Sagon (2026-01-13T16:36:20.498Z)

Dear @Laure.Moinat[@Laure.Moinat](https://hpc-community.unige.ch/u/laure.moinat) I have the same question for you as I did for Maura:)


## Post 5 by @Laure.Moinat (2026-01-14T10:29:45.066Z)

Dear Yann,
This is how we launch our simulations ‘srun --ntasks=25 --partition=public-cpu --time=4-00:00:00 --multi-prog P280.conf > std_outp 2>&1’. I currently have no job ID that is running with this issue.  I can launch one if necessary.
Best,
Laure


## Post 6 by @Yann.Sagon (2026-01-14T10:52:24.079Z)

Laure.Moinat:
> This is how we launch our simulations
Thanks for the feedback. Unless used for a good reason, do not use srun to launch long running tasks. If you disconnect from the login node the job is killed. If the login node has an issue and we restart it (this happened multi time since beginning of this year) the job is killed also.
Create an sbatch script and launch it using sbatch `<your script>`.
Example to create `<your script>`
```
#SBATCH --ntasks=25 
#SBATCH --partition=public-cpu 
#SBATCH --time=4-00:00:00 
#SBATCH --multi-prog

srun P280.conf > std_outp 2>&1’.
```


## Post 7 by @Laure.Moinat (2026-01-14T11:00:54.361Z)

Thanks for the clarification and for the example, I will change my scripts!
Best,
Laure
