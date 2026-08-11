# CPU/GPU job failing

- Source: https://hpc-community.unige.ch/t/3432

- Created: 2024-04-30T12:09:33.207Z

- Posts: 4

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Debajyoti.Sengupta (2024-04-30T12:09:33.263Z)

Hello all, I am trying to get a cpu node to mount my vs code session on - for which I use this slurm script.
```
#!/bin/sh
#SBATCH --job-name=bgcpu
#SBATCH --cpus-per-task=1
#SBATCH --time=00-12:00:00
#SBATCH --partition=private-dpnc-cpu,shared-cpu
#SBATCH --output=/home/users/s/senguptd/logs/slurm-%A.out
#SBATCH --mem=25GB

export XDG_RUNTIME_DIR=""
echo $HOSTNAME
srun sleep 12h
```
Once I get a node I login to that node and mount my VSCode for the day.
However I am getting this error when I try to launch the job
image
image990×704 27.8 KB
[image990×704 27.8 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/bac8898eeab98d7928038e6c2129ab47409a4cf2.png)
Any idea what the raisedsignal:53 refers to, and what the solution to this is?


## Post 2 by @Debajyoti.Sengupta (2024-05-01T10:45:30.299Z)

It looks like none of my scripts are able to launch jobs - all of them fail with the same error: RaiseSignal:53 (they don’t seem to create any log file either).
Any help is appreciated as my entire workflow is currently disrupted.


## Post 3 by @Debajyoti.Sengupta (2024-05-01T11:03:41.051Z)

Okay I did some digging and it looks like it was because I had hit disk quota. I cleared my tmp (all I could) and it looks like things are running for now.


## Post 4 by @Yann.Sagon (2024-05-02T13:33:43.717Z)

Debajyoti.Sengupta:
> Okay I did some digging and it looks like it was because I had hit disk quota
Dear @Debajyoti.Sengupta[@Debajyoti.Sengupta](https://hpc-community.unige.ch/u/debajyoti.sengupta) , my2cents: you can check your quota usage with the command listed [here].(hpc:storage_on_hpc [eResearch Doc][hpc:storage_on_hpc [eResearch Doc]](https://doc.eresearch.unige.ch/hpc/storage_on_hpc#check_disk_usage_on_home_and_scratch))
