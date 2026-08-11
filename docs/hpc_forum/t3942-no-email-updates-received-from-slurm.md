# No email updates received from SLURM

- Source: https://hpc-community.unige.ch/t/3942

- Created: 2025-04-29T09:17:12.039Z

- Tags: baobab, slurm

- Posts: 7

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Marco.Froelich (2025-04-29T09:17:12.087Z)

Username: froelicm
Cluster: Baobab
I am not recieving any email notification from jobs submitted to Baobab GPUs, although the jobs are running and terminate error-free. Last week this was working: since yesterday, I am had no more updates.
This is the relevant SBATCH instruction in my script being submitted (unchanged since last week):
#SBATCH --mail-user=marco.froelich@unige.ch
#SBATCH --mail-type=ALL


## Post 2 by @Camille.Serquet (2025-04-29T12:13:45.446Z)

Both me and my colleague have the same issue !


## Post 3 by @Adrien.Albert (2025-04-30T07:39:22.950Z)

Hi @Marco.Froelich[@Marco.Froelich](https://hpc-community.unige.ch/u/marco.froelich),
I take a look, I keep you inform


## Post 4 by @nicolas.clairis (2025-04-30T07:43:42.810Z)

Same issue here since monday, still not solved apparently


## Post 5 by @Adrien.Albert (2025-04-30T08:16:17.343Z)

@nicolas.clairis[@nicolas.clairis](https://hpc-community.unige.ch/u/nicolas.clairis), Yes still on it :slight_smile:


## Post 6 by @nicolas.clairis (2025-04-30T08:18:39.137Z)

All the pending emails arrived at 10:11 so I am guessing this is fixed now


## Post 7 by @Adrien.Albert (2025-04-30T08:30:49.526Z)

The issue has been fixed
[2025] Current issues on HPC Cluster[[2025] Current issues on HPC Cluster](https://hpc-community.unige.ch/t/2025-current-issues-on-hpc-cluster/3788/12) HPC ChangeLog[HPC ChangeLog](https://hpc-community.unige.ch/c/hpc-announce/hpc-changelog/9)
> Description
Cluster: baobab 
An unexpected behavior temporarily impacted the automatic sending of outgoing emails on one of our servers. The issue has been resolved, and email delivery is functioning normally. No email loss has been detected during the incident. 
Status : solved green_circle
start: 2025-04-29T09:00:00Z (UTC) 
end : 2025-04-30T08:21:00Z (UTC)
