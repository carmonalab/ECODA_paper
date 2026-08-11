# Second user on GPU

- Source: https://hpc-community.unige.ch/t/3991

- Created: 2025-06-25T08:10:26.684Z

- Tags: baobab

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Manuel.Hentschel (2025-06-25T08:10:26.727Z)

## Primary informations
Username: hentsche
Cluster: baobab
Jobid: 837908
## Description
I requested a GPU and 80GB of VRAM and got allocated an `NVIDIA A100 80GB`.
There was another user’s process running on the GPU, causing my process to fail for not being able to allocate enough VRAM.
## Steps to Reproduce
Request a tunnel using the following sbatch command (linebreaks added here for legibility), and connect to it with vscode:
```
sbatch
  --job-name=codeTunnel
  --time=08:00:00
  --error=/home/users/h/hentsche/slurm/%j/stderr
  --output=/home/users/h/hentsche/slurm/%j/stdout
  --mail-type=ALL
  --mem=300000
  --cpus-per-task=32
  --gres=gpu:1,VramPerGpu:80G
  --partition=shared-gpu
  /home/users/h/hentsche/slurm_tunnel/tunnel.sh
```
In the created session, run a command that requires <80GB of VRAM.
## Expected Result
The command runs and has enough VRAM.
## Actual Result
My process failed because it ran out of VRAM at ca 62GB usage.
Below is a screen shot of the VRAM usage from while my process was running and one from after my process failed. The other process belongs to another user.
image
image1009×571 35.3 KB
[image1009×571 35.3 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/f767a27beedfaa5724f7ee36d6564791a3357ff8.png)
image
image1273×703 60.9 KB
[image1273×703 60.9 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/5dc2c9fe5378ac45871fa3adf4a3be0c90f9e1e4.png)


## Post 2 by @Yann.Sagon (2025-06-25T12:01:39.031Z)

Dear @Manuel.Hentschel[@Manuel.Hentschel](https://hpc-community.unige.ch/u/manuel.hentschel) thanks for the notification, this is now solved.
[2025] Current issues on HPC Cluster[[2025] Current issues on HPC Cluster](https://hpc-community.unige.ch/t/2025-current-issues-on-hpc-cluster/3788/15) HPC ChangeLog[HPC ChangeLog](https://hpc-community.unige.ch/c/hpc-announce/hpc-changelog/9)
> Description
Cluster: multi cluster 
Under specific circumstances, it is possible that you will be able to use a resource such as a GPU that is already associated with another user. This may happen when you request resources using salloc and then connect to the compute node using ssh later. Alternatively, Slurm may allocate you a resource that is already in use by a previous job. This issue has already been resolved, but it will only be fully resolved once all resources are released by the runnin…
