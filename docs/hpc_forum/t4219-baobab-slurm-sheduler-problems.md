# Baobab slurm sheduler problems

- Source: https://hpc-community.unige.ch/t/4219

- Created: 2026-02-10T14:01:01.849Z

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Ivan.Oleksiyuk (2026-02-10T14:01:01.908Z)

## Primary informations
Username: oleksiyu
Cluster: baobab
## Description
Something seems to be not right with the SLURM scheduler on baobab:
- All jobs that I submit start and the “timer” goes, but there is no output of any job, even the slurm log file is not created.
- The old cancelled jobs hang in the list indefinitely
- When I try to ssh to a node of an “running” job I get this:
- Access denied by pam_slurm_adopt: you have no active jobs on this node
- Connection closed by 192.168.103.8 port 22
Could this be fixed?
image
image1225×288 19.5 KB
[image1225×288 19.5 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/fcb8f0a9c371a8aa7905c32dfd7ca23973bc1c68.png)
## Steps to Reproduce
Submit any small job on Baobab with some log output path (e.g. below) and cancel it after some time.
#SBATCH --output=/home/users/o/oleksiyu/output.out
## Expected Result
What did you expect to happen when running the steps above?
I expect the log to appear at a specified location and after cancelling the job to be gone from “squeue --me“ in less than 1 min
## Actual Result
No log appears. Job still hangs in the queue after cancelling for a long time


## Post 2 by @Ivan.Oleksiyuk (2026-02-10T14:11:38.156Z)

Seems to be working again for now (15:11 10.02.2026)
