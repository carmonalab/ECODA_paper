# Jobs taking much longer than before/failing before they are even ran

- Source: https://hpc-community.unige.ch/t/3647

- Created: 2024-09-24T09:59:56.239Z

- Tags: yggdrasil

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @William.Ceva (2024-09-24T09:59:56.287Z)

Hello,
I am running jobs, that are (I am almost certain) practically identical to jobs I ran ~2 weeks ago.  However, these same jobs are now taking twice or three times as long to run.  I am submitting to the following partitions:
private-astro-cpu, public-cpu, shared-cpu, public-bigmem, shared-bigmem
I have also had (beginning this morning) jobs fail immediately after submitting them.  No output or error files are produced, they just fail immediately once they are submitted with no way of determining what happened.  Examples (from the emails I did receive notifying me the jobs failed):
Slurm Job_id=35559853 Name=HD8049_2A Failed, Run time 00:00:00, FAILED
Slurm Job_id=35559839 Name=HD7449_16B Failed, Run time 00:00:00, FAILED
Slurm Job_id=35559833 Name=HD7449_16A Failed, Run time 00:00:00, FAILED
~40 jobs had this issue.
I am wondering if these issues, particularly my jobs being much slower now, are related to the recent electrical issues on yggdrasil: [2024] Current issues on HPC Cluster - #23 by Gael.Rossignol[[2024] Current issues on HPC Cluster - #23 by Gael.Rossignol](https://hpc-community.unige.ch//hpc-community.unige.ch/t/2024-current-issues-on-hpc-cluster/3245/23)?
Is it possible that the cpus on yggdrasil are now, unknowingly/unintentionally after the fixes to the electrical system, being “undervolted” and are therefore slower?  I am naive about hpc in general, so apologies if “undervolting” is irrelevant here.


## Post 2 by @Adrien.Albert (2024-09-24T16:00:45.059Z)

Hi @William.Ceva[@William.Ceva](https://hpc-community.unige.ch/u/william.ceva)
I’ve checked and the node where your job was run doesn’t have the beegfs Home and scratch mounted. I’ve drained them from production for analysis.
