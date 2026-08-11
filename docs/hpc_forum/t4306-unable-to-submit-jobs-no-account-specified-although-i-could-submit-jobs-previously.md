# Unable to submit jobs: "No account specified" although I could submit jobs previously

- Source: https://hpc-community.unige.ch/t/4306

- Created: 2026-06-09T18:30:09.929Z

- Tags: slurm

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Myriam.Trigilia (2026-06-09T18:30:09.995Z)

Hello,
I am currently unable to submit jobs on Baobab.
When running:
```
sbatch script_HC_simnibs_charm_array.sbatch
```
I receive:
```
sbatch: error: No account specified.
sbatch: error: Batch job submission failed: Invalid account or account/partition combination specified
```
I also tried:
- submitting a previous SimNIBS/CHARM script that worked successfully in the past,
- specifying an account explicitly:
```
sbatch --account=share_guedjc_teach script_HC_simnibs_charm_array.sbatch
```
- using another partition (`public-cpu` instead of `shared-cpu`).
In all cases, I obtained the same error.
I checked the cluster status using:
```
sinfo
```
and the partitions appear to be operational.
I also ran:
```
sshare -U -u $USER
sacctmgr show user $USER withassoc
```
Both commands returned only the table headers and no account information.
I have successfully submitted jobs on Baobab in the past, so I am wondering whether my Slurm account association may have changed?
Could someone please advise on how to check which account I should use, or whether my user account (`trigilia`) currently has a valid Slurm account association?
Thank you very much for your help.
Best regards,
Myriam Trigilia


## Post 2 by @Yann.Sagon (2026-06-10T08:30:28.809Z)

Dear @Myriam.Trigilia[@Myriam.Trigilia](https://hpc-community.unige.ch/u/myriam.trigilia) I’ve answered to your issue by email, I’m closing this thread on the forum.
