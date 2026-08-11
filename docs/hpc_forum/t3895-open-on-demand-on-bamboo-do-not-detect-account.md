# Open-on-demand on bamboo do not detect account

- Source: https://hpc-community.unige.ch/t/3895

- Created: 2025-03-31T10:55:27.169Z

- Tags: bamboo, openondemand

- Posts: 3

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Lucille.Delisle1 (2025-03-31T10:55:27.221Z)

Dear HPC team,
On the open-on-demand of bamboo all the apps that propose the field ‘auto_accounts’ do not propose any option while on yggdrazil and baobab it works fine…
On baobab:
image
image586×230 11.6 KB
[image586×230 11.6 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/d7b178fee1f4b3a797f6f8e630993c1825f56cc1.png)
On bamboo:
image
image604×167 6.56 KB
[image604×167 6.56 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/e1198e93c46a1c401d473a4bdb4cbc2cc3edf717.png)


## Post 2 by @Yann.Sagon (2025-03-31T14:23:49.551Z)

Dear @Lucille.Delisle1[@Lucille.Delisle1](https://hpc-community.unige.ch/u/lucille.delisle1)
Lucille.Delisle1:
> On the open-on-demand of bamboo all the apps that propose the field ‘auto_accounts’ do not propose any option
Thanks for the notification, this is solved by @Gael.Rossignol[@Gael.Rossignol](https://hpc-community.unige.ch/u/Gael.Rossignol)
Best


## Post 3 by @Rut.Gabarrosolanas (2026-03-26T16:27:23.222Z)

Hello, I also can’t indicate the account on open-on-demand bamboo and get this error message:
Failed to submit session with the following error:
```
sbatch: error: Batch job submission failed: Invalid account or account/partition combination specified
```
- If this job failed to submit because of an invalid job name please ask your administrator to configure OnDemand to set the environment variable OOD_JOB_NAME_ILLEGAL_CHARS.
- The RStudio session data for this session can be accessed under the staged root directory[staged root directory](https://openondemand.bamboo.hpc.unige.ch/pun/sys/dashboard/files/fs/home/users/g/gabarros/ondemand/data/sys/dashboard/batch_connect/sys/ug_RStudio/output/477fea76-70cd-4783-a952-68d66fb41db4).
Could you please help? Thanks!
Best,
Rut
