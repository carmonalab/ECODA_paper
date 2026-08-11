# Slurm scheduling: dependency and requeueing

- Source: https://hpc-community.unige.ch/t/3587

- Created: 2024-08-07T15:57:06.492Z

- Tags: all

- Posts: 4

- Category: 10

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Matthias.Kruckow (2024-08-07T15:57:06.529Z)

# When does SLURM set the date to count the age for the scheduler?
Simply it sets it, whenever a job gets the OK, for NODELIST(REASON)=(Priority).
# But is that the way it should be?
Here two examples:
- A node has a failure and a job gets requeued.
- A job is dependent on other job(s).
In both cases, the age will start at the moment the job can run and not at the time it was submitted.
# What are the consequences?
Some jobs get delayed behind others, because the ageing starts not on submission.
Here an example, when submitting jobs J1 - J100 with depending jobs D1to10 - D91to100. Even you submit J1 - J10, D1to10, J11 - 20, D11to20, … SLURM will schedule them J1 - J100, D1to10 - D91to100. Hence, you need to wait until your last submitted job J100 finished before the jobs D* start and therefore you won’t get results on an already finished part, because a dependent post-processing is still pending to late in the queue. (I concentrate in the details on this example because that is reproduceable with a single account)
In the case of a requeue, it is even worse, that you can even not do anything about getting suddenly a run being done later, even this might be a more important one among yours (but SLURM decides to submit one of your others in the queue before doing the requeued one). Sure you can cancel and resubmit your other jobs, which will cause a waste of your and CPU time. Additionally, it won’t help when other users submitted a lot of jobs, too (Hence, a very limited workaround).
# A workaround, one shouldn’t use.
For the dependent jobs one could get around wrong age, by not using the dependency option of SLURM, but instead put it as a job, which till wait/sleep until the others finish, but can run directly after, because it is already running. This would block CPUs with code just waiting, hence I recommend no one to do that.
# What I’d like to have changed?
Is it possible to configure SLURM such that the ageing always uses the difference between submission and now, instead of ready to run and now?
(This would help stopping the need for people using dirty workarounds as stated before.)


## Post 2 by @Adrien.Albert (2024-08-07T16:50:51.091Z)

Hi @Matthias.Kruckow[@Matthias.Kruckow](https://hpc-community.unige.ch/u/matthias.kruckow),
I am having trouble understanding everything. Could you please send me your sbatch script?
I need to see how you are managing your job dependencies.
If you run a job array with IDs [0-100] (YYYYY) and you request a dependency on YYYYY for job XXXXX, it is normal for job XXXXX to wait until the entire job array finishes before starting. However, if you request a dependency on specific job array IDs (e.g., YYYYY_0, YYYYY_1, etc.), job XXXXX should start after each corresponding job array element finishes, not after the entire job array.
But it seems difficult to create a dependency of JobArrayID_Y on JobArrayID_X.


## Post 3 by @Matthias.Kruckow (2024-08-09T14:56:36.718Z)

I have sent you a minimum working example via email.


## Post 4 by @Matthias.Kruckow (2024-08-28T14:17:34.092Z)

Just to summarize the email conversation:
Adjusting the SLURM setting ACCRUE_ALWAYS, might be a solution to keep the finishing order of jobs as close as possible to the submission order, when dependencies or requeues are involved.
As far as @Gael.Rossignol[@Gael.Rossignol](https://hpc-community.unige.ch/u/Gael.Rossignol) told me, making this change would need to be discussed internally for the HPC team, because it might has other side effects.
