# No start time listed in squeue

- Source: https://hpc-community.unige.ch/t/3768

- Created: 2024-12-17T09:10:45.716Z

- Posts: 4

- Category: 13

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Cody.Cardenas (2024-12-17T09:10:45.760Z)

The past few weeks I have noticed that start my start times are almost nonexistent and might not show up for over a day. I usually use the start time to consider whether or not I will start another batch of jobs immediately or work on something else in the mean time. For example, my latest array job:
```
squeue --format="%.18i %.9P %.40j %S %.8T %.10M %.9l %.6D %R" -u $USER
             JOBID PARTITION                                     NAME START_TIME    STATE       TIME TIME_LIMI  NODES NODELIST(REASON)
   14069541_[1-16] public-cp                              BEAST_rerun N/A  PENDING       0:00   4:00:00      1 (Priority)
```
I am curious if there is something I am doing that is preventing my jobs from running? Or is it just that my priority is so low that the job scheduler is uncertain of the start time?
I know I’ve asked for and submitted quite a few jobs lately, so I wonder if the later is the case.
Thanks in advance!
(PS - I know asking for a full four days is not encouraged, but these beast runs can take quite a while to run)


## Post 2 by @Yann.Sagon (2024-12-18T15:36:31.171Z)

Dear @Cody.Cardenas[@Cody.Cardenas](https://hpc-community.unige.ch/u/cody.cardenas)
Cody.Cardenas:
> Or is it just that my priority is so low that the job scheduler is unsure of the start time?
That is probably the reason. In fact, it is not your priority that is low, it is that you are using the public-cpu partition, which has only 8 compute nodes.
Cody.Cardenas:
> I know asking for a full four days is not encouraged, but these beast runs can take quite a while to run.
Is it possible to use some checkpointing (restart) perhaps?
Another suggestion is to use Bamboo: 29 compute nodes (128 cores each) are sitting idle waiting for your job. Bonus, the home storage is on SSD and super fast.


## Post 3 by @Cody.Cardenas (2024-12-19T08:56:37.222Z)

Ive not had access to bamboo since it was live, I’ll look into how to access it.
I didn’t think that Beast had checkpointing, but I’ll have a look into it and reconsider how to format my jobs so I can meet the 100-200 million MCMC chain iterations necessary with less demanding times.
Thanks for the feedback @Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon)!


## Post 4 by @Yann.Sagon (2024-12-19T09:28:06.999Z)

Cody.Cardenas:
> Ive not had access to bamboo since it was live, I’ll look into how to access it.
Every user has access to the three clusters:
https://doc.eresearch.unige.ch/hpc/access_the_hpc_clusters#cluster_connection[https://doc.eresearch.unige.ch/hpc/access_the_hpc_clusters#cluster_connection](https://doc.eresearch.unige.ch/hpc/access_the_hpc_clusters#cluster_connection)


## Post 5 by @Yann.Sagon (2024-12-19T10:45:00.589Z)

A post was split to a new topic: Use Beast checkpointing with slurm job array[Use Beast checkpointing with slurm job array](https://hpc-community.unige.ch/t/use-beast-checkpointing-with-slurm-job-array/3771)
