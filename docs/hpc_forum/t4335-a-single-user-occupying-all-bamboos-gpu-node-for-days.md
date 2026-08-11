# A single user occupying all Bamboo's GPU node for days

- Source: https://hpc-community.unige.ch/t/4335

- Created: 2026-07-02T07:23:32.364Z

- Posts: 5

- Category: 8

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Njiva.Andrianarivelo (2026-07-02T07:23:32.425Z)

Hello,
user hajianf3 has occupied all 10 public GPU nodes on bamboo for days now, bocking the access to other users. I’m just asking if it’s normal or authorized to occupy 10 GPU nodes for that long ?
In my lab we try to use max two node at the same time for heavy works to not bother other people’s work.
Best,
Andry


## Post 2 by @Adrien.Albert (2026-07-02T09:38:59.139Z)

Hello @Njiva.Andrianarivelo[@Njiva.Andrianarivelo](https://hpc-community.unige.ch/u/njiva.andrianarivelo),
Thank you for raising your concern.
First, I would like to kindly remind  that the forum is intended for technical discussions and user support, not for publicly singling out individual users. If you have identified a specific user and believe there may be an issue, please send us the relevant details by email, we will be happy to review the situation.
Regarding the situation you observed, the current behavior appears to be consistent with Slurm’s scheduling policies.
Slurm does not schedule jobs purely on a first-come, first-served basis. Job priority is computed from several factors, including:
- Fairshare usage
- Job age (time spent waiting in the queue)
- Partition and QoS settings
- Requested resources
- Other site-specific scheduling parameters
The following command may help to better understand how Slurm computes job priorities:
```
(bamboo)-[root@admin1 ~]$ sprio -l
          JOBID PARTITION     USER    ACCOUNT   PRIORITY       SITE        AGE      ASSOC  FAIRSHARE    JOBSIZE  PARTITION    QOSNAME        QOS        NICE                 TRES
        3991335 shared-gp    vador    Palpatine   430376          0      38599          0      16768          9     375000     normal          0           0                     
        4026427 shared-gp     maul    Palpatine   376508          0       1475          0          0         34     375000     normal          0           0                     
        4026478 shared-gp Assaj_Ven+  Palpatine   378424          0       3416          0          0          9     375000     normal          0           0
```
---
One important factor is job age: jobs that remain in the queue for a longer period gradually gain priority. After reviewing the queue, I observed that this user’s jobs had been waiting longer than several other pending jobs, which contributed to their higher priority when resources became available.
It is also worth noting that the user is not exclusively occupying all GPU resources, as GPU jobs from other users are running at the same time.
In addition, Slurm uses a backfill scheduler. This mechanism continuously re-evaluates the queue and starts jobs whenever it can do so without delaying higher-priority reservations. The objective is to maximize cluster utilization while preserving scheduling fairness. As a result, it is not unusual to see multiple jobs from the same user start simultaneously when resources become available.
Another important aspect is the requested walltime and resources. The backfill scheduler can generally place jobs more easily when they request realistic resource requirements and an accurate walltime estimate. Conversely, jobs requesting significantly more resources than needed, or requesting the maximum walltime allowed by the partition, are often more difficult to fit into the schedule and may remain queued longer. Therefore, accurately estimating resource requirements and walltime is beneficial both for cluster efficiency and for reducing scheduling delays.
Finally, many users submit workloads through sbatch, often using job arrays. These jobs can remain in the queue for extended periods and continue to accumulate priority while waiting. This is expected behavior and is one of the reasons why batch submission is generally recommended for large-scale workloads.
Best regards,


## Post 3 by @Njiva.Andrianarivelo (2026-07-02T10:46:17.587Z)

Thanks for the clarification and sorry for the rude approach.
Best,


## Post 4 by @Adrien.Albert (2026-07-02T13:19:26.543Z)

Hi @Njiva.Andrianarivelo[@Njiva.Andrianarivelo](https://hpc-community.unige.ch/u/njiva.andrianarivelo),
I also forgot to mention that you seems have access to a private partition on Baobab :slightly_smiling_face:.
You can run your jobs on this partition and benefit from its scheduling priority.
```
(baobab)-[XXXX@admin1 ~]$ sinfo |grep gpu
private-bellone-gpu       up 7-00:00:00      1    mix gpu046
```
You may also consider using other cluster (Baobab or Yggdrasil), where GPU resources may be available depending on the current cluster load.
Best


## Post 5 by @Njiva.Andrianarivelo (2026-07-02T15:44:03.025Z)

That’s what we are doing moste of the time :slight_smile:
Best,
