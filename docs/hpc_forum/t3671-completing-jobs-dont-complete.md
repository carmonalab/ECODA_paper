# Completing jobs don't complete

- Source: https://hpc-community.unige.ch/t/3671

- Created: 2024-10-03T09:26:53.589Z

- Tags: yggdrasil

- Posts: 3

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Matthias.Kruckow (2024-10-03T09:26:53.638Z)

## Primary informations
Username: kruckow (and others)
Cluster: Yggdrasil
## Description
I was checking the jobs which were run recently and I was missing some, which should have finished by now but aren’t. I finally found them in completing state CG. Interestingly, they have still 0 time and the SLURM output files are not created yet either. Here a list, when looking for completing jobs:
```
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
       35596853_98 private-a alpha1e-  nigioni CG       0:00      1 cpu145
       35596853_97 private-a alpha1e-  nigioni CG       0:00      1 cpu145
       35596853_96 private-a alpha1e-  nigioni CG       0:00      1 cpu145
       35596853_95 private-a alpha1e-  nigioni CG       0:00      1 cpu145
       35596853_94 private-a alpha1e-  nigioni CG       0:00      1 cpu145
       35596853_93 private-a alpha1e-  nigioni CG       0:00      1 cpu145
       35596853_92 private-a alpha1e-  nigioni CG       0:00      1 cpu145
       35596853_91 private-a alpha1e-  nigioni CG       0:00      1 cpu145
       35596853_86 private-a alpha1e-  nigioni CG       0:00      1 cpu145
       35596853_84 private-a alpha1e-  nigioni CG       0:00      1 cpu145
       35596853_83 private-a alpha1e-  nigioni CG       0:00      1 cpu145
       35596853_82 private-a alpha1e-  nigioni CG       0:00      1 cpu146
       35596853_79 private-a alpha1e-  nigioni CG       0:00      1 cpu146
       35596853_65 private-a alpha1e-  nigioni CG       0:00      1 cpu145
       35596853_64 private-a alpha1e-  nigioni CG       0:00      1 cpu145
       35596853_63 private-a alpha1e-  nigioni CG       0:00      1 cpu145
       35596853_54 private-a alpha1e-  nigioni CG       0:00      1 cpu146
       35596853_53 private-a alpha1e-  nigioni CG       0:00      1 cpu146
       35596853_51 private-a alpha1e-  nigioni CG       0:00      1 cpu146
       35596853_50 private-a alpha1e-  nigioni CG       0:00      1 cpu146
       35596853_28 private-a alpha1e-  nigioni CG       0:00      1 cpu145
       35596853_22 private-a alpha1e-  nigioni CG       0:00      1 cpu145
       35596853_21 private-a alpha1e-  nigioni CG       0:00      1 cpu145
       35596853_20 private-a alpha1e-  nigioni CG       0:00      1 cpu145
       35596853_18 private-a alpha1e-  nigioni CG       0:00      1 cpu145
       35596853_17 private-a alpha1e-  nigioni CG       0:00      1 cpu145
       35596853_16 private-a alpha1e-  nigioni CG       0:00      1 cpu145
       35596853_15 private-a alpha1e-  nigioni CG       0:00      1 cpu145
       35596853_14 private-a alpha1e-  nigioni CG       0:00      1 cpu145
       35596853_13 private-a alpha1e-  nigioni CG       0:00      1 cpu145
       35596853_12 private-a alpha1e-  nigioni CG       0:00      1 cpu145
       35596811_81 private-a alpha1e-  nigioni CG       0:00      1 cpu145
       35596811_80 private-a alpha1e-  nigioni CG       0:00      1 cpu145
       35596811_79 private-a alpha1e-  nigioni CG       0:00      1 cpu145
       35596811_77 private-a alpha1e-  nigioni CG       0:00      1 cpu145
       35596811_76 private-a alpha1e-  nigioni CG       0:00      1 cpu145
       35596811_74 private-a alpha1e-  nigioni CG       0:00      1 cpu145
       35596811_60 private-a alpha1e-  nigioni CG       0:00      1 cpu146
       35596811_57 private-a alpha1e-  nigioni CG       0:00      1 cpu145
       35596811_56 private-a alpha1e-  nigioni CG       0:00      1 cpu145
       35596811_55 private-a alpha1e-  nigioni CG       0:00      1 cpu145
       35596811_54 private-a alpha1e-  nigioni CG       0:00      1 cpu145
       35596811_53 private-a alpha1e-  nigioni CG       0:00      1 cpu145
       35596811_20 private-a alpha1e-  nigioni CG       0:00      1 cpu145
       35596811_18 private-a alpha1e-  nigioni CG       0:00      1 cpu145
       35596811_17 private-a alpha1e-  nigioni CG       0:00      1 cpu145
       35596811_16 private-a alpha1e-  nigioni CG       0:00      1 cpu145
        35596811_4 private-a alpha1e-  nigioni CG       0:00      1 cpu145
        35596811_3 private-a alpha1e-  nigioni CG       0:00      1 cpu145
        35596811_2 private-a alpha1e-  nigioni CG       0:00      1 cpu145
        35596811_1 private-a alpha1e-  nigioni CG       0:00      1 cpu145
       35659774_24 private-a perLTarr   pedros CG       0:00      1 cpu146
       35659774_22 private-a perLTarr   pedros CG       0:00      1 cpu146
       35659774_21 private-a perLTarr   pedros CG       0:00      1 cpu146
       35659774_12 private-a perLTarr   pedros CG       0:00      1 cpu146
       35659774_11 private-a perLTarr   pedros CG       0:00      1 cpu146
       35659774_10 private-a perLTarr   pedros CG       0:00      1 cpu146
        35537989_8 private-a mesa_gri  kruckow CG       0:00      1 cpu146
        35537989_5 private-a mesa_gri  kruckow CG       0:00      1 cpu146
        35537989_4 private-a mesa_gri  kruckow CG       0:00      1 cpu146
        35537989_3 private-a mesa_gri  kruckow CG       0:00      1 cpu146
        35537989_0 private-a mesa_gri  kruckow CG       0:00      1 cpu146
```
They look to be all concentrated on two nodes: cpu145 and cpu146. According to the output of `scontrol` SLURM has last evaluated them on Sep. 30 `LastSchedEval=2024-09-30`. Thus, I think there is a manual action needed to let SLURM finish those jobs. The timing looks to coincide with restoring the node which were shut down, see Full partitions are down on yggdrasil[Full partitions are down on yggdrasil](https://hpc-community.unige.ch/t/full-partitions-are-down-on-yggdrasil/3653).
I’m not sure about others, but at least my jobs have dependencies on them, hence I’d prefer a solution, which let them run as they should have already. If there is no option to restore them, I guess they need to get canceled (in case, please inform all the effected users, which of their jobs got canceled).


## Post 2 by @Matthias.Kruckow (2024-10-07T11:07:20.434Z)

An update: it looks like that SLURM was able to repair the issue itself:
An hour ago the jobs got requeued.
I think, the key to trigger the requeuing was the fact, that today we reached the 7 days time limit of the partition private-astro-cpu (those jobs submitted to).
P.S.: There are again 2 nodes down on private-astro-cpu. I hope the rest won’t follow soon again.


## Post 3 by @Gael.Rossignol (2024-10-09T07:58:09.080Z)

Dear Matthias,
Thanks for your message and excellent investigations. There are no interresting informations regarding those jobs in the log for example :
```
[2024-10-08T13:53:41.784] sched: Allocate JobId=35537989_0(35659888) NodeList=cpu146 #CPUs=4 Partition=private-astro-cpu
[2024-10-09T01:54:32.260] _job_complete: JobId=35537989_0(35659888) WEXITSTATUS 0
[2024-10-09T01:54:32.262] _job_complete: JobId=35537989_0(35659888) done
```
If you encourter again the issue you can post it here.
Best regards,
