# Module fails to load on some arrays on Baobab

- Source: https://hpc-community.unige.ch/t/3335

- Created: 2024-02-24T08:45:09.746Z

- Posts: 6

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Simon.Hug (2024-02-24T08:45:09.803Z)

Dear all
I have, since yesterday (the same code produced no errors before), an
issue on the yggdrasil cluster. I submit a slurm-job with 20
arrays on public-cpu. The slurm-job loads the following modules:
module load GCC/11.3.0  OpenMPI/4.1.4 R/4.2.1
module load rgdal/1.6-6
While some jobs run perfectly fine, for some others I get the
following error message:
Lmod has detected the following error: The following module(s) are unknown:
“GCC/11.3.0”
It seems that on some cpus the GCC module can not be loaded. Any idea
how to address this issues?
thanks for your hlpe and best wishes, simon


## Post 2 by @Adrien.Albert (2024-02-26T12:41:08.790Z)

Dear @Simon.Hug[@Simon.Hug](https://hpc-community.unige.ch/u/simon.hug)
By creating a problem in a category on HPC issues - HPC Community[HPC issues - HPC Community](https://hpc-community.unige.ch/c/hpc-support/hpc-issues/) you automatically have a Template to fill, to ensure a good understanding for all and a quick answer, try to fill it as much as possible :pray:t3:
Could you give me more information about the node/jobID where you have this problem?


## Post 3 by @Simon.Hug (2024-02-26T12:41:19.285Z)

ooppss, this actually happened on baobab and not yggdrasil. best, simon


## Post 4 by @Simon.Hug (2024-02-27T07:25:50.796Z)

The job-id was 7644778: of the twenty (1-20) arrays that were started more than half ended with this error message, while a few went through. So for instance 7644778_18 ended with the error message, while 7644778_13 worked perfectly fine.
best wishes


## Post 5 by @Adrien.Albert (2024-02-27T10:21:00.365Z)

Dear Simon,
Thank you for this information. Indeed, 2 nodes did not have the module available, this has been fixed.
We apologize for the inconvenience.


## Post 6 by @Simon.Hug (2024-02-27T10:38:19.880Z)

thanks for looking into this and fixing this problem. best wishes, simon
