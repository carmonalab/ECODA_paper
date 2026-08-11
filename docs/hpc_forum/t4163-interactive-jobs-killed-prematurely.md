# Interactive jobs killed prematurely

- Source: https://hpc-community.unige.ch/t/4163

- Created: 2025-12-08T14:55:26.105Z

- Tags: baobab

- Posts: 20

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Vilius.Cepaitis (2025-12-08T14:55:26.153Z)

Dear HPC experts,
Ever since the baobab maintenance finished last week, my interactive jobs started to get killed before they are meant to expire.
For example, today I submitted jobs with 8 hours runtime:
`salloc -n1 -c16 --partition=public-cpu,private-dpnc-cpu,shared-cpu --time=08:00:00 --mem=20G`
But each time the job got killed after about 1-2 hours:
```
➜ weakly-supervised-search (cathode-improvements) ✗ srun: Job step aborted: Waiting up to 92 seconds for job step to finish.
salloc: Job allocation 6129035 has been revoked.
                                                [2025-12-08T15:44:50.841] error: *** STEP 6129035.interactive ON cpu089 CANCELLED AT 2025-12-08T15:44:50 DUE TO TIME LIMIT ***
srun: error: cpu089: task 0: Killed
```
This happened at least twice on different nodes. Could you please help me with this?
Best regards,
Vilius


## Post 2 by @Adrien.Albert (2025-12-09T09:10:54.025Z)

Dear @Vilius.Cepaitis[@Vilius.Cepaitis](https://hpc-community.unige.ch/u/vilius.cepaitis),
Thank you very much for the details provided. :folded_hands:
The issue you are experiencing appears to be similar to the one discussed here: Job timeout despite not hitting the timelimit - #8 by Michael.Sonner[Job timeout despite not hitting the timelimit - #8 by Michael.Sonner](https://hpc-community.unige.ch/t/job-timeout-despite-not-hitting-the-timelimit/3042/8).
Based on the logs:
```
[2025-12-08T14:56:12.460] sched: _slurm_rpc_allocate_resources JobId=6129035 NodeList=(null) usec=5978
[2025-12-08T14:56:13.005] sched: Allocate JobId=6129035 NodeList=cpu089 #CPUs=16 Partition=private-dpnc-cpu
[2025-12-08T15:41:21.003] job_time_limit: inactivity time limit reached for JobId=6129035
[2025-12-08T15:46:20.880] cleanup_completing: JobId=6129035 completion process took 299 seconds
```
To clarify, was the job terminated while a process was still running under `salloc`?
We are also considering that this may be related to network issues, which could cause similar behaviour, and we are currently investigating this as part of the ongoing cluster incident:
[2025] Current issues on HPC Cluster[[2025] Current issues on HPC Cluster](https://hpc-community.unige.ch/t/2025-current-issues-on-hpc-cluster/3788/27) HPC ChangeLog[HPC ChangeLog](https://hpc-community.unige.ch/c/hpc-announce/hpc-changelog/9)
> Description
Cluster: Yggdrasil 
We are currently experiencing issues with Slurm on Yggdrasil, which are causing unexpected behavior (delays when starting tasks, failures, etc.). We are working to resolve this issue as quickly as possible. 
– 
HPC Team 
Status : In progress orange_circle
start: 2025-12-05T14:30:00Z (UTC) 
end:


## Post 3 by @Vilius.Cepaitis (2025-12-09T10:07:14.930Z)

Hi @Adrien.Albert[@Adrien.Albert](https://hpc-community.unige.ch/u/adrien.albert),
Thank you very much for your response and for looking into this.
I was running a vscode server and snakemake on the node in question, but perhaps the scheduler thought the job was idle. Also, I ran salloc within tmux when it got killed.
Best regards,
Vilius


## Post 4 by @Vilius.Cepaitis (2026-01-14T20:03:11.541Z)

Hello,
I’m still experiencing this issue.
I need my jobs to run non-interrupted as I have snakemake running on the allocated job managing my workflow, which can take up to a day to finish. Do you have any suggestions what to do?
Best regards,
Vilius
> salloc: Job allocation 6547165 has been revoked.
> srun: Job step aborted: Waiting up to 92 seconds for job step to finish.
> [2026-01-14T17:37:06.137] error: *** STEP 6547165.interactive ON cpu239 CANCELLED AT 2026-01-14T17:37:06 DUE TO TIME LIMIT ***
> srun: error: cpu239: task 0: Killed


## Post 5 by @maciej.falkiewicz (2026-01-15T08:36:43.692Z)

+1 with the same type of issue


## Post 6 by @Adrien.Albert (2026-01-15T08:37:54.434Z)

Hello;
Could you please give me your salloc/srun cmd ?


## Post 7 by @Vilius.Cepaitis (2026-01-15T08:39:02.738Z)

Hi @Adrien.Albert[@Adrien.Albert](https://hpc-community.unige.ch/u/adrien.albert),
It’s typically something like this
> salloc -n1 -c16 --partition=public-cpu,private-dpnc-cpu,shared-cpu --time=24:00:00 --mem=30G
Cheers,
Vilius


## Post 8 by @maciej.falkiewicz (2026-01-15T11:31:32.234Z)

```
srun -p shared-gpu,private-kalousis-gpu,private-cui-gpu --gres=gpu:4,VramPerGpu:48G --time=7-00:00:00 --mem=128000 --cpus-per-task=32 [...]
```


## Post 9 by @Adrien.Albert (2026-01-15T13:45:24.238Z)

Hello Vilius;
```
[2026-01-14T17:30:45.011] job_time_limit: inactivity time limit reached for JobId=6547165
[2026-01-14T17:38:35.992] cleanup_completing: JobId=6547165 completion process took 470
 seconds
```
Your job has been killed for ‘inactivity’.
The SLURM FAQ:
> Why is my job killed prematurely?
> Slurm has a job purging mechanism to remove inactive jobs (resource allocations) before reaching its time limit, which could be infinite. This inactivity time limit is configurable by the system administrator. You can check its value with the command
> > scontrol show config | grep InactiveLimit
> The value of InactiveLimit is in seconds. A zero value indicates that job purging is disabled. A job is considered inactive if it has no active job steps or if the srun command creating the job is not responding. In the case of a batch job, the srun command terminates after the job script is submitted. Therefore batch job pre- and post-processing is limited to the InactiveLimit. Contact your system administrator if you believe the InactiveLimit value should be changed.
On Unige Cluster:
```
(baobab)-[root@login1 ~]$ scontrol show config | grep InactiveLimit
InactiveLimit           = 300 sec
```
We have test with screen and we have no issue. However you are using tmux, we suspect it pause/disconnect the ssh process.
@maciej.falkiewicz[@maciej.falkiewicz](https://hpc-community.unige.ch/u/maciej.falkiewicz) are you using Tmux or another tool instead shell ?
It could be a network issue beetween your machine and the cluster impacting the ssh process.


## Post 10 by @Adrien.Albert (2026-01-15T13:46:33.886Z)

we also test with salloc during one hour doing nothing and the job continue without being killed for inactivity.


## Post 11 by @maciej.falkiewicz (2026-01-15T18:34:01.144Z)

> @maciej.falkiewicz[@maciej.falkiewicz](https://hpc-community.unige.ch/u/maciej.falkiewicz) are you using Tmux or another tool instead shell ?
> It could be a network issue beetween your machine and the cluster impacting the ssh process.
yes @Adrien.Albert[@Adrien.Albert](https://hpc-community.unige.ch/u/adrien.albert) , I do use tmux. There can be a network issue between login node and compute node, but the termination wouldn’t happen with TIMEOUT, no?
I can try the same with sbatch if you reserve the node for me :slight_smile:


## Post 12 by @Yann.Sagon (2026-01-16T10:44:05.895Z)

@Vilius.Cepaitis[@Vilius.Cepaitis](https://hpc-community.unige.ch/u/vilius.cepaitis)
my two cents to optimize your job submission: the shared-cpu max time is 12h and you are requesting 24h. Do not add this partition in the list of requested partition OR use 12h00 as time limit.
This won’t fix your issue but it is better for the scheduling.


## Post 13 by @Adrien.Albert (2026-01-16T10:56:41.290Z)

Yes, the job status shows TIMEOUT, but the Slurm logs indicate also the cause is inactivity, meaning Slurm considered the interactive session (`salloc` / `srun`) as not responding.
We also know this does not happen with `sbatch`, because batch jobs do not depend on the stability of your SSH connection.
Another user recently had exactly the same issue, and migrating to `sbatch` fully solved his problem.
For this reason, we very strongly recommend running jobs via `sbatch` whenever possible.
With `salloc` and `srun`, any network disruption can cause Slurm to kill the job.
Have you observed the same behavior when using OpenOnDemand?
OOD keeps the session on the server side, so it is more tolerant to network instability.
Also keep in mind: a laptop that goes into sleep mode or closes its lid often loses network connectivity, which can trigger this situation.
## My test:
To confirm the root cause, I performed the following test as user:
- SSH connection to login node
- Start a `screen` session:`screen -S test`
- Run `salloc` inside this `screen`
- Detach from the screen (`Ctrl-A D`)
- Disconnect from the login node entirely
The `salloc` is doing “nothing” (just an idle shell waiting for input).
Result:
The `salloc` allocation continued running for more than 30 minutes without interruption, even though the SSH connection was closed.
Example (`sacct` output):
```
(baobab)-[root@login1 ~]$ sac -u alberta
          JobID    JobName    Account      User        NodeList  ReqCPUS     ReqMem   NTasks               Start                 End    Elapsed      State 
--------------- ---------- ---------- --------- --------------- -------- ---------- -------- ------------------- ------------------- ---------- ---------- 
        6562597 interacti+      burgi   alberta          cpu349        1      3000M          2026-01-16T11:02:56             Unknown   00:34:13    RUNNING 
```
This test orients my suspicion on network interruptions between your laptop and the login node.


## Post 14 by @maciej.falkiewicz (2026-01-16T13:09:56.741Z)

Yann.Sagon:
> my two cents to optimize your job submission: the shared-cpu max time is 12h and you are requesting 24h. Do not add this partition in the list of requested partition OR use 12h00 as time limit.
actually, this sounds like a possible explanation. My job was scheduled for 7days, but I had shared-gpu in the partition’s list. But the TIMEOUT was after ~24 hours, not 12 :slight_smile:
After the TIMEOUT took place, I scheduled the same job but only with private partitions, and the job is running 1day +
Update: the new job died after 1-18:06:59 :joy:
Adrien.Albert:
> Another user recently had exactly the same issue, and migrating to `sbatch` fully solved his problem.
@Adrien.Albert[@Adrien.Albert](https://hpc-community.unige.ch/u/adrien.albert) thanks for the info, switching to sbatch :saluting_face:
Thanks for your time on investigations!


## Post 15 by @Adrien.Albert (2026-01-16T13:45:14.016Z)

@maciej.falkiewicz[@maciej.falkiewicz](https://hpc-community.unige.ch/u/maciej.falkiewicz)
Just to be sure how are you using tmux?
- Do you connect to the login node first and then start tmux there?
- Or do you start tmux on your laptop and then connect to the login node afterward?
My entire analysis was based on scenario #2.


## Post 16 by @maciej.falkiewicz (2026-01-16T13:58:47.403Z)

I understand that you have to ask such questions, but it is a bit offensive :slight_smile:
I was running tmux on the login node.


## Post 17 by @Vilius.Cepaitis (2026-01-16T14:35:23.985Z)

Hi @Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon) and @Adrien.Albert[@Adrien.Albert](https://hpc-community.unige.ch/u/adrien.albert)
Thank you very much for your investigations and suggestions. I will also switch to using `sbatch`.
Perhaps as a useful datapoint, my setup was also connecting to the login node first, starting a tmux session there and running `salloc` from inside tmux when the timeouts happened.
Cheers,
Vilius


## Post 18 by @Malte.Algren (2026-01-16T17:02:17.125Z)

Hi,
Same as Vilius here.
I’m using tmux and salloc, and ever since december i have also experienced this premature timeout.
I understand that yes, we could also use sbatch for this, but I think a lot of people have been using tmux + salloc and this has suddenly become a problem (without us modifying our workflow)
This happened to me today (and yes it is a tmux session)
```
(baobab)-\[algren@gpu023 resources\]$ ctsalloc: Job allocation 6561472 has been revoked.
srun: Job step aborted: Waiting up to 92 seconds for job step to finish.
\[2026-01-16T14:02:17.857\] error: \*\*\* STEP 6561472.interactive ON gpu023 CANCELLED AT 2026-01-16T14:02:17 DUE TO TIME LIMIT \*\*\*
srun: error: gpu023: task 0: Killed
(baobab)-\[algren@login1 \~\]$
(baobab)-\[algren@login1 \~\]$
(baobab)-\[algren@login1 \~\]$
…
(baobab)-\[algren@login1 \~\]$  salloc -c4 --partition=private-dpnc-gpu,shared-gpu, --time=00-12:00:00 --mem=32GB --gres=gpu:1,VramPerGpu:2G
salloc: Pending job allocation 6570106
salloc: job 6570106 queued and waiting for resources
salloc: job 6570106 has been allocated resources
salloc: Granted job allocation 6570106
salloc: Waiting for resource configuration
salloc: Nodes gpu044 are ready for job
(baobab)-\[algren@gpu044 \~\]$ srun: Job step aborted: Waiting up to 92 seconds for job step to finish.
salloc: Job allocation 6570106 has been revoked.
\[2026-01-16T14:19:31.959\] error: \*\*\* STEP 6570106.interactive ON gpu044 CANCELLED AT 2026-01-16T14:19:31 DUE TO TIME LIMIT \*\*\*
srun: error: gpu044: task 0: Killed
(baobab)-\[algren@login1 \~\]$  salloc -c4 --partition=private-dpnc-gpu,shared-gpu, --time=01-12:00:00 --mem=32GB --gres=gpu:1,VramPerGpu:2G

(baobab)-\[algren@gpu023 \~\]$ diffgae
```


## Post 19 by @Adrien.Albert (2026-01-19T12:05:52.820Z)

Dear All,
We have openned an issue on SchedMD (slurm provider). Waiting for their answer !


## Post 20 by @Yann.Sagon (2026-01-19T13:29:51.686Z)

[2026] Current issues on HPC Cluster[[2026] Current issues on HPC Cluster](https://hpc-community.unige.ch/t/2026-current-issues-on-hpc-cluster/4185/8) HPC ChangeLog[HPC ChangeLog](https://hpc-community.unige.ch/c/hpc-announce/hpc-changelog/9)
> Description
Cluster: all 
Since the latest slurm update, some interactive jobs (using srun or salloc) are killed prematurely. For the user it appears as if the job reached its timelimit, but in the admin logs it is indicated that the job was killed due to inactivity timeout. We opened a case at schedMD. 
Status : Ongoing orange_circle
start: 2025-12-19T10:00:00Z (UTC) 
end:
