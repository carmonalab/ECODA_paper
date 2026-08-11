# Baobab keeps crashing for no reason

- Source: https://hpc-community.unige.ch/t/3842

- Created: 2025-02-26T15:08:27.347Z

- Posts: 15

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @nicolas.clairis (2025-02-26T15:08:27.391Z)

## Primary informations
Username: clairis
Cluster: Baobab
## Description
I tried a lot of different version of my batch, but it keeps crashing all the time. My goal is just to run fMRI preprocessing in Matlab for 74 subjects in parallel. I initially tried with a parfor loop in Matlab, then abandoned it and replaced it by a loop inside the sbatch using the --array command but whatever I do the batch crashes after around 1h without any error message but without finishing the job it was supposed to do. I don’t understand why it keeps failing, especially since I don’t get any error message in Matlab or in the batch suggesting that the problem lies in the way Slurm interacts with my script (RAM issue?)
## Steps to Reproduce
#!/bin/bash
#SBATCH --job-name=preprocessing
#SBATCH --partition=public-cpu
#SBATCH --mem=8G # define RAM to use in total
#SBATCH --array=1-73
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=3-00:00:00
#SBATCH --chdir=/home/users/c/clairis/scripts/
#SBATCH --output=logs/name=preprocessing_matlab_output_%j_%a.txt
#SBATCH --error=logs/name=preprocessing_matlab_output_%j_%a.txt
module load MATLAB
SUBJECTS=(“E002” “E004” “E006” “E007” “E010” “E011” “E012” “E014” “E015” “E017” “E019” “E020” “E021” “E022” “E023” “E028” “E035” “E036” “E037” “E043” “E044” “E048” “E050” “E052” “E055” “E056” “E059” “E061” “E063” “E064” “E068” “E070” “E076” “E078” “E084” “E090” “E092” “E093” “E100” “E101” “E107” “E108” “E109” “E111” “E112” “E114” “E115” “E117” “E123” “E124” “E125” “E126” “E128” “E130” “E131” “E133” “E135” “E136” “E137” “E138” “E139” “E140” “E141” “E142” “E145” “E146” “E148” “E149” “E150” “E151” “E153” “E158” “E161”);
SUBJECT=${SUBJECTS[$SLURM_ARRAY_TASK_ID-1]}
echo “Processing subject $SUBJECT on node $SLURMD_NODENAME”
matlab -nodesktop -softwareopengl -nosplash -nodisplay -r “addpath(‘/home/users/c/clairis/scripts’); preprocessing_loop(‘$SUBJECT’); exit;” # adds script folder path + launches the script
## Expected Result
The batch should have launched the 73 corresponding subjects in parallel and I should have my preprocessing done. Instead, it starts doing it but crashes randomly about 1h or so after the start after performing 4-6 subjects (depending on the number of available nodes) and then crashes with no further explanation right in the middle…
## Actual Result
The log files indicate that the script crashes randomly but before finishing.
I managed to make it finish one subject once when asking for just 2 nodes to work in parallel as if by reducing the number of jobs it would be more efficient. I’m currently considering to remove the parallelizing (which is the main interest of the cluster) or to switch back to a local analysis as this has already cost me weeks without understanding what is going wrong so any help would be more than welcome!


## Post 2 by @Adrien.Albert (2025-02-26T15:51:47.615Z)

Hi @nicolas.clairis[@nicolas.clairis](https://hpc-community.unige.ch/u/nicolas.clairis)
```
(baobab)-[root@login1 ~]$ scontrol show job 15203285_73
JobId=15203285 ArrayJobId=15203285 ArrayTaskId=73 JobName=wololo_preprocessing
   UserId=clairis(452609) GroupId=hpc_users(5000) MCS_label=N/A
   Priority=1003854 Nice=0 Account=cpiguet QOS=normal
   JobState=FAILED Reason=RaisedSignal:53(Real-time_signal_19) Dependency=(null)
```
is this topic means something to you ? :
stackoverflow.com[stackoverflow.com](https://stackoverflow.com/questions/77642713/all-slurm-jobs-fail-silently-with-exit-code-053)
Cornelius Roemer
#### All slurm jobs fail silently with exit code 0:53[All slurm jobs fail silently with exit code 0:53](https://stackoverflow.com/questions/77642713/all-slurm-jobs-fail-silently-with-exit-code-053)
slurm, sbatch
asked by
  
  
    Cornelius Roemer
  [Cornelius Roemer](https://stackoverflow.com/users/7483211/cornelius-roemer)
  on 10:30PM - 11 Dec 23 UTC[10:30PM - 11 Dec 23 UTC](https://stackoverflow.com/questions/77642713/all-slurm-jobs-fail-silently-with-exit-code-053)


## Post 3 by @nicolas.clairis (2025-02-26T16:07:17.854Z)

Hi
Indeed the message I receive in my email box is the following:
“Slurm Array Summary Job_id=15203285_* (15203285) Name=preprocessing Ended, Mixed, MaxSignal [53]”
so I also have the code 53. However, the path is correct in my case (there is a writtable scripts/logs folder where the output log files are written without any issue). Could the problem come from the fact that in my batch the output and error log files have the same name?
In the thread you shared, someone also says “the same thing happens if you exceed your quota” which would be consistent with the fact that my batch often crashes around 1h after starting but I don’t know why I would have such a limited quota though?
I will try again by changing the error logfile name, but I doubt it comes from it as I think I already had a batch like that in another institute and it was not a problem


## Post 4 by @nicolas.clairis (2025-02-26T16:35:14.451Z)

It’s ongoing with the job number 15203934


## Post 5 by @nicolas.clairis (2025-02-26T16:43:59.050Z)

Not sure if it will crash or not, but I already see that this led to the creation of the corresponding error files while they were not created before. Thanks for the tip. Let’s see if that solves the issue


## Post 6 by @nicolas.clairis (2025-02-26T17:12:39.298Z)

So it seems that there is one subject which crashes leading to the whole script to crash (not sure why Slurm doesn’t compute the others when one crashes) but now I have the error log file to see it as well


## Post 7 by @Adrien.Albert (2025-02-27T09:33:00.784Z)

Hi @nicolas.clairis[@nicolas.clairis](https://hpc-community.unige.ch/u/nicolas.clairis)
nicolas.clairis:
> Could the problem come from the fact that in my batch the output and error log files have the same name?
Multiple processes writing to the same file can cause data corruption and inconsistencies, so it’s not recommended. But I think SLURM just splits `stderr` and `stdout`, so it should not be an issue.
nicolas.clairis:
> I don’t know why I would have such a limited quota though
https://doc.eresearch.unige.ch/hpc/storage_on_hpc#check_disk_usage_on_the_clusters[https://doc.eresearch.unige.ch/hpc/storage_on_hpc#check_disk_usage_on_the_clusters](https://doc.eresearch.unige.ch/hpc/storage_on_hpc#check_disk_usage_on_the_clusters)
```
(baobab)-[root@login1~]$ beegfs-get-quota-home-scratch.sh -u !$
beegfs-get-quota-home-scratch.sh -u clairis
home dir: /home/users/c/clairis
scratch dir: /srv/beegfs/scratch/users/c/clairis

          user/group                 ||           size          ||    chunk files
  storage     |   name        |  id  ||    used    |    hard    ||  used   |  hard
  ----------------------------|------||------------|------------||---------|---------
home        |        clairis|452609||  673.56 GiB| 1024.00 GiB||  1741262|unlimited
scratch     |        clairis|452609||      0 Byte|   unlimited||        0| 10000000
```
Your quota seems OK.
nicolas.clairis:
> not sure why Slurm doesn’t compute the others when one crashes
- Is there result dependency beetween previous jobs ?
- Could you give the path of the log directory to check the errors ?


## Post 8 by @nicolas.clairis (2025-02-27T09:47:40.895Z)

Adrien.Albert:
> Could you give the path of the log directory to check the errors ?
I just ran another batch (15211163_) which also crashed despite fixing the two previously mentionned issues (one subject had an issue with a path definition + the error file output now has a different name). Again everything crashed ~30min after starting. All the error logs are empty (can be found in /home/users/c/clairis/scripts/logs) and the batch script shows that all subjects were interrupted in the beginning of the analysis with no subject actually completed.
Adrien.Albert:
> Is there result dependency beetween previous jobs ?
No, all jobs (i.e. all subjects) should be completely independent.
I didn’t re-run a test recently, but the only subject that could be run until the end was due to the fact that I reduced the number of subjects processed in parallel to 2. Could there be an issue related to the fact that I give the full array of 73 individuals? Should I specify a smaller number via `--array=1-73%5` for example to reduce the number of subjects being processed in parallel? I don’t know why that would be an issue though as the cluster should be precisely optimized to deal with that kind of situations right?


## Post 9 by @nicolas.clairis (2025-02-27T09:48:35.652Z)

the email message is the following:
`Slurm Array Summary Job_id=15211163_* (15211163) Name=wololo_preprocessing Ended, Mixed, MaxSignal [53]`
so again code 53 strangely


## Post 10 by @nicolas.clairis (2025-02-27T11:00:05.463Z)

This is my current version of the batch:
```
#!/bin/bash
#SBATCH --job-name=wololo_preprocessing
#SBATCH --mail-user=xxx@unige.ch
#SBATCH --mail-type=END
#SBATCH --partition=public-cpu
#SBATCH --mem-per-cpu=32G
#SBATCH --array=1-73
#SBATCH --ntasks-per-node=1 # number of tasks (i.e. subjects) per node
#SBATCH --time=3-00:00:00
#SBATCH --chdir=/home/users/c/clairis/scripts/
#SBATCH --output=logs/name=preprocessing_matlab_output_%j_%a.txt 
#SBATCH --error=logs/name=preprocessing_matlab_error_%j_%a.txt 
module load MATLAB
SUBJECTS=("E002" "E004" "E006" "E007" "E010" "E011" "E012" "E014" "E015" "E017" "E019" "E020" "E021" "E022" "E023" "E028" "E035" "E036" "E037" "E043" "E044" "E048" "E050" "E052" "E055" "E056" "E059" "E061" "E063" "E064" "E068" "E070" "E076" "E078" "E084" "E090" "E092" "E093" "E100" "E101" "E107" "E108" "E109" "E111" "E112" "E114" "E115" "E117" "E123" "E124" "E125" "E126" "E128" "E130" "E131" "E133" "E135" "E136" "E137" "E138" "E139" "E140" "E141" "E142" "E145" "E146" "E148" "E149" "E150" "E151" "E153" "E158" "E161");
SUBJECT=${SUBJECTS[$SLURM_ARRAY_TASK_ID-1]}
echo "Processing subject $SUBJECT on node $SLURMD_NODENAME"
matlab -nodesktop -softwareopengl -nosplash -nodisplay -r "addpath('/home/users/c/clairis/scripts'); preprocessing_loop('$SUBJECT'); exit;" # adds script folder path + launches the script
```
Any idea where the error comes from ?


## Post 11 by @nicolas.clairis (2025-02-27T14:27:32.356Z)

I’m going to try to reduce the number of jobs performed in parallel with the following command:
`#SBATCH --array=1-73%5`
and I removed any constraint on the RAM so that Slurm allocates it adaptively. Let’s see if these changes improve things, but would be nicer to be maximize the number of participants performed in parallel though.


## Post 12 by @nicolas.clairis (2025-02-27T15:37:40.214Z)

Sadly, the script crashed again despite these changes and still no error message in the logs or in the matlab outputs.


## Post 13 by @Adrien.Albert (2025-02-28T09:56:44.297Z)

Waiting for our meeting with Nicolas to get more details.


## Post 14 by @nicolas.clairis (2025-03-05T08:47:57.348Z)

After discussing with Adrien, I moved all my files from the /home to the /scratch folder and redid all the analysis after changing all the paths and everything seems to have worked this time. It seems that the issue was coming from the fact that the analysis was saturating the allowed space in the home folder (1To) and my script was crashing whenever that limit was reached (hence the apparent random timing of when it was crashing).


## Post 15 by @Adrien.Albert (2025-03-05T09:26:20.503Z)

@nicolas.clairis[@nicolas.clairis](https://hpc-community.unige.ch/u/nicolas.clairis)
I would like to add a small clarification: from my point of view, I suspect that home storage is almost full rather than your quota limit, (both are possible):
```
(baobab)-[root@admin1 ~]$ df /home -h
Filesystem      Size  Used Avail Use% Mounted on
beegfs_home     138T  137T  1.9T  94% /home
```
We’ve sent an email to the main users (>500GB) inviting them to empty their homes as far as possible.
