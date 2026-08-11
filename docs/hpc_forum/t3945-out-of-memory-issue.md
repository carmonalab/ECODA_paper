# OUT OF MEMORY issue

- Source: https://hpc-community.unige.ch/t/3945

- Created: 2025-05-02T08:38:42.451Z

- Tags: baobab

- Posts: 12

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @nicolas.clairis (2025-05-02T08:38:42.498Z)

## Primary informations
Username: $clairis
Cluster: $Baobab
## Description
With my colleague Camille Serquet we need to run some pretty heavy scripts in Matlab but we keep being stuck because we are getting an OUT OF MEMORY issue. Here is the last email I got for example:
`Slurm Job_id=17557933 Name=CAPS_1b Ended, Run time 1-05:26:08, OUT_OF_MEMORY`
And the matlab error log file is the following:
```
/var/spool/slurmd/job17557933/slurm_script: line 27: 1428884 Killed                  matlab -nodesktop -softwareopengl -nosplash -nodisplay -r "addpath('/home/users/c/clairis/scripts'); eSmile_CAPS_part1b; exit;"
slurmstepd: error: Detected 1 oom_kill event in StepId=17557933.batch. Some of the step tasks have been OOM Killed.
```
I set up the batch to use 250Go of RAM on the public-bigmem partition. Is there any way to use more amounts of RAM than that? When I tried at some point I was always getting an error message saying that it was not possible to go higher. I see that in principle shared-bigmem can go up to 500Go but since the script appears pretty long it may crash before the 12h limit… Is there any other cluster/partition where I could try with more RAM? Do you know what makes the script crash? Is it because it is automatically killed when the RAM being required is too high or is it because some other people also require RAM at the same time? It’s been a while me and my colleague Camille Serquet are stuck with this and are running our scripts over and over without being able to get to the end so we would like to find a solution to stop making the cluster run for nothing… Any help/advice would be welcome!
## Steps to Reproduce
My batch was:
```
#!/bin/bash
#SBATCH --job-name=CAPS_1b

# Send an email when the job is completed
#SBATCH --mail-user=nicolas.clairis@unige.ch
#SBATCH --mail-type=END

# Select partition
#SBATCH --partition=public-bigmem
#SBATCH --mem=250G

# Request CPU resource for a serial job
#SBATCH --ntasks=1

# Set your minimum acceptable walltime, format: day-hours:minutes:seconds
#SBATCH --time=4-00:00:00

# Set the working directory - put your MATLAB script there
#SBATCH --chdir=/home/users/c/clairis/scripts/
#SBATCH --output=logs/CAPS1/eSmile_CAPS1b_matlab_output_%j_%a.txt 
#SBATCH --error=logs/CAPS1/eSmile_CAPS1b_matlab_err_%j_%a.txt 

# load relevant softwares on Baobab (Matlab + SPM)
module load MATLAB

# lancer le script
matlab -nodesktop -softwareopengl -nosplash -nodisplay -r "addpath('/home/users/c/clairis/scripts'); eSmile_CAPS_part1b; exit;"
```


## Post 2 by @nicolas.clairis (2025-05-02T08:47:24.988Z)

This is the output when I write `sacct --format=Start,AveCPU,State,MaxRSS,JobID,NodeList,ReqMem --units=G -j 17557933`  to check the details of the memory usage:
```
Start     AveCPU      State     MaxRSS JobID               NodeList     ReqMem
------------------- ---------- ---------- ---------- ------------ --------------- ----------
2025-04-30T15:42:24            OUT_OF_ME+            17557933              cpu245       250G
2025-04-30T15:42:24 1-05:15:47 OUT_OF_ME+    234.70G 17557933.ba+          cpu245
2025-04-30T15:42:24   00:00:00  COMPLETED      0.00G 17557933.ex+          cpu245
```
so it seems that my batch did not even use the full 250Go of RAM that I required but then I struggle to understand why it crashed?


## Post 3 by @Gael.Rossignol (2025-05-02T09:21:49.265Z)

Dear Nicolas,
Since `sacct` samples memory usage, it only shows the last recorded value, which is why you don’t see the full 250G.
You could consider using the bigmem partition on Bamboo, as those servers have more memory available.
Best regards,


## Post 4 by @nicolas.clairis (2025-05-02T09:46:54.959Z)

Thanks for the advice. Is there a way to transfer the data easily from Baobab to Bamboo? I would rather avoiding copy-pasting everything back to my local pc and then back to Bamboo if there is a way to avoid doing that


## Post 5 by @Gael.Rossignol (2025-05-02T12:20:38.332Z)

Dear Nicolas,
You can directly copy data from each login to avoid moving on your computer. Some documentation are available on the link :
doc.eresearch.unige.ch[doc.eresearch.unige.ch](https://doc.eresearch.unige.ch/nasac/users_documentation?s%5B%5D=rsync)
### nasac:users_documentation [eResearch Doc][nasac:users_documentation [eResearch Doc]](https://doc.eresearch.unige.ch/nasac/users_documentation?s%5B%5D=rsync)
Best regards,


## Post 6 by @nicolas.clairis (2025-05-02T12:47:02.426Z)

Which paragraph are you referring to exactly? Isn’t there a way to do it directly via filezilla?


## Post 7 by @nicolas.clairis (2025-05-02T15:21:49.143Z)

Do you mean the " Data transfer ### Migrate data to UNIGE tape solution" section? Am I supposed to zip all the folder?


## Post 8 by @nicolas.clairis (2025-05-06T13:37:01.065Z)

Could you provide some more guidance please? The links referred to in the “data transfer” paragraph lead to a “404 page not found” page[page](https://gitlab.unige.ch/storage/nas/ug-nas/-/blob/main/misc/ug-nas-create-zip-archive-from-folder.sh)…
The data is very heavy (several To) so I would like to avoid losing several more weeks by transferring it locally, especially because I am afraid I won’t have enough space to do so in my local pc so if you can share your tips to make the transfer directly from Baobab to Bamboo I would be very glad to follow your advice!


## Post 9 by @Nicolas.Clairis1 (2025-05-07T09:13:07.637Z)

I managed to have access locally and to get the ug-nas-create-zip-archive-from-folder.sh function but still don’t get exactly how this works to give the Baobab vs Bamboo so I would be glad if you can walk me through a bit please.


## Post 10 by @Nicolas.Clairis1 (2025-05-07T15:30:00.978Z)

Isn’t it hpc:best_practices [eResearch Doc][hpc:best_practices [eResearch Doc]](https://doc.eresearch.unige.ch/hpc/best_practices?s%5B%5D=rsync#rsync) what I should do? I saw you refer to that in this conversation Mounting from one cluster to another - #2 by Adrien.Albert[Mounting from one cluster to another - #2 by Adrien.Albert](https://hpc-community.unige.ch/t/mounting-from-one-cluster-to-another/3425/2)
I’m going to give it a try


## Post 11 by @Gael.Rossignol (2025-05-28T15:12:36.395Z)

Dear Nicolas,
I was thinking about the paragraph :
“# Transfer data from cluster to another with”
Let me know if helps you.
Best regards,


## Post 12 by @Yann.Sagon (2025-05-30T10:57:44.014Z)

I’ve added to the faq the link to the rsync how to as it may be helpful for other people : hpc:faq [eResearch Doc][hpc:faq [eResearch Doc]](https://doc.eresearch.unige.ch/hpc/faq#how_can_i_copy_data_from_one_c)
Best
