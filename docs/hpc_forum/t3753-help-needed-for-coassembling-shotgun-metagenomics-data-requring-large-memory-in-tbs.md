# Help needed for coassembling Shotgun metagenomics data requring large memory (in TBs)

- Source: https://hpc-community.unige.ch/t/3753

- Created: 2024-12-02T13:14:09.066Z

- Tags: baobab

- Posts: 7

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Jaspreet.Saini (2024-12-02T13:14:09.109Z)

Dear HPC,
My username is sainij. I have around 1TB of microbial datasets collected from the gut of participants associated with Alzheimer disease. This dataset is generated from Shotgun DNA sequencing technology, and I would like to make effort to process these samples together. However, this step is intensive and requires crazy amount of RAM because of the size of the datasets.
My initial strategy if to use tools (example MEGAHIT[MEGAHIT](https://academic.oup.com/bioinformatics/article/31/10/1674/177884)) which had low RAM requirements. I initially used 500GB of ram with 8 CPUs per task using shared-bigmem, however, I still had out of memory issue. I tried to use public-bigmem to give higher RAM, but it seems to be down?
I kindly ask you for your advice on the matter. Thank you for your assistance.
Kind regards,
Jaspreet
Following is the script I am working;
```
#!/bin/bash
#SBATCH -e M-%A-error
#SBATCH -o M-%A-out
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH --partition=public-bigmem
#SBATCH --mem 995GB
#SBATCH --mail-type=ALL
#SBATCH --time=48:00:00
#SBATCH --mail-user=jaspreet.saini@unige.ch

megahit -1 Gmad_cat_R1.BBnorm.fastq.gz -2 Gmad_cat_R2.BBnorm.fastq.gz  -o megahit_output_gmad -t 16 --continue
```


## Post 2 by @Adrien.Albert (2024-12-03T21:24:27.766Z)

Hi @Jaspreet.Saini[@Jaspreet.Saini](https://hpc-community.unige.ch/u/jaspreet.saini)
Which cluster are you trying to run your job on?
Have you checked that your partition has enough memory to meet your requirements?
On Boabab, the nodes are quite old and some nodes have more memory than the public bigmem partition. The bigmem partition certainly no longer lives up to its name. :no_mouth:
```
(baobab)-[alberta@login1 ~]$ sinfo -p public-bigmem -o %n,%m
HOSTNAMES,MEMORY
cpu246,224000
cpu245,256000
```
But let’s see the shared-bigmem partition:
```
(baobab)-[root@login1 alberta]$ sinfo -p shared-bigmem -o  %n,%m
HOSTNAMES,MEMORY
cpu246,224000
cpu154,256000
cpu203,512000
cpu218,512000
cpu219,512000
cpu245,256000
cpu312,1024000    <========= 1TB
cpu313,1024000    <========= 1TB
```
So My following jobs is working on shared-bigmem:
```
(baobab)-[root@login1 alberta]$ srun --mem=995GB -p shared-bigmem --time=00:00:02 hostname 
srun: job 13894880 queued and waiting for resources
^Csrun: Job allocation 13894880 has been revoked
```
However, you are limited to 12 Hours
BUT on bamboo you have public partition (and so 4 days walltime):
```
(bamboo)-[root@admin1 ~]$ sinfo -p public-bigmem -o %n,%m
HOSTNAMES,MEMORY
cpu044,1024000
cpu045,1024000
```
I didn’t check on Yggdrasil :slight_smile:


## Post 3 by @Jaspreet.Saini (2024-12-05T12:34:21.591Z)

Dear Adrien,
Thank you for getting in touch and for useful insights. I tried to schedule job with shared-bigmem using 995GB memory. However, the job seems to be stuck  `         13897446 shared-bi megahit.   sainij PD       0:00      1 (ReqNodeNotAvail, UnavailableNodes:cpu[218-219,246,312-313])`.   Would you suggest me to move to Bamboo once it gets online? Thank you.


## Post 4 by @Gael.Rossignol (2024-12-06T08:47:57.716Z)

Hello,
Bamboo is now online and two bigmem nodes (with 1Tb RAM) are available on partition : shared-bigmem
In addition I populate MEGAHIT software on Bamboo using Easybuild, to facilitate your migration :
```
(base) (bamboo)-[rossigng@login1 ~]$ ml GCCcore/12.3.0 MEGAHIT/1.2.9
(base) (bamboo)-[rossigng@login1 ~]$ megahit --version
MEGAHIT v1.2.9
```
Best regards,


## Post 5 by @Jaspreet.Saini (2024-12-06T16:01:19.210Z)

Super! Thank you so much.


## Post 6 by @Jaspreet.Saini (2024-12-12T11:37:57.909Z)

Follow up on the task requiring intensive memory,
I ran a job using the shared-bigmem partition on Yggdrasil, and MEGAHIT utilized approximately 1.33TB of memory. Below are the memory consumption details:
```
seff 36974920  
perl: warning: Setting locale failed.  
perl: warning: Please check that your locale settings:  
  LANGUAGE = (unset),  
  LC_ALL = (unset),  
  LC_CTYPE = "UTF-8",  
  LANG = "en_US.UTF-8"  
perl: warning: Falling back to a fallback locale ("en_US.UTF-8").  

Job ID: 36974920  
Cluster: yggdrasil  
User/Group: sainij/hpc_users  
State: TIMEOUT (exit code 0)  
Nodes: 1  
Cores per node: 32  
CPU Utilized: 00:00:00  
CPU Efficiency: 0.00% of 16-00:11:12 core-walltime  
Job Wall-clock time: 12:00:21  
Memory Utilized: 1.33 TB  
Memory Efficiency: 90.64% of 1.46 TB  
```
Unfortunately, the job was terminated due to the time limit. The `--continue` option in MEGAHIT is not very effective as it restarts the job from the beginning instead of resuming from the last checkpoint. The following logs highlight this issue:
```
Continue mode activated. Ignore all options except for -o/--out-dir.  
2024-12-09 13:20:26 - MEGAHIT v1.2.9  
2024-12-09 13:20:26 - Using megahit_core with POPCNT and BMI2 support  
2024-12-09 13:20:26 - passing checkpoint 0  
2024-12-09 13:20:26 - passing checkpoint 1  
2024-12-09 13:20:26 - k-max reset to: 141  
2024-12-09 13:20:26 - Start assembly. Number of CPU threads: 16  
2024-12-09 13:20:26 - k list: 21,29,39,59,79,99,119,141  
2024-12-09 13:20:26 - Memory used: 1459801862553  
2024-12-09 13:20:26 - Extract solid (k+1)-mers for k = 21  
slurmstepd: error: *** JOB 36974920 ON cpu120 CANCELLED AT 2024-12-10T01:20:46 DUE TO TIME LIMIT ***  
```
I also checked the public-bigmem partition, but it supports only up to 770GB of memory, which is insufficient for my requirements:
```
sinfo -p public-bigmem -o %n,%m,%cpu  
HOSTNAMES, MEMORY, CPUS  
cpu115, 770000, 16  
cpu112, 770000, 16  
cpu114, 770000, 16  
cpu113, 770000, 16
```
Could you please advise on how I can proceed to run this job successfully? I require access to a partition that supports up to 1.5TB of memory with atleast 24-48 hours of run time.
Thank you in advance for your assistance.


## Post 7 by @Yann.Sagon (2024-12-18T15:56:53.602Z)

Jaspreet.Saini:
> `CPU Utilized: 00:00:00`
This is weird. It seems your job is using CPU very poorly, maybe this is why 12h isn’t enough? Is the job I/O bounded? I meant lots of read/write on disks? If yes, try to use local /scratch that are faster for temporary files?
Jaspreet.Saini:
> I require access to a partition that supports up to 1.5TB of memory with atleast 24-48 hours of run time.
You need to ask the owner of the private partition `private-wesolowski-bigmem` if you can use their nodes.
