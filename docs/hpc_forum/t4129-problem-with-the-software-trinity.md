# Problem with the software Trinity

- Source: https://hpc-community.unige.ch/t/4129

- Created: 2025-10-24T12:21:34.304Z

- Posts: 8

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Nathan.Evan (2025-10-24T12:21:34.388Z)

Hello, I am currently trying to assemble RNA-seq data with the software Trinity. I got an error related to samtools installation. Here is the script i use :
#!/bin/bash
#SBATCH --job-name trinity_dpla            # this is a parameter to help you sort your job when listing it
#SBATCH --mem 100G
#SBATCH --cpus-per-task 12
#SBATCH --partition shared-cpu         # the partition to use. By default debug-cpu
#SBATCH --time 12:00:00                  # maximum run time.
ml GCC/8.3.0 SAMtools/1.10 OpenMPI/3.1.4 Trinity/2.9.1-Python-3.7.4
srun Trinity --seqType fq
–left “Rep1_1.sp.fastq.gz,Rep2_1.sp.fastq.gz”
–right “Rep1_2.sp.fastq.gz,Rep2_2.sp.fastq.gz”
–CPU 12
–max_memory 90G
–output /home/users/e/evann/scratch/01_Trinity_Assembly.
and here is the error i get
Trinity version: Trinity-v2.9.1
-ERROR: couldn’t run the network check to confirm latest Trinity software version.
Argument “error” isn’t numeric in numeric ne (!=) at /opt/ebsofts/Trinity/2.9.1-foss-2019b-Python-3.7.4/trinityrnaseq-v2.9.1/Trinity line 3865.
Error, need samtools installed that is at least as new as version 1.3 at /opt/ebsofts/Trinity/2.9.1-foss-2019b-Python-3.7.4/trinityrnaseq-v2.9.1/Trinity line 3866.
srun: error: cpu331: task 0: Exited with exit code 127
I didn’t find a solution to this. Can somebody help me?
Thank you in advance


## Post 2 by @Yann.Sagon (2025-10-31T12:52:02.684Z)

Dear @Nathan.Evan[@Nathan.Evan](https://hpc-community.unige.ch/u/nathan.evan)
I’m compiling a recent version of Trinity as the one is using outdated dependencies. I’ll let you know once done, it takes some time because there is a lot of dependencies such as R.


## Post 3 by @Nathan.Evan (2025-11-03T13:18:21.364Z)

Thank you very much @Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon) ^^


## Post 4 by @Nathan.Evan (2025-11-06T10:56:35.294Z)

Dear @Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon)
I understand you probably have a lot to do, but we need this analysis relatively urgently. Do you have a time frame for when it will be ready ?


## Post 5 by @Yann.Sagon (2025-11-06T16:16:32.576Z)

Hi, I had a couple of issue compiling this software as it has a bunch of dependencies and I was stuck when compiling the ~1300 extensions of R. Finally I found out that the compute node I was using to do the compilation wasn’t the one I expected (i.e. not the correct CPU, incompatible with some other nodes).
Well, ongoing process!


## Post 6 by @Nathan.Evan (2025-11-07T10:31:57.430Z)

Okey, thank you very much and please tell me when I can run my analysis ^^


## Post 7 by @Yann.Sagon (2025-11-07T14:34:52.586Z)

Please check: New software installed: Trinity version 2.15.2[New software installed: Trinity version 2.15.2](https://hpc-community.unige.ch/t/new-software-installed-trinity-version-2-15-2/4141)


## Post 8 by @Nathan.Evan (2025-11-07T15:26:50.872Z)

Thank you so much ! Have a nice weekend
