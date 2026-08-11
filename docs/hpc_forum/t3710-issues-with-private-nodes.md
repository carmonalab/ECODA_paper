# Issues with private nodes?

- Source: https://hpc-community.unige.ch/t/3710

- Created: 2024-11-03T17:44:24.361Z

- Tags: yggdrasil

- Posts: 3

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Emeline.Bolmont (2024-11-03T17:44:24.411Z)

Hello dear HPC team,
I’m a bit puzzled as to what is happening to me. Here are the details below.
- what did I try:
- Here’s my launch script:
```
#!/bin/bash

#SBATCH --job-name=H2_1bar_095AU
#SBATCH --mail-type=START,END,FAIL
#SBATCH --mail-user=emeline.bolmont@unige.ch
#SBATCH --ntasks=24
#SBATCH --partition=private-astro-dynaclim-cpu
#SBATCH --time=60:00:00

source /home/bolmonte/2024_PCM/trunk/LMDZ.COMMON/arch/arch-BAO_YGG.env

srun 20241101_gcm_64x48x30_phystd_b40x35_para.e > gcm_64x48x30_phystd_b40x35_H2_1bar_H2Ovar_sma_095_AU_simu_0.out 2>&1 # parallele
```
- what didn’t work:
- It didn’t launch on the partition I asked for. It launches on debug-cpu.
- what was the expected result:
- It should have launched on the partition private-astro-dynaclim-cpu.
- what was the error message:
- In my .out file, I have the following:
```
srun: job 13520767 queued and waiting for resources
srun: job 13520767 has been allocated resources
srun: error: Nodes cpu001 are still not ready
srun: error: Something is wrong with the boot of the nodes.
```
- path to the relevant files (logs, sbatch script, etc):
- Everything is here: /home/bolmonte/2024_PCM/H2_simus/H2_0.95AU_1p
- My launch script: launch_script.sh
- The log : gcm_64x48x30_phystd_b40x35_H2_1bar_H2Ovar_sma_095_AU_simu_0.out
Thanks very much for your help!
Best,
Emeline Bolmont
PS: Do you now know when will my private nodes be migrated on Yggdrasil?


## Post 2 by @Emeline.Bolmont (2024-11-13T10:42:11.556Z)

I just wanted to say that somehow doing the following worked out. It’s not very satisfying, but it does the trick!
```
#!/bin/bash
source /home/bolmonte/2024_PCM/trunk/LMDZ.COMMON/arch/arch-BAO_YGG.env

srun --job-name=H2_1bar_095AU --mail-type=START,END,FAIL --mail-user=emeline.bolmont@unige.ch --ntasks=24 --partition=private-astro-dynaclim-cpu --time=60:00:00 20241101_gcm_64x48x30_phystd_b40x35_para.e > gcm_64x48x30_phystd_b40x35_H2_1bar_H2Ovar_sma_095_AU_simu_0.out 2>&1 # parallele
```


## Post 3 by @Yann.Sagon (2024-11-14T14:13:32.985Z)

Dear @Emeline.Bolmont[@Emeline.Bolmont](https://hpc-community.unige.ch/u/emeline.bolmont) ,
sorry, we missed this message. I’ve moved it to the category support issue which is more appropriate.
Emeline.Bolmont:
> It didn’t launch on the partition I asked for. It launches on debug-cpu.
This sounds like you submitted your job using `bash` and not `sbatch`. If this is the case, all the Slum pragma are treated as comment and have no effect.
Best
Yann
