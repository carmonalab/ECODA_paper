# Baobab file writing on scratch from compute node very slow?

- Source: https://hpc-community.unige.ch/t/3767

- Created: 2024-12-16T17:16:42.029Z

- Posts: 8

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Quentin.Vagne (2024-12-16T17:16:42.073Z)

## Primary informations
Username: vagne
Cluster: baobab
Hi, I am currently running a code that write files after simulation steps. The files are not very big, they are typically 1.2M in size. For a reason that I don’t understand, my program sits for a very long time, using very minimal cpu, while progressively writing the file. It somehow takes several minutes for each file to be written which puzzles me since sending and receiving files from/to the login node is as fast as usual. I’m running my code in the scratch folder and writing the files there.
Is there any known problem that would explain why file writing from a compute node into the scratch folder is slow? Or should I look at my code to see if it’s doing something wrong?
Thanks for the help,
Quentin


## Post 2 by @Ludovic.Dumoulin (2024-12-17T12:58:13.969Z)

Hello Quentin,
The scratch storage has not been working properly for the past 4–6 months. If you need to save small data files, I recommend saving them to your home directory for now (even though it’s not ideal). Otherwise, you may experience significant slowdowns or crashes.
We are currently waiting for the HPC team to implement a permanent solution. The delay is due to the need to establish new usage rules for the scratch space, which is taking some time.
Hopefully, the situation will be resolved soon (before having to pay for simulations that crash mid-way due to storage issue).
Best,
Ludovic


## Post 3 by @Quentin.Vagne (2024-12-17T13:01:53.424Z)

Thanks for the reply!! Also now it seems that slurm is down? At least on my side?


## Post 4 by @Quentin.Vagne (2024-12-18T08:20:41.151Z)

Alright, slurm is back up, thanks to the person who fixed it!
And yeah, I moved my code to the normal home folder and file writing is clearly faster, even though it’s still quite slow.


## Post 5 by @Yann.Sagon (2024-12-18T15:24:53.222Z)

Ludovic.Dumoulin:
> The scratch storage has not been working properly for the past 4–6 months. If you need to save small data files, I recommend saving them to your home directory for now (even though it’s not ideal).
Thank you for your message. I’d just like to clarify that your message is not an official response. We (hpc team) are not aware that the scratch storage hasn’t been working properly for the last 4-5 months. What we have been facing is users filling up the storage and forgetting to clear all their temporary files. So yes, what has happened is that the storage has been full a couple of times, this isn’t the case yet as we have deleted many inactive users with their associated data (up to 30TB for some of them!). So we’ll be forced to enforce quotas next year.
Ludovic.Dumoulin:
> The delay is due to the need to establish new usage rules for the scratch space, which is taking some time.
We have three clusters to manage and plenty of other issues to solve as well, thanks for your patience.
Ludovic.Dumoulin:
> Hopefully, the situation will be resolved soon (before having to pay for simulations that crash mid-way due to storage issue)
We hope to get a lot of money to buy faster storage :rofl:


## Post 6 by @Yann.Sagon (2024-12-18T15:28:02.460Z)

Quentin.Vagne:
> Is there any known problem that would explain why file writing from a compute node into the scratch folder is slow? Or should I look at my code to see if it’s doing something wrong?
Can you please share your sbatch file so we can check which storage location you are using and what resource you use?
We remind you that you can as well request fast storage and use the local scratch.
https://doc.eresearch.unige.ch/hpc/storage_on_hpc#cluster_storage[https://doc.eresearch.unige.ch/hpc/storage_on_hpc#cluster_storage](https://doc.eresearch.unige.ch/hpc/storage_on_hpc#cluster_storage)


## Post 7 by @Quentin.Vagne (2025-02-19T10:14:12.409Z)

Hi,
Sorry for the extremely late answer. I didn’t see that there were new messages in the thread since I’m not regularly browsing the forum. (Side note: I looked at ways to enable to receive forum notifications by email, but I didn’t find any?).
I think that my answer now will probably be not so relevant anymore with the new policy on deleting scratch files. But basically my sbatch file is always something like that:
```
#!/bin/sh
#SBATCH --job-name XXX
#SBATCH --error out_err
#SBATCH --output out_stdout
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 3
#SBATCH --partition shared-cpu
#SBATCH --time 2:00:00

srun ./hlEmptyApp 
```
I was running a finite element code from my scratch folder “/home/users/v/vagne/scratch/”. It was writing small files (1.2MB) every simulation step, and the simulation steps were quite slow (so it was not trying to write 15 files per second). I took care to login on the compute node and check the cpu usage of my code, which was sitting at almost 0 percent during file writing, indicating that the program was limited by the file writing step. I also checked using “ls” and could see the file size slowly increase, showing that the writing process was slow, I guess?
I would assume ultimately, that at the time I was running the simulation, maybe the scratch storage was particularly full or heavily being read from/written on?
I will start running simulations again soon, after the current maintenance is over. I can try to run them on scratch again to see if the situation is different.
Thanks for the help.


## Post 8 by @Yann.Sagon (2025-02-20T13:13:54.629Z)

Dear Quentin,
Quentin.Vagne:
> Side note: I looked at ways to enable to receive forum notifications by email, but I didn’t find any
I’ve added a quick how to here: hpc:faq [eResearch Doc][hpc:faq [eResearch Doc]](https://doc.eresearch.unige.ch/hpc/faq#hpc_community_forum)
Quentin.Vagne:
> I took care to login on the compute node and check the cpu usage of my code
You can check with something like “htop” if the process appears with a  “D” state it is waiting on I/O.
Can you try to write the file to /scratch instead which is a local storage on the compute node. My suggestion is that you synchronize the file every N steps to your central scratch space? Anyway, maybe the issue was temporary, let us know if it is working as expected.
Best
