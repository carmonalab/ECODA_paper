# [HPC][baobab][slurm] Error allocating cpu

- Source: https://hpc-community.unige.ch/t/3974

- Created: 2025-06-06T12:04:45.882Z

- Posts: 5

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Mario.AlvesCardoso (2025-06-06T12:04:45.926Z)

Hey,
Username: alvescar
Cluster: private-dpnc-cpu,shared-bigmem,shared-cpu
Subject: failure allocating cpu
jobid: 3421 3472 3507 3508
I’m trying to allocate a cpu in baobab using the following command:
salloc -t 05:00:00 -p private-dpnc-cpu,shared-bigmem,shared-cpu --mem 128G
And I get the following message
image
image1077×186 30 KB
[image1077×186 30 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/7b69a84153cde8204d08f1be7f566959e1c26680.png)
I’ve tried changing the mem values and different clusters, but I get the same issue. I can still allocate jobs using sbatch but seems that I cant log to a cpu and can only stay in the login node.
I was wondering if there is any change I need to do to the salloc command to be able to log to a cpu node.
Best regards,
Mário


## Post 2 by @Adrien.Albert (2025-06-06T12:57:11.869Z)

Hi @Mario.AlvesCardoso[@Mario.AlvesCardoso](https://hpc-community.unige.ch/u/mario.alvescardoso)
Could you try again and let me know how it’s going ?


## Post 3 by @francois-xavier.meyer (2025-06-06T13:17:56.162Z)

Hi @Adrien.Albert[@Adrien.Albert](https://hpc-community.unige.ch/u/Adrien.Albert),
I had the same issue and now it is working for me.
Best regards,
François-Xavier


## Post 4 by @Mario.AlvesCardoso (2025-06-06T13:56:34.316Z)

Hi @Adrien.Albert[@Adrien.Albert](https://hpc-community.unige.ch/u/Adrien.Albert), it’s working now! Thanks so much.


## Post 5 by @Dimitrios.Proios (2025-06-06T16:00:04.628Z)

It is working thank you @Adrien.Albert[@Adrien.Albert](https://hpc-community.unige.ch/u/Adrien.Albert)
