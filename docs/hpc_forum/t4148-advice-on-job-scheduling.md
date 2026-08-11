# Advice on job scheduling

- Source: https://hpc-community.unige.ch/t/4148

- Created: 2025-11-18T10:14:58.814Z

- Tags: bamboo

- Posts: 4

- Category: 13

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Anna.Ramisch (2025-11-18T10:14:58.875Z)

Dear HPC team,
I have a question regarding job scheduling for a large analysis. I am running a VEP (Variant Effect Predictor) analysis across all chromosomes, which is quite time-consuming. To improve throughput, I split each chromosome into many smaller chunks. This results in ~20,000 jobs, each requesting 4 hours of runtime and 15 GB of memory. I then rerun the small number of failed chunks (typically <5%) with increased resources.
Each chunk contains the same number of variants, but the processing time varies quite a lot, so assigning a single reliable resource request to all chunks is difficult. Ideally, I would run one job per chromosome, but that would require a large amount of memory, many CPU cores, and long runtimes, and would likely result in extremely long queueing times (if even possible at all).
From the information on the website, it seems that renting private compute nodes is only possible for a minimum period of six months, right? Since I only need this full analysis roughly once per year, renting a node for that long may not be practical. With that in mind, I wanted to ask whether you might have any recommendations or best practices for running this kind of workload more efficiently on the existing clusters.
I also wanted to check whether there is a maximum number of jobs that can be submitted to the scheduling system at once, so that I remain within any submission or array limits.
Thank you very much for your help and any advice you may have.
Best regards,
Anna


## Post 2 by @Anna.Ramisch (2025-11-19T09:59:56.219Z)

Dear HPC team,
I have a follow-up question. After my post, I just tried to run the smallest chromosome on the public-cpu partition, and the job actually finished in about 14 hours. This is surprising because earlier I processed the data in chunks of 100,000 variants on the shared-cpu partition (the full chromosome has ~24.5 million variants), and some of those chunks took more than 2 hours each. There is no substantial preparation step in the code, so that shouldn’t explain the difference.
Additionally, when I reran the same slow chunks, they sometimes finished in only 30 minutes.
Do you have any idea what might cause such large variability in run time (both chunking vs. not chunked, but also for the same file)? I can share the code if needed, but I’m wondering whether there is a general explanation for this behavior. Due to chunking, I ended up spending much much more total CPU time than expected, 20 times more (rough estimate), which seems unusual, and unfortunately very expensive :confused: Since I’m not a computer scientist, maybe there’s a very obvious answer to this that I just can’t think of.
Thank you already,
Anna


## Post 3 by @Gael.Rossignol (2025-11-27T14:27:29.213Z)

Dear Anna,
I feel that the issue in this is the storage used.
Where are stocked files? How many access does it needs?
Sometime the more efficient process is to copy the file on the local scratch (`/scratch`) of the node, process files in the same space (`/scratch`) then copy result on the shared storage (scratch /srv/beegfs/scratch or home /home).
Can you test to do this and check if the time of processing remains stable?
Best regards,


## Post 4 by @Anna.Ramisch (2025-12-11T12:27:35.583Z)

Thank you for your answer, and sorry for not replying earlier - I totally missed it.
I found a way with way fewer batches, so I don’t have the same issue for now. But it makes sense what you’re saying. From what I observed the more of the chunked jobs were running on the same cpu, the slower it got. The VEP tool has to more or less permanently read 500 Gb files located on shared storage, so I guess when this happens a lot in parallel, this could explain it.
I will still see if the solution with the local scratch can be helpful for the rest of my jobs. Thanks again!
All the best,
Anna
