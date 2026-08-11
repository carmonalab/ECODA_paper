# Jupyter kernel memory limit?

- Source: https://hpc-community.unige.ch/t/3516

- Created: 2024-06-27T12:29:41.482Z

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Max.Briel (2024-06-27T12:29:41.524Z)

Hi!
I’m trying to run a Jupyter notebook server on an interactive node and am running into issues with the memory/kernel dieing.
I’m able to setup the linking to the correct node and ssh tunnel the session to my local machine and open it in the browser or visual studio code.
I’ve submitted the session with a single cpu and 4GB of RAM.
However, when I try to load a significant amount of data (~500MB), the kernel dies, but the slurm job stays active.
I am able to load smaller amount of data and run non-data-loading cells.
Since the kernel dies but the slurm job remains active, I was wondering if there’s an additional memory limit for Jupyter notebooks?
If so, can this limit be increased?


## Post 2 by @Gael.Rossignol (2024-06-27T14:06:59.625Z)

Max.Briel:
> If so, can this limit be increased?
Hello,
You have to check with “seff” with job id if the memory is not too light for this kind of job.
Best regards,
