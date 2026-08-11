# Too many open files error

- Source: https://hpc-community.unige.ch/t/3388

- Created: 2024-03-22T09:16:57.868Z

- Posts: 3

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Debajyoti.Sengupta (2024-03-22T09:16:57.928Z)

Hi I am trying to run an export job (involves pytorch) with the following job specifications
```
#!/bin/sh
#SBATCH --job-name=export
#SBATCH --cpus-per-task=12
#SBATCH --time=1:00:00
#SBATCH --partition=private-dpnc-gpu,shared-gpu
#SBATCH --output=/home/users/s/senguptd/UniGe/astro/skycurtains/logs/%A_%a.out
#SBATCH --chdir=/home/users/s/senguptd/UniGe/astro/skycurtains/
#SBATCH --mem=15GB
#SBATCH --gpus=1
```
I get the following error `Too many open files. Communication with the workers is no longer possible.` But this seems to happen only on certain nodes gpu023, gpu024, gpu038. This worked on gpu040 for instance. Is the number of cpu requests too high for the former nodes and I need to turn those down?
What I can confirm is that within the code there are no files that are opened outside of a context manager.
Any help figuring this issue would be appreciated. If you require additional information, please let me know.
Thanks,
Deb


## Post 2 by @Adrien.Albert (2024-03-26T17:44:24.156Z)

Hi @Debajyoti.Sengupta[@Debajyoti.Sengupta](https://hpc-community.unige.ch/u/debajyoti.sengupta)
This problem does not concern nodes or software, but shared storage. At the time you needed to open/read a file, there were “Too many files open” at the same time.
This behavior is often due to a user (or several users) whose jobs perform a lot of I/O.
Can you tell me which jobID is affected? I’ll check the logs for any matches that might explain your problems.


## Post 3 by @Yann.Sagon (2024-03-27T07:56:44.114Z)

Debajyoti.Sengupta:
> I get the following error `Too many open files. Communication with the workers is no longer possible.`
My 2 cents: this error seems to be a known output/issue from PyTorch Too many open files error · Issue #11201 · pytorch/pytorch · GitHub[Too many open files error · Issue #11201 · pytorch/pytorch · GitHub](https://github.com/pytorch/pytorch/issues/11201) . As the `sbatch` you posted is incomplete, we don’t know which PyTorch you are running with what parameters. Please share that with us.
