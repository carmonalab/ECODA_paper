# BlockingIOError Yggdrasil

- Source: https://hpc-community.unige.ch/t/3693

- Created: 2024-10-21T07:15:45.035Z

- Tags: yggdrasil

- Posts: 7

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @felix.vecchi (2024-10-21T07:15:45.082Z)

Hi!
When I try to run a job on Yggdrasil I get an BlockingIOError which I cant sort out. Ive pasted the full error message below. The strange thing is that this error is not raised when I dont go via the slurm queue but just use ‘salloc’ to run it on a computing node. Does someone know what is going on?
Felix
```
Traceback (most recent call last):
  File "/srv/beegfs/scratch/users/f/fvecchi/Run/shear_fit_zs_LastShell_noise.py", line 204, in <module>
    result = search.fit(model=model, analysis=analysis)
             ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/users/f/fvecchi/cosma_home/Code/AllProjects/PyAutoFit/autofit/non_linear/search/abstract_search.py", line 594, in fit
    result = self.start_resume_fit(
             ^^^^^^^^^^^^^^^^^^^^^^
  File "/home/users/f/fvecchi/cosma_home/Code/AllProjects/PyAutoFit/autofit/non_linear/search/abstract_search.py", line 115, in decorated
    return func(self, *args, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/users/f/fvecchi/cosma_home/Code/AllProjects/PyAutoFit/autofit/non_linear/search/abstract_search.py", line 721, in start_resume_fit
    search_internal = self._fit(
                      ^^^^^^^^^^
  File "/home/users/f/fvecchi/cosma_home/Code/AllProjects/PyAutoFit/autofit/non_linear/search/nest/nautilus/search.py", line 149, in _fit
    search_internal = self.fit_multiprocessing(
                      ^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/users/f/fvecchi/cosma_home/Code/AllProjects/PyAutoFit/autofit/non_linear/search/nest/nautilus/search.py", line 265, in fit_multiprocessing
    return self.call_search(search_internal=search_internal, model=model, analysis=analysis)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/users/f/fvecchi/cosma_home/Code/AllProjects/PyAutoFit/autofit/non_linear/search/nest/nautilus/search.py", line 322, in call_search
    search_internal.run(
  File "/home/users/f/fvecchi/cosma_home/env3.11/lib/python3.11/site-packages/nautilus/sampler.py", line 448, in run
    self.write_shell_update(self.filepath, -1)
  File "/home/users/f/fvecchi/cosma_home/env3.11/lib/python3.11/site-packages/nautilus/sampler.py", line 1319, in write_shell_update
    fstream = h5py.File(Path(filepath), 'r+')
              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/users/f/fvecchi/cosma_home/env3.11/lib/python3.11/site-packages/h5py/_hl/files.py", line 561, in __init__
    fid = make_fid(name, mode, userblock_size, fapl, fcpl, swmr=swmr)
          ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/users/f/fvecchi/cosma_home/env3.11/lib/python3.11/site-packages/h5py/_hl/files.py", line 237, in make_fid
    fid = h5f.open(name, h5f.ACC_RDWR, fapl=fapl)
          ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "h5py/_objects.pyx", line 54, in h5py._objects.with_phil.wrapper
  File "h5py/_objects.pyx", line 55, in h5py._objects.with_phil.wrapper
  File "h5py/h5f.pyx", line 102, in h5py.h5f.open
BlockingIOError: [Errno 11] Unable to synchronously open file (unable to lock file, errno = 11, error message = 'Resource temporarily unavailable')
srun: error: cpu073: tasks 0-1: Exited with exit code 1
srun: error: cpu073: task 2: Exited with exit code 1
srun: First task exited 30s ago
srun: StepId=35948812.0 task 3: running
srun: StepId=35948812.0 tasks 0-2: exited abnormally
srun: Terminating StepId=35948812.0
slurmstepd: error: *** STEP 35948812.0 ON cpu073 CANCELLED AT 2024-10-18T14:23:47 ***
srun: Job step aborted: Waiting up to 92 seconds for job step to finish.
```


## Post 2 by @Yann.Sagon (2024-10-21T08:16:51.930Z)

Hi,
please share your sbatch with us.
Best
Yann


## Post 3 by @felix.vecchi (2024-10-21T08:38:59.998Z)

Sure!
```
#!/bin/bash -l

#SBATCH --ntasks 4
#SBATCH -J Nzs1.6
#SBATCH -o standard_output_file.%J.out
#SBATCH -e standard_error_file.%J.err
#SBATCH -p public-cpu
#SBATCH -t 24:00:00
#SBATCH --mail-type=END                          # notifications for job done & fail
#SBATCH --mail-user=<felix.vecchi@epfl.ch>

srun python shear_fit_zs_LastShell_noise.py 1.6 1
```


## Post 4 by @felix.vecchi (2024-10-21T08:40:07.876Z)

And before I sbatch this file I load python3.11 and my virtual environment


## Post 5 by @felix.vecchi (2024-10-21T09:43:10.261Z)

I think I posted my answer on the issue level and not in a direct reply to you. idk if that matters haha


## Post 6 by @Yann.Sagon (2024-10-21T09:49:33.941Z)

Hi,
You are requesting 4 tasks in your sbatch and your python code is probably not aware of how to handle this and so executed 4 times, which explains the blocking IO.
You should check our doc to see the different job types: hpc:slurm [eResearch Doc][hpc:slurm [eResearch Doc]](https://doc.eresearch.unige.ch/hpc/slurm#single_thread_vs_multi_thread_vs_distributed_jobs)
In your case it is probably either single core or multi core but not distributed.
Best


## Post 7 by @Yann.Sagon (2024-10-21T09:52:33.971Z)

felix.vecchi:
> I think I posted my answer on the issue level and not in a direct reply to you.
I think the result is the same. Anyway, remember you can still edit your post if needed, for example to add
felix.vecchi:
> And before I sbatch this file I load python3.11 and my virtual environment
To your previous post wher you sent your sbatch.
