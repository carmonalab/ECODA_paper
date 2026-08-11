# Maximum number of jobs

- Source: https://hpc-community.unige.ch/t/3514

- Created: 2024-06-25T15:17:11.608Z

- Posts: 12

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Mara.Attia (2024-06-25T15:17:11.672Z)

## Primary informations
Username: oattia
Cluster: baobab
## Description
I need to perform ~ 100k jobs, which are variations of one job but with different parameters. I know the preferred way to do this would be to use job arrays, but I already coded it differently before I knew the existence of job arrays. The way I am currently doing it is to launch a “master” bash script, which loops over individual bash scripts, each one of them launching a job.
I checked the maximum number of jobs I could simultaneously have in the queue by doing `scontrol show config | grep 'MaxJobCount`’. The result is 60k, but in reality I can only launch 10k at the same time. Indeed, after launching my master bash script and waiting for the individual scripts to be launched, running `squeue | grep -c oattia` stays blocked on 9999. What is the reason of this discrepancy?
Plus, according to the slurm documentation, using job arrays won’t make a difference because a job array doesn’t count as just one job, but as many jobs defined in the array (so 10k tops, as returned by `scontrol show config | grep 'MaxArraySize'` I suppose). What is the benefit of using job arrays if they cannot bypass MaxJobCount, aside from good practice and readability (and I/O issues related to having many individual bash scripts)?
## Steps to Reproduce
Master bash script:
```
#! /bin/bash
#SBATCH --partition=private-astro-cpu
declare -a StringArray=("file_1.bash" "file_2.bash" [...] "file_100000.bash")
for script in "${StringArray[@]}"; do
    sbatch $script
done
```
Every file (`file_k.bash`) is a pretty basic one, running a python script with a specified set of parameters.
## Expected Result
Bash files from 1 to 60000 should run, as specified by `MaxJobCount`.
## Actual Result
Only the first 9999 actually run. The others are just ignored.


## Post 2 by @Mara.Attia (2024-06-26T13:00:27.919Z)

Hi,
I did some investigation. If I understand correctly, the 10k job limit actually comes from the QOS. Indeed, no `GrpJobs`, `MaxJobs`, `GrpSubmit`, or `MaxSubmit` are defined at the association/account level with:
```
sacctmgr show associations user=oattia
```
On the other hand, the previous command tells me my QOS is `normal`, and typing:
```
sacctmgr show qos name=normal
```
reveals that `MaxSubmitPU` is set to `10000`. So my question now is: is there a way to bypass this limit? How do I launch, e.g., 100k jobs without doing it manually 10 times for blocks of 10k jobs after waiting for the completion of each block (either naively or with a job array, the result should be the same with respect to the QOS limit)?
Many thanks.


## Post 3 by @Yann.Sagon (2024-06-28T10:23:25.295Z)

Mara.Attia:
> What is the benefit of using job arrays if they cannot bypass MaxJobCount, aside from good practice and readability (and I/O issues related to having many individual bash scripts)?
Handling job arrays is indeed good practice, only one script to write, easier to manage: you can cancel the job array etc.
Mara.Attia:
> reveals that `MaxSubmitPU` is set to `10000`.
This is done on purpose as having to many jobs in the queue is an issue for Slurm (scheduling and memory). So it is up to you to find another way.
Please let us know what is the duration of each of your job and what are the resources needed so we may suggest you other solutions.
Best


## Post 4 by @Mara.Attia (2024-07-03T09:55:34.476Z)

Hi, thank you for your answer.
In terms of duration, the median wall time is around 5 minutes per job. However, they can strongly vary from just a few seconds to a couple of days (the simulations are chaotic). In terms of resources, I need one CPU per job.
Do you have a viable solution for me? What is the best way to handle high throughput computations?
Cheers.


## Post 5 by @Yann.Sagon (2024-07-03T11:31:57.694Z)

Maybe you should request resources for one big job (for example 50 cpus) and then you use parallel[parallel](https://www.gnu.org/software/parallel/) to distribute the load, having always 50 running jobs.


## Post 6 by @Mara.Attia (2024-07-04T13:50:56.126Z)

Thanks for your answer. I see `parallel` is not installed on baobab, but `xargs` is. Would the following script work for launching 100k jobs, 50 by 50 (supposing there are 100k input files in `/path/to/input/`)?
```
#!/bin/bash
#SBATCH --job-name=test
#SBATCH --cpus-per-task=50
#SBATCH --partition=private-astro-cpu
#SBATCH --hint=nomultithread

find /path/to/input/*.txt -print0 | xargs -0 -n1 -P 50 srun -n1 --exclusive python -u python_script.py
```
Also, how to generate one distinct .out and .err per job (with the same name as the parsed input file)?
Finally, how are exceptions handled? If a job crashes, will it cancel the following ones?
Many thanks.


## Post 7 by @Yann.Sagon (2024-07-04T14:01:01.359Z)

Mara.Attia:
> I see `parallel` is not installed on baobab
You can still use the one we provide through easybuild:
```
(baobab)-[sagon@login1 ~]$ ml GCCcore/11.2.0 parallel/20210722
```
Mara.Attia:
> Would the following script work for launching 100k jobs, 50 by 50
Do not forget to request enough time to run the 100k jobs! If this is too big, split the job.
Mara.Attia:
> Also, how to generate one distinct .out and .err per job
It seems you can add `--output` to `srun` to specify the log file.


## Post 8 by @Mara.Attia (2024-07-04T15:42:20.998Z)

Thanks. I will try the following, which seems to be doing what I want:
```
#!/bin/bash
#SBATCH --job-name=test
#SBATCH --cpus-per-task=50
#SBATCH --partition=private-astro-cpu

xargs -0 -n1 -P 50 -I{} -a input.inp -d '\n' \
srun -n1 -c1 -o {}.out -e {}.err --exclusive \
python -u python_script.py /path/to/input_folder/{}.txt
```
Here, `input.inp` contains the names of my input files (one per line) without any extension, so that I can easily append any extension using `-I{}`.
Before testing the above, small question about resources allocation. Will that script indeed allocate 50 cpus (via `#SBATCH --cpus-per-task=50`), and launch the jobs 50 by 50 (via `xargs  -P 50`) with one cpu per job (via `srun -c1`)? I am not sure to know how to actually check that.


## Post 9 by @Yann.Sagon (2024-07-05T07:13:59.358Z)

Mara.Attia:
> Will that script indeed allocate 50 cpus (via `#SBATCH --cpus-per-task=50`), and launch the jobs 50 by 50 (via `xargs  -P 50`) with one cpu per job (via `srun -c1`)? I am not sure to know how to actually check that.
Technically, an srun inside an sbatch is a step[step](https://slurm.schedmd.com/job_launch.html#step_allocation). After re reading this documentation, it appears that many options supported by srun when allocating resources are not supported by srun when launching step. Thus, you have to give a try, note sure for example that the option `-e` would work.
Once your job is running (you can always start with less CPUs and use the debug partition for this purpose, you can connect to the allocated compute node using ssh, then you can launch `htop` to see if you jobs are running as you would expect.


## Post 10 by @Mara.Attia (2024-07-05T10:04:17.192Z)

So I tested this solution, and it works exactly as it should.
Here is the bash script in question:
```
#!/bin/bash
#SBATCH --job-name=test
#SBATCH --cpus-per-task=3
#SBATCH --partition=debug-cpu
#SBATCH --time=00:15:00

xargs -0 -n 1 -P 3 -I {} -a input.inp -d '\n' \
srun -n 1 -c 1 -o {}.out -e {}.err --exclusive \
python -u test.py {}
```
The python script `test.py`:
```
import sys
import psutil
import time

if __name__ == '__main__':
    arg = float(sys.argv[1])
    print('The square of {} is {}.'.format(arg, arg**2))
    print('My cpu id is {:d}.'.format(psutil.Process().cpu_num()))
    time.sleep(60)
    print('End.')
```
The input file `input.inp`:
```
01
02
03
04
exception
06
07
08
09
10
11
12
```
The job runs by blocks of 3 steps as anticipated, taking a wall-time of 4 minutes (as each step has to `time.sleep(60)` for a minute, and there are 12 steps).
All the output/error files are correctly created (`01.out`, `01.err`, etc).
The `cpu id` shown in the `.out` files are distinct (either `0`, `1` or `2`) within each block of 3 steps, as it should.
The forced error (`exception` within the input file) is reported in `exception.err` as well as in the wrapper job output `slurm-11200736.out` (`srun: error: cpu002: task 0: Exited with exit code 1`). It does not prevent the launch of any of the other steps.
All good, many thanks for your help!


## Post 11 by @Yann.Sagon (2024-07-05T11:32:19.916Z)

Hi, well done, happy that it worked as expected!
Best


## Post 12 by @Matthias.Kruckow (2024-08-28T14:11:31.454Z)

Just a remark: at some point your list in the bunch of jobs run in the master job is empty and your overall usages of the requested resources will drop from that point until it converges towards 1 script running, while there there is a lot more reserved for your SLURM job. Hence, I suggest you to find a good balance on how many resources your wrapping job should request to avoid having too many of your assigned resources being idle on average.
