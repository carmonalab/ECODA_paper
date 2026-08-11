# Problem launching Python job on Baobab

- Source: https://hpc-community.unige.ch/t/3400

- Created: 2024-04-03T09:55:59.943Z

- Posts: 4

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Monika.Avila (2024-04-03T09:56:00.036Z)

## Primary informations
Username: avilamar
Cluster: Baobab
## Description
I want to submit a Python job on Baobab using a batch file and the command sbatch, but I get the following error:
Running main.py on cpu001.baobab
slurmstepd: error: execve(): Python: No such file or directory
srun: error: cpu001: task 0: Exited with exit code 2
## Steps to Reproduce
## My batch file is called main.sh and contains:
#!/bin/bash
#SBATCH --mem=3G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --licenses=matlab@matlablm.unige.ch
#SBATCH --time=00:15:00
#SBATCH --mail-user=monika.avila@unige.ch
#SBATCH --partition=debug-cpu
module load GCCcore/8.3.0
module load Python
BASE_MFILE_NAME=main
unset DISPLAY
echo “Running ${BASE_MFILE_NAME}.py on $(hostname)”
## srun Python ${BASE_MFILE_NAME}.py
## Expected Result
I expect to obtain a job that is running
## Actual Result
I get an error
Running main.py on cpu001.baobab
slurmstepd: error: execve(): Python: No such file or directory
srun: error: cpu001: task 0: Exited with exit code 2
I also replaced the las line of my batch file by
srun ${BASE_MFILE_NAME}.py
And in that case I get the mistake
Running main.py on cpu001.baobab
slurmstepd: error: execve(): main.py: Permission denied
srun: error: cpu001: task 0: Exited with exit code 13


## Post 2 by @Gael.Rossignol (2024-04-03T11:44:45.486Z)

Monika.Avila:
> Python
Dear Monika,
It seems you just have a typo in your sbatch scripts. You call “Python” but the software name is “python” lowercase.
Best regards,


## Post 3 by @Monika.Avila (2024-05-13T09:45:33.156Z)

Thanks, Gael! You are very kind for taking the time to answer my post. However, it was not the upper or lower case that made a difference in my case.  Actually, I solved the issue with the following line:
module load GCC/8.2.0-2.31.1 OpenMPI/3.1.3 Python/3.7.2 SciPy-bundle/2019.03


## Post 4 by @Adrien.Albert (2024-05-13T11:22:38.940Z)

Monika.Avila:
> module load GCC/8.2.0-2.31.1 OpenMPI/3.1.3 Python/3.7.2 SciPy-bundle/2019.03
Hi Monika,
the module command lets you manage the software versions you want to use. In your case, you’re only loading a different version of python installed by default on the server.
According to the following tests, it’s `python` in lower case that works (regardless of the version used)
```
## check command
(baobab)-[alberta@login2 ~]$ Python --version
-bash: Python: command not found
(baobab)-[alberta@login2 ~]$ python --version
Python 3.6.8
(baobab)-[alberta@login2 ~]$ module load GCC/8.2.0-2.31.1 OpenMPI/3.1.3 Python/3.7.2 SciPy-bundle/2019.03
(baobab)-[alberta@login2 ~]$ Python --version
-bash: Python: command not found
(baobab)-[alberta@login2 ~]$ python --version
Python 3.7.2

## Launching job on debug-cpu
(baobab)-[alberta@login2 ~]$ srun Python --version
srun: job 9927929 queued and waiting for resources
srun: job 9927929 has been allocated resources
slurmstepd: error: execve(): Python: No such file or directory
srun: error: cpu001: task 0: Exited with exit code 2
(baobab)-[alberta@login2 ~]$ srun python --version
srun: job 9927931 queued and waiting for resources
srun: job 9927931 has been allocated resources
Python 3.7.2
```
Best Regards
