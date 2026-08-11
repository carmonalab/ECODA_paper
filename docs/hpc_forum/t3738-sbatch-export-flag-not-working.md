# `sbatch --export` flag not working?

- Source: https://hpc-community.unige.ch/t/3738

- Created: 2024-11-21T11:04:13.418Z

- Posts: 5

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Berk.Gercek (2024-11-21T11:04:13.460Z)

I’ve been attempting to use the `--export` flag with `sbatch` to specify per-run environment variables for individual scripts. The function, however, does not seem to be working as it is specified in the SLURM manual for that command[SLURM manual for that command](https://slurm.schedmd.com/sbatch.html).
A minimal example (ran on `debug-cpu` on Bamboo) is:
script.sh
```
#!/bin/bash
echo ${A}
echo ${B}
```
I then, in a terminal with `A` and `B` unset, run the following command per the man page:
```
$ sbatch script.sh --export=A=1,B=2
```
The expected output is that the file `slurm-XXXX.out` should contain 1 and 2, but instead the file is empty with no contents.
I have found that running the following command instead will work:
```
$ export A=1
$ sbatch ./script.sh --export=A=1
```
but trying to change the value in the `--export` flag from 1 to another value will retain the original exported value of 1.
The command should work as specified in my first example according to the documentation, allowing me to bind a variable using the export flag only. Am I missing something, or is this a bug?


## Post 2 by @Yann.Sagon (2024-11-21T12:38:43.066Z)

Dear @Berk.Gercek[@Berk.Gercek](https://hpc-community.unige.ch/u/berk.gercek)
Berk.Gercek:
> I then, in a terminal with `A` and `B` unset, run the following command per the man page:
> `$ sbatch script.sh --export=A=1,B=2`
The issue here is the argument order: you are specifying the flag `--export` as argument for your script `script.sh` which is not what you want. You need to do something like the following:
```
$ sbatch --export=A=1,B=2 script.sh 
```
By the way, are you sure you need to that? By default all the users’s env var are propagated. If you specify a variable to export using `--export` then no other variables are exported.


## Post 3 by @Berk.Gercek (2025-01-12T11:39:08.547Z)

Hey Yann,
Thank you for the response! I was planning on using the env variables to specify certain parameters to the job script that change between simultaneously submitted jobs. I.e SUBJECT=“A” for one and SUBJECT=“B” for another. It seemed like an easy way to approach the problem, though it seems to be difficult in practice.
In the example you provided the issue is that B is still not passed on to the job script environment. A is successfully propagated but B is empty when echoed. Do you have an idea of why that might be?
Berk


## Post 4 by @Yann.Sagon (2025-01-13T07:40:00.720Z)

Berk.Gercek:
> I was planning on using the env variables to specify certain parameters to the job script that change between simultaneously submitted jobs. I.e SUBJECT=“A” for one and SUBJECT=“B” for another.
Then, no need to use the export flag. As I said all then env variables are propagated by default. See the example below.
```
(baobab)-[sagon@login1 ~]$ SUBJECT=A sbatch --wrap "srun env | grep SUBJECT"
Submitted batch job 14275216
(baobab)-[sagon@login1 ~]$ cat slurm-14275216.out
SUBJECT=A
```
The previous solution I wrote is working too, A and B are correctly exported.
```
(baobab)-[sagon@login1 ~]$ sbatch --export=A=1,B=2 --wrap "srun /usr/bin/env | grep -e 'B=' -e 'A='"
Submitted batch job 14275228
(baobab)-[sagon@login1 ~]$ cat slurm-14275228.out
A=1
B=2
SLURM_EXPORT_ENV=A=1,B=2
```
But better avoid to use export flag as it removes other important env variables.


## Post 5 by @Berk.Gercek (2025-01-19T09:26:18.920Z)

I’ll find a workaround, probably using an env file defining job constants! Thanks.
