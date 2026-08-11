# CPU Resource: Timeout after a few Hours?

- Source: https://hpc-community.unige.ch/t/3470

- Created: 2024-06-03T16:07:56.390Z

- Posts: 2

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Imahn.Shekhzadeh1 (2024-06-03T16:07:56.424Z)

Dear HPC Team,
I want to run a bash file on the HPC cluster for some longer time, and for this I request a CPU node like this:
```
salloc --partition=public-longrun-cpu --time=1-12:00:00 --nodes=1
```
However, after two hours, this happened:
```
(base) (baobab)-[shekhza2@cpu008 benchmark]$ ./pipeline.sh
client_loop: send disconnect: Broken pipe
```
This is not good, as my bash file now stopped executing. Can you please explain why this happens and how I can prevent this from happening?
Best,
Imahn Shekhzadeh


## Post 2 by @Imahn.Shekhzadeh1 (2024-06-07T11:41:21.012Z)

solved by creating `tmux` session
