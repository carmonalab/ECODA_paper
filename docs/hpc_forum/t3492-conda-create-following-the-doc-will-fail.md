# Conda create following the doc will fail?

- Source: https://hpc-community.unige.ch/t/3492

- Created: 2024-06-14T10:24:03.638Z

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Matthieu.Stigler (2024-06-14T10:24:03.681Z)

Hi
I am trying to create a conda environment, and I am following the steps in the HCP UNIGE documentation (hpc:applications_and_libraries [eResearch Doc][hpc:applications_and_libraries [eResearch Doc]](https://doc.eresearch.unige.ch/hpc/applications_and_libraries#conda)), I run:
```
module load Anaconda3
conda create --prefix $HOME/scratch/test-env --name environment_name
```
Running this, however, leads to the error:
> conda create: error: argument -n/–name: not allowed with argument -p/–prefix
It seems that indeed `--name` and `--prefix` cannot be used simultaneously? How should I do then? If I skip the `--name` part, the resulting name environment will be the whole prefix?


## Post 2 by @Matthieu.Stigler (2024-06-25T09:18:18.775Z)

Following up on this, it seems indeed the HPC documentation is wrong, looking at `conda create --help`, which returns:
> This command requires either the -n NAME or -p PREFIXoption
Could you kindly advice how I should do if I want to create an environment with name “name” yet on scratch?
