# [tutorial] Run a GitLab Runner that submits SLURM jobs on Baobab

- Source: https://hpc-community.unige.ch/t/4239

- Created: 2026-03-10T12:15:09.015Z

- Posts: 7

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yann.Sagon (2026-03-10T12:15:09.154Z)

# How to run a GitLab Runner that submits SLURM jobs
This document explains how to create a GitLab runner and configure it so that your CI/CD jobs are executed on the cluster using SLURM.
---
## 1. Create a token for registering a new GitLab runner
In your GitLab project, go to:
Settings → CI/CD → Runners → New project runner
Generate a new runner token:
image
image485×367 25.3 KB
[image485×367 25.3 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/dcfab1e2ab72db849c934295219e1a62d8cf996c.png)
image
image1872×161 21.3 KB
[image1872×161 21.3 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/77a515310fca711d1af09dbfe2a673a74e2ff207.png)
Specify a runner tag.
Jobs that specify this tag in `.gitlab-ci.yml` will be executed by this runner:
image
image1198×842 62.7 KB
[image1198×842 62.7 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/f071f9af360af59ccab0e2fb2e06210ce1745368.png)
Be sure to copy the token displayed at the end of the runner creation process (2).
This token is shown only once and will be required later:
image
image772×765 41.4 KB
[image772×765 41.4 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/9485583a262fa911ee762d92e87c65e9ec36e848.png)
---
## 2. Clone the `slurm-gitlab-executor` project
On login1, run:
```
git clone https://github.com/Algebraic-Programming/slurm-gitlab-executor.git slurm-gitlab-executor
cd slurm-gitlab-executor
```
---
## 3. Register your GitLab runner
From the cluster terminal, run (replace the token):
```
gitlab-runner register --name slurm-gitlab-executor --url https://gitlab.unige.ch --token glrt-<token-from-step-above>
```
---
## 4. Generate a new runner configuration
Run:
```
./generate-config.sh "/home/sagon/.gitlab-runner/config.toml" > config.toml
```
Then adjust the generated config file as needed. For example:
```
[[...]
    builds_dir = "/home/sagon/slurm-gitlab-executor/wd/builds"
    cache_dir = "/home/sagon/slurm-gitlab-executor/wd/cache"
[...]
```
---
## 5. Start the runner on login1
Run:
```
gitlab-runner run --config config.toml
```
---
## 6. Add a .gitlab-ci.yml to your project
In your GitLab repository, add a `.gitlab-ci.yml` file.
An example using SLURM is available here[here](https://gitlab.unige.ch/hpc/softs/-/blob/master/s/slurm-gitlab-executor/gitlab-ci.yml).
---
## 7. Launch a pipeline
In GitLab:
Build • Pipelines • New pipeline
You can now trigger a pipeline and configure SLURM parameters depending on your workflow.
Comments welcome!


## Post 2 by @Maroun.BouSleiman (2026-03-16T15:21:08.069Z)

Hi Yann,
Thanks a lot for the great tutorial.
Some of my CI jobs need to build container images on compute nodes using buildah inside Apptainer with --fakeroot. Currently Apptainer falls back to a root-mapped namespace because my user has no subuid/subgid entries:
INFO: User not listed in /etc/subuid, trying root-mapped namespace
This causes lchown failures during image layer extraction since only UID/GID 0 is mapped:
lchown /etc/gshadow: invalid argument (requested 0:42)
Request: Can you please add subuid/subgid entries for my user or all users on login and compute nodes?
# /etc/subuid and /etc/subgid
username:100000:65536
I saw you had already done something similar in Feb 2022 on login2.baobab (see forum thread about podman subuid errors: https://hpc-community.unige.ch/t/solved-cannot-login-to-login2/260[https://hpc-community.unige.ch/t/solved-cannot-login-to-login2/260](https://hpc-community.unige.ch/t/solved-cannot-login-to-login2/260)).
Thanks!
Maroun


## Post 3 by @Yann.Sagon (2026-03-18T12:56:49.908Z)

Dear @Maroun.BouSleiman[@Maroun.BouSleiman](https://hpc-community.unige.ch/u/maroun.bousleiman) thanks for the notification. We just figured out that it wasn’t done automatically by our scripts. This is now the case for subuid and subgid.
Best regards
Yann


## Post 4 by @Maroun.BouSleiman (2026-03-20T08:00:33.017Z)

Thank you Yann.
It works perfectly now. Docker builds are now possible on the HPC with buildah.
Cheers,
Maroun


## Post 5 by @Maroun.BouSleiman (2026-05-04T12:16:52.554Z)

Hi again Yann,
Is this setup specific for baobab? I was trying to get a runner working for builds on bamboo and encountering the same subuid and subgid issue. Didn’t try yggdarsil.
It would be great to be able to use another cluster when one is undergoing maintenance.
Thanks,
Maroun


## Post 6 by @Yann.Sagon (2026-05-04T13:13:53.223Z)

Hi,
Maroun.BouSleiman:
> Is this setup specific for baobab?
No, but it wasn’t yet applied as the change was planed for the next maintenance. I’ve forced the synchronization of the file now and it should be good.
Best regards
Yann


## Post 7 by @Maroun.BouSleiman (2026-05-04T13:25:19.966Z)

Thanks a lot Yann for the fast response!
