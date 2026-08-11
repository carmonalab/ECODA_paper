# Bamboo GPUs: Drained / Draining

- Source: https://hpc-community.unige.ch/t/4350

- Created: 2026-07-17T10:52:15.366Z

- Posts: 3

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Marco.Froelich (2026-07-17T10:52:15.425Z)

Hi,
Seems that all GPUs on Bamboo are in draining or drained state due to:
```
sinfo -R
health_cuda___nvidia root      2026-07-17T12:13:31 gpu[001-004,006-008,010-011]
health_ps___blocked  slurm     2026-07-17T11:45:20 gpu009
```
Kind regards,
froelicm


## Post 2 by @Yann.Sagon (2026-07-20T07:38:38.795Z)

Dear @Marco.Froelich[@Marco.Froelich](https://hpc-community.unige.ch/u/marco.froelich) we had again to reboot the nodes upgrade the kernel for security reason. Sorry for the inconvenience.
Best regards
Yann


## Post 3 by @Marco.Froelich (2026-07-20T07:57:02.928Z)

Hi @Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon), thanks for your response!
No worries at all, thanks for keeping the systems safe!
