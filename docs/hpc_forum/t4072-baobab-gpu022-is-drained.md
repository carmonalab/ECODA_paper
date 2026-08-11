# [baobab] gpu022 is drained

- Source: https://hpc-community.unige.ch/t/4072

- Created: 2025-08-29T13:08:40.064Z

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Ramon.CalvoGonzalez (2025-08-29T13:08:40.119Z)

## Primary informations
Username: calvogon
Cluster: Baobab
## Description
After a failed CUDA job, the GPUs became unusable and I ended the job. The node is now in `drain` state, with the reason: `health_cuda___nvidia-smi`. Could it be that it needs to have manual intervention? Running `nvidia-smi –gpu-reset`might help, although I’m not sure (I can’t run it myself since it needs su).
Is there a way I can reset the GPUs myself if this ever happens again? Just so I don’t have to bother the HPC team.
Best, Ramon.


## Post 2 by @Gael.Rossignol (2025-09-01T11:08:18.323Z)

Dear Ramon,
Gpu has been rebooted and set in production again, thanks for notify us about this trouble.
“gpu-burn” software didn’t detect any issue.
Best regards,
