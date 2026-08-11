# New compute node on Bamboo gpu[009-010]

- Source: https://hpc-community.unige.ch/t/4256

- Created: 2026-03-18T14:15:51.027Z

- Tags: bamboo

- Posts: 1

- Category: 9

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yann.Sagon (2026-03-18T14:15:51.228Z)

Dear users, please welcome a new family member: gpu[009-010]
```
sinfo -N --format="%10N %.6D %20P %.4c %.8z %.6m %.8d %G %.57f" -n gpu[009-010]
NODELIST    NODES PARTITION            CPUS    S:C:T MEMORY TMP_DISK GRES                                            AVAIL_FEATURES
gpu009          1 private-kalousis-gpu   96   4:24:1 384000  7123000 gpu:nvidia_geforce_rtx_5090:4(S:0-3),VramPerGpu:no_consume:32G EPYC-9654,V12,COMPUTE_CAPABILITY_12_0,COMPUTE_TYPE_BLACKW
gpu009          1 shared-gpu             96   4:24:1 384000  7123000 gpu:nvidia_geforce_rtx_5090:4(S:0-3),VramPerGpu:no_consume:32G EPYC-9654,V12,COMPUTE_CAPABILITY_12_0,COMPUTE_TYPE_BLACKW
gpu010          1 private-kalousis-gpu   96   4:24:1 384000  7123000 gpu:nvidia_geforce_rtx_5090:4(S:0-3),VramPerGpu:no_consume:32G EPYC-9654,V12,COMPUTE_CAPABILITY_12_0,COMPUTE_TYPE_BLACKW
gpu010          1 shared-gpu             96   4:24:1 384000  7123000 gpu:nvidia_geforce_rtx_5090:4(S:0-3),VramPerGpu:no_consume:32G EPYC-9654,V12,COMPUTE_CAPABILITY_12_0,COMPUTE_TYPE_BLACKW
```
