# New compute node on Bamboo: gpu004.bamboo

- Source: https://hpc-community.unige.ch/t/4009

- Created: 2025-07-17T15:36:05.751Z

- Posts: 1

- Category: 9

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yann.Sagon (2025-07-17T15:36:05.830Z)

Dear users, please welcome a new family member: gpu004.bamboo this is the first compute node with an H100 card!
```
(bamboo)-[root@gpu004 ~]$ sinfo -N --format="%10N %.6D %20P %.4c %.8z %.6m %.8d %G %.57f" -n gpu004
NODELIST    NODES PARTITION            CPUS    S:C:T MEMORY TMP_DISK GRES                                            AVAIL_FEATURES
gpu004          1 shared-gpu             64   1:64:1 773000  7000000 gpu:nvidia_h100_nvl:1(S:0),VramPerGpu:no_consume:96G EPYC-9554,V11,DOUBLE_PRECISION_GPU,COMPUTE_CAPABILITY_9_0
gpu004          1 private-carmona-gpu    64   1:64:1 773000  7000000 gpu:nvidia_h100_nvl:1(S:0),VramPerGpu:no_consume:96G EPYC-9554,V11,DOUBLE_PRECISION_GPU,COMPUTE_CAPABILITY_9_0
```
