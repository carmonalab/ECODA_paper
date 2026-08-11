# Dear users, please welcome a new family member: gpu006.bamboo

- Source: https://hpc-community.unige.ch/t/4087

- Created: 2025-09-16T15:39:11.629Z

- Tags: bamboo

- Posts: 1

- Category: 9

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yann.Sagon (2025-09-16T15:39:11.736Z)

Dear users, please welcome a new family member: gpu006.bamboo
```
sinfo -N --format="%10N %.6D %20P %.4c %.8z %.6m %.8d %G %.57f" -n gpu006
NODELIST    NODES PARTITION            CPUS    S:C:T MEMORY TMP_DISK GRES                                            AVAIL_FEATURES
gpu006          1 shared-gpu             96   1:96:1 773000  7000000 gpu:nvidia_h200_nvl:4(S:0),VramPerGpu:no_consume:141G EPYC-9654,V11,DOUBLE_PRECISION_GPU,COMPUTE_CAPABILITY_9_0
gpu006          1 private-drim-gpu       96   1:96:1 773000  7000000 gpu:nvidia_h200_nvl:4(S:0),VramPerGpu:no_consume:141G EPYC-9654,V11,DOUBLE_PRECISION_GPU,COMPUTE_CAPABILITY_9_0
```
