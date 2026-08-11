# New compute node on Bamboo: gpu005.bamboo

- Source: https://hpc-community.unige.ch/t/4086

- Created: 2025-09-16T15:36:41.640Z

- Tags: bamboo

- Posts: 1

- Category: 9

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yann.Sagon (2025-09-16T15:36:41.758Z)

Dear users, please welcome a new family member: gpu005.bamboo
```
sinfo -N --format="%10N %.6D %20P %.4c %.8z %.6m %.8d %G %.57f" -n gpu005
NODELIST    NODES PARTITION            CPUS    S:C:T MEMORY TMP_DISK GRES                                            AVAIL_FEATURES
gpu005          1 shared-gpu            128   2:64:1 773000  7000000 gpu:nvidia_h200_nvl:4(S:0),VramPerGpu:no_consume:141G EPYC-9554,V11,DOUBLE_PRECISION_GPU,COMPUTE_CAPABILITY_9_0
gpu005          1 private-burgi-gpu     128   2:64:1 773000  7000000 gpu:nvidia_h200_nvl:4(S:0),VramPerGpu:no_consume:141G EPYC-9554,V11,DOUBLE_PRECISION_GPU,COMPUTE_CAPABILITY_9_0
```
