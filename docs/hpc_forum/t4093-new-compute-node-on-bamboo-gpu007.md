# New compute node on Bamboo gpu007

- Source: https://hpc-community.unige.ch/t/4093

- Created: 2025-09-18T13:45:33.432Z

- Tags: bamboo

- Posts: 1

- Category: 9

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yann.Sagon (2025-09-18T13:45:33.548Z)

Dear users, please welcome a new family member: gpu007.bamboo
```
(bamboo)-[root@login1 ~]$ sinfo -N --format="%10N %.6D %20P %.4c %.8z %.6m %.8d %G %.57f" -n gpu007
NODELIST    NODES PARTITION            CPUS    S:C:T MEMORY TMP_DISK GRES                                            AVAIL_FEATURES
gpu007          1 shared-gpu            128  1:128:1 773000  7000000 gpu:nvidia_rtx_pro_6000_blackwell:4(S:0),VramPerGpu:no_consume:97G EPYC-9754,V11,DOUBLE_PRECISION_GPU,COMPUTE_CAPABILITY_12_
gpu007          1 private-engelke-gpu   128  1:128:1 773000  7000000 gpu:nvidia_rtx_pro_6000_blackwell:4(S:0),VramPerGpu:no_consume:97G EPYC-9754,V11,DOUBLE_PRECISION_GPU,COMPUTE_CAPABILITY_12_
```
