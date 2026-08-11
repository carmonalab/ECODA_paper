# New compute installed gpu048.baobab

- Source: https://hpc-community.unige.ch/t/3640

- Created: 2024-09-19T09:16:33.936Z

- Tags: baobab

- Posts: 1

- Category: 9

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yann.Sagon (2024-09-19T09:16:34.050Z)

Dear users,
We have installed a new gpu server: gpu048 with 8 x A6000 48GB NVLINK.
```
(baobab)-[root@gpu048 ~]$ sinfo -N --format="%10N %.6D %20P %.4c %.8z %.6m %.8d %G %.76f" -n gpu048
NODELIST    NODES PARTITION            CPUS    S:C:T MEMORY TMP_DISK GRES                                                               AVAIL_FEATURES
gpu048          1 shared-gpu            128   2:64:1 512000   300000 gpu:ampere:8,VramPerGpu:no_consume:48G EPYC-7763,V10,COMPUTE_CAPABILITY_8_6,COMPUTE_TYPE_AMPER,SIMPLE_PRECISION_GPU
gpu048          1 private-kalousis-gpu  128   2:64:1 512000   300000 gpu:ampere:8,VramPerGpu:no_consume:48G EPYC-7763,V10,COMPUTE_CAPABILITY_8_6,COMPUTE_TYPE_AMPER,SIMPLE_PRECISION_GPU
```
Enjoy!
