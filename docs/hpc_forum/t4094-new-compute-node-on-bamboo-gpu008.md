# New compute node on Bamboo gpu008

- Source: https://hpc-community.unige.ch/t/4094

- Created: 2025-09-18T15:24:14.650Z

- Tags: bamboo

- Posts: 1

- Category: 9

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yann.Sagon (2025-09-18T15:24:14.762Z)

Dear users, please welcome a new family member: gpu008.bamboo
```
(bamboo)-[root@admin1 ~]$ sinfo -N --format="%10N %.6D %20P %.4c %.8z %.6m %.8d %G %.57f" -n gpu008
NODELIST    NODES PARTITION            CPUS    S:C:T MEMORY TMP_DISK GRES                                            AVAIL_FEATURES
gpu008          1 shared-gpu            128   2:64:1 773000  7000000 gpu:nvidia_rtx_pro_6000_blackwell:4(S:0-1),VramPerGpu:no_consume:97G EPYC-9554,V11,DOUBLE_PRECISION_GPU,COMPUTE_CAPABILITY_12_
gpu008          1 private-hepia-gpu     128   2:64:1 773000  7000000 gpu:nvidia_rtx_pro_6000_blackwell:4(S:0-1),VramPerGpu:no_consume:97G EPYC-9554,V11,DOUBLE_PRECISION_GPU,COMPUTE_CAPABILITY_12_
```
