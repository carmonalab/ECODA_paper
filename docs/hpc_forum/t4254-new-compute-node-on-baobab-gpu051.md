# New compute node on Baobab gpu051

- Source: https://hpc-community.unige.ch/t/4254

- Created: 2026-03-18T10:56:25.875Z

- Tags: baobab

- Posts: 1

- Category: 9

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yann.Sagon (2026-03-18T10:56:26.053Z)

Dear users, please welcome a new family member: gpu051.baobab
```
sinfo -N --format="%10N %.6D %20P %.4c %.8z %.6m %.8d %G %.57f" -n gpu051
NODELIST    NODES PARTITION            CPUS    S:C:T MEMORY TMP_DISK GRES                                            AVAIL_FEATURES
gpu051          1 private-gsem-gpu      128   8:16:1 384000  7123000 gpu:nvidia_rtx_pro_4500_blackwell:8(S:0-7),VramPerGpu:no_consume:32G EPYC-9554,V11,COMPUTE_CAPABILITY_12_0,COMPUTE_TYPE_BLACKW
gpu051          1 shared-gpu            128   8:16:1 384000  7123000 gpu:nvidia_rtx_pro_4500_blackwell:8(S:0-7),VramPerGpu:no_consume:32G EPYC-9554,V11,COMPUTE_CAPABILITY_12_0,COMPUTE_TYPE_BLACKW
```
