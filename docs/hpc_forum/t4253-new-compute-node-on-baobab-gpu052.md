# New compute node on Baobab gpu052

- Source: https://hpc-community.unige.ch/t/4253

- Created: 2026-03-18T09:17:24.288Z

- Tags: baobab

- Posts: 1

- Category: 9

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yann.Sagon (2026-03-18T09:17:24.481Z)

Dear users, please welcome a new family member: gpu052.baobab
```
sinfo -N --format="%10N %.6D %20P %.4c %.8z %.6m %.8d %G %.57f" -n gpu052
NODELIST    NODES PARTITION            CPUS    S:C:T MEMORY TMP_DISK GRES                                            AVAIL_FEATURES
gpu052          1 private-kruse-gpu      64   4:16:1 768000  7123000 gpu:nvidia_h100_nvl:2(S:1,3),VramPerGpu:no_consume:94G EPYC-9554,V11,COMPUTE_CAPABILITY_9_0,COMPUTE_TYPE_HOPPER,
gpu052          1 shared-gpu             64   4:16:1 768000  7123000 gpu:nvidia_h100_nvl:2(S:1,3),VramPerGpu:no_consume:94G EPYC-9554,V11,COMPUTE_CAPABILITY_9_0,COMPUTE_TYPE_HOPPER,
```
