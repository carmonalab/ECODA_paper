# New compute node on Bamboo gpu011

- Source: https://hpc-community.unige.ch/t/4255

- Created: 2026-03-18T13:41:37.409Z

- Posts: 1

- Category: 9

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yann.Sagon (2026-03-18T13:41:37.569Z)

Dear users, please welcome a new family member: .
```
sinfo -N --format="%10N %.6D %20P %.4c %.8z %.6m %.8d %G %.57f" -n gpu011
NODELIST    NODES PARTITION            CPUS    S:C:T MEMORY TMP_DISK GRES                                            AVAIL_FEATURES
gpu011          1 shared-gpu             64   4:16:1 773000  3546000 gpu:nvidia_rtx_pro_6000_blackwell:4(S:0-3),VramPerGpu:no_consume:96G EPYC-9554,V11,COMPUTE_CAPABILITY_9_0,COMPUTE_TYPE_BLACKWE
```
