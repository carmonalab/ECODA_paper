# New gpu node installed gpu049

- Source: https://hpc-community.unige.ch/t/3779

- Created: 2025-01-08T15:39:37.164Z

- Posts: 1

- Category: 9

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Gael.Rossignol (2025-01-08T15:39:37.249Z)

We have installed a new gpu server: gpu049 with 8 x RTX4090 24GB. This is a new architecture with latest gpu hardware.
```
NODELIST    NODES PARTITION            CPUS    S:C:T MEMORY TMP_DISK GRES                                            AVAIL_FEATURES
gpu049          1 shared-gpu            128   2:64:1 384000  7000000 gpu:nvidia_geforce_rtx_4090:8(S:0-1),VramPerGpu:no_consume:24G EPYC-9554,V11,COMPUTE_CAPABILITY_8_9,COMPUTE_TYPE_AMPERE,
gpu049          1 private-dpnc-gpu      128   2:64:1 384000  7000000 gpu:nvidia_geforce_rtx_4090:8(S:0-1),VramPerGpu:no_consume:24G EPYC-9554,V11,COMPUTE_CAPABILITY_8_9,COMPUTE_TYPE_AMPERE,
```
Best regards,
