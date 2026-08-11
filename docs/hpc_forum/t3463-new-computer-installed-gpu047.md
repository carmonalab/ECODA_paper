# New computer installed gpu047

- Source: https://hpc-community.unige.ch/t/3463

- Created: 2024-05-29T09:51:48.623Z

- Posts: 1

- Category: 9

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Gael.Rossignol (2024-05-29T09:51:48.690Z)

Dear users,
We have installed a new gpu server: gpu047 with 8 x A5000 25GB.
```
NODELIST    NODES PARTITION            CPUS    S:C:T MEMORY TMP_DISK GRES                                            AVAIL_FEATURES
gpu047          1 shared-gpu            128   2:64:1 512000  1500000 gpu:ampere:8,VramPerGpu:no_consume:25G EPYC-7763,V8,COMPUTE_CAPABILITY_8_6,COMPUTE_TYPE_AMPERE,S
gpu047          1 private-dpnc-gpu      128   2:64:1 512000  1500000 gpu:ampere:8,VramPerGpu:no_consume:25G EPYC-7763,V8,COMPUTE_CAPABILITY_8_6,COMPUTE_TYPE_AMPERE,S
```
Best regards,
