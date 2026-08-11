# New compute node installed gpu050.baobab

- Source: https://hpc-community.unige.ch/t/3994

- Created: 2025-06-27T12:09:08.611Z

- Posts: 1

- Category: 9

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Gael.Rossignol (2025-06-27T12:09:08.689Z)

Hello,
I pleased to announce we have installed a new gpu server: gpu050 :slight_smile:
```
NODELIST    NODES PARTITION            CPUS    S:C:T MEMORY TMP_DISK GRES                                            AVAIL_FEATURES
gpu050          1 shared-gpu             96   1:96:1 768000  7000000 gpu:nvidia_rtx_5000_ada_generation:4(S:0),VramPerGpu:no_consume:32G EPYC-9654,V12,COMPUTE_CAPABILITY_8_9,COMPUTE_TYPE_LOVELAC
gpu050          1 private-dpnc-gpu       96   1:96:1 768000  7000000 gpu:nvidia_rtx_5000_ada_generation:4(S:0),VramPerGpu:no_consume:32G EPYC-9654,V12,COMPUTE_CAPABILITY_8_9,COMPUTE_TYPE_LOVELAC
```
Best regards,
