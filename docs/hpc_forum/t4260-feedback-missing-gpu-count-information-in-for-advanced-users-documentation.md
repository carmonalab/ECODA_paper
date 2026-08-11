# Feedback: Missing GPU count information in "for_advanced_users" documentation

- Source: https://hpc-community.unige.ch/t/4260

- Created: 2026-03-27T08:27:00.703Z

- Posts: 3

- Category: 10

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @maciej.falkiewicz (2026-03-27T08:27:00.760Z)

Hi @HPCTeam,
I hope you’re all doing well.
I’ve been consulting the “for_advanced_users” section of our cluster documentation recently and noticed that the specific breakdown of GPU counts per node seems to be missing or no longer visible.
Having this information readily available is incredibly helpful for optimizing job scripts and ensuring we’re requesting resources efficiently. Would it be possible to restore those details to the documentation?
Best regards,
Maciej


## Post 2 by @Gael.Rossignol (2026-03-27T14:18:24.458Z)

maciej.falkiewicz:
> breakdown of GPU counts per node seems to be missing o
Dear Maciej,
You are right this informations is no longer visible because when this was provided the documentation was not enough trivial (many lines for same kind of gpu).
But I think this information can be easily find using :
```
(bamboo)-[root@slurm1 ~]$ sinfo -o "%N %G"
NODELIST GRES
cpu[001-052] (null)
gpu004 gpu:nvidia_h100_nvl:1(S:2),VramPerGpu:no_consume:94G
gpu006 gpu:nvidia_h200_nvl:4(S:0,2),VramPerGpu:no_consume:141G
gpu[007,011] gpu:nvidia_rtx_pro_6000_blackwell:4(S:0-3),VramPerGpu:no_consume:96G
gpu008 gpu:nvidia_rtx_pro_6000_blackwell:4(S:3-4),VramPerGpu:no_consume:96G
gpu[009-010] gpu:nvidia_geforce_rtx_5090:4(S:0-3),VramPerGpu:no_consume:32G
gpu003 gpu:nvidia_a100_80gb_pcie:4(S:0-3),VramPerGpu:no_consume:80G
gpu[001-002] gpu:nvidia_geforce_rtx_3090:8(S:1-3,5-7),VramPerGpu:no_consume:24G
```
Best regards,


## Post 3 by @maciej.falkiewicz (2026-05-08T08:21:26.019Z)

Dear Gael,
Thank you for the useful command!
Best regards,
Maciej
