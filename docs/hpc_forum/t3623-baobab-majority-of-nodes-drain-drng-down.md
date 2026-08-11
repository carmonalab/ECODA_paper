# [baobab] Majority of nodes drain/drng/down

- Source: https://hpc-community.unige.ch/t/3623

- Created: 2024-09-02T10:02:48.894Z

- Tags: baobab

- Posts: 3

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @maciej.falkiewicz (2024-09-02T10:02:48.934Z)

## Primary informations
Username: falkiewi
Cluster: Baobab
Dear @support[@support](https://hpc-community.unige.ch/groups/support) ,
The majority of nodes on the cluster are in drain, drng or down state.
```
$ sinfo -R
REASON               USER      TIMESTAMP           NODELIST
issue-6277 : Not res slurm     2024-08-26T14:28:56 cpu[002,005,007-008,045-056,058-059,061-062]
issue-6277           root      2024-08-22T10:44:06 cpu004
health_BEEGFS__tcp_c root      2024-09-02T11:40:07 cpu[156,173-183,187-188,190-193,196-198,208,220-221,240,242,270-271,274,287,290,292-293,295,303,308-311]
Not responding       slurm     2024-08-30T21:25:06 cpu285
health_BEEGFS__tcp_c root      2024-09-02T10:52:06 cpu[247,253]
health_BEEGFS__tcp_c root      2024-09-02T10:52:07 cpu251
health_BEEGFS__tcp_c root      2024-09-02T10:52:08 cpu252
health_BEEGFS__tcp_c root      2024-09-02T11:43:07 cpu326,gpu017
health_BEEGFS__tcp_c root      2024-09-02T10:52:05 cpu[065-066,194,223,262,273,275-276,298,301,304,317]
health_BEEGFS__tcp_c root      2024-09-02T10:58:05 cpu[157-160,162,164-169,171-172,184,189,195,200-202,216-217,225,263-264,294,306-307,332],gpu005
health_BEEGFS__tcp_c root      2024-09-02T10:58:05 cpu161
health_BEEGFS__tcp_c root      2024-09-02T10:58:06 cpu185
health_BEEGFS__tcp_c root      2024-09-02T11:40:08 cpu[199,205-206,222]
health_BEEGFS__tcp_c root      2024-09-02T10:52:05 cpu[207,249-250,254-261,299-300,305,315,319,321]
health_BEEGFS__tcp_c root      2024-09-02T10:52:08 cpu224
issue-6277 : Not res slurm     2024-08-26T14:27:16 cpu186
health_cuda___nvidia root      2024-09-02T10:55:08 gpu034
health_cuda___nvidia root      2024-09-02T00:21:47 gpu002
health_cuda___nvidia root      2024-08-30T00:46:53 gpu008
Kill task failed     root      2024-08-28T03:30:07 gpu012
Kill task failed     root      2024-08-28T21:17:07 gpu018
Kill task failed     root      2024-08-31T03:39:02 gpu025
issue-6299           root      2024-08-29T15:14:10 gpu047
```
Kind regards,
Maciej Falkiewicz


## Post 2 by @Gael.Rossignol (2024-09-02T15:28:44.981Z)

Dear Maciej,
Thank you very much to notice that. I see the problem and all nodes are back in production state. We plan to upgrade storage fabric this month in order to solve those kind of issues.
Best regards,


## Post 3 by @Raphael.Rubino (2024-09-08T15:55:41.592Z)

Hello,
Lots of GPUs on baobab again in drain.
Best regards
