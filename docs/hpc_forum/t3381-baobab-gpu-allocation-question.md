# Baobab GPU Allocation Question

- Source: https://hpc-community.unige.ch/t/3381

- Created: 2024-03-20T11:56:38.297Z

- Posts: 3

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Raphael.Rubino (2024-03-20T11:56:38.360Z)

Hello,
I have a question about GPU allocation on Baobab.
When I look at the overall GPU usage with the `sinfo` command, I get the allocated nodes, for instance `gpu002` is allocated.
However, when I check the node with the `scontrol` command, I get the following:
> CfgTRES=cpu=10,mem=257000M,billing=10,gres/gpu=6,gres/gpu:titan=6
> AllocTRES=cpu=10,mem=146G,gres/gpu=4,gres/gpu:titan=4
Which indicates that 2 GPUs are still available, but all CPUs are being used, is that right or am I missing something?
Thanks for your help!


## Post 2 by @Yann.Sagon (2024-03-21T09:58:33.145Z)

Dear @Raphael.Rubino[@Raphael.Rubino](https://hpc-community.unige.ch/u/raphael.rubino)
indeed, this seems to be the case. The issue with this node is that it is old and has only 12 CPU cores and two of them are reserved Slurm Workload Manager - Core Specialization[Slurm Workload Manager - Core Specialization](https://slurm.schedmd.com/core_spec.html). Using `sinfo` it appears allocated as for Slurm all the CPUs are allocated, thus no more jobs can use this node.


## Post 3 by @Raphael.Rubino (2024-03-21T10:08:53.748Z)

Thank you for the info @Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon) !
