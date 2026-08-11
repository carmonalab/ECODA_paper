# [baobab] Issue with gpu027 Resource Allocation

- Source: https://hpc-community.unige.ch/t/4261

- Created: 2026-03-27T08:30:41.242Z

- Posts: 4

- Category: 1

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @maciej.falkiewicz (2026-03-27T08:30:41.298Z)

Hi @HPCTeam,
I’d like to report a recurring issue I’ve encountered with node `gpu027` regarding GPU memory allocation.
It appears that `gpu027` has a non-homogeneous GPU configuration where some cards have less than 80GB of VRAM. Even when specifically using the constraint:
`--gres="gpu:2,VramPerGpu:80G"`
The Slurm scheduler is still assigning jobs to this node. Unfortunately, because some of the allocated GPUs don’t meet the 80GB requirement, the jobs are hitting Out-of-Memory (OOM) errors shortly after starting.
Could you please look into this and see if there is a way to ensure jobs requesting 80GB of VRAM are not assigned to the lower-capacity cards on that node?
Thank you for your time and assistance!
Best regards,
Maciej


## Post 2 by @Gael.Rossignol (2026-03-27T14:08:36.306Z)

Dear Maciej,
You are right, the config is not currently configured to manage correctly both type of gpu.
But if you need to 80G gpu you can directly use the option : nvidia_a100_80gb_pcie.
We will try to update config later by using the gres.
Best regards,


## Post 3 by @maciej.falkiewicz (2026-03-28T20:30:59.599Z)

Dear Gael,
Thank you for the quick reply.
Gael.Rossignol:
> you can directly use the option : nvidia_a100_80gb_pcie.
There are multiple GPU models with at least 80GB VRAM in the cluster; the A100 would be overly restrictive. In the world of limited resources, we always have to make the least restrictive configurations possible!
Best regards,
Maciej


## Post 4 by @Adrien.Albert (2026-03-30T15:36:57.942Z)

Hello @maciej.falkiewicz[@maciej.falkiewicz](https://hpc-community.unige.ch/u/maciej.falkiewicz),
maciej.falkiewicz:
> In the world of limited resources, we always have to make the least restrictive configurations possible!
I agree that we should aim for the least restrictive configuration possible, and in a fully homogeneous environment, this would be much easier to achieve. :wink:
In this particular case, the situation is somewhat unusual. This node received a different GPU model for historical reasons, the card was added, later, at the request of its owner, and this creates a technical limitation rather than a configuration choice. All other GPU nodes use the same model, so this issue does not appear elsewhere in the cluster.
At the moment, Slurm does not provide any mechanism to prioritize or select specific GPUs within a node, which is what leads to this constraint.
Best Regards,
