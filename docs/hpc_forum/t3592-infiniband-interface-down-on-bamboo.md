# Infiniband interface down on Bamboo

- Source: https://hpc-community.unige.ch/t/3592

- Created: 2024-08-12T11:54:34.574Z

- Tags: bamboo

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Berk.Gercek (2024-08-12T11:54:34.613Z)

Hello all,
I’m trying to run a dask scheduler/worker system on bamboo. I usually use infiniband (`ib0`) on both the scheduler node and the worker nodes. However in this case I was getting timeout when the workers attempted to contact the scheduler.
Switching to ethernet (`eno1`) was enough to resolve the problem.
Is there a current infiniband issue on Bamboo?
Regards,
Berk


## Post 2 by @Gael.Rossignol (2024-08-13T08:47:40.784Z)

Dear Berk,
Actually IB is working fine on bamboo we have done a lot of bench.
As beegfs is using rdma, we can assume that all seems to be ok.
Could you please share your config files to check?
Best regards,
