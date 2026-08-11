# Change associated PI

- Source: https://hpc-community.unige.ch/t/4113

- Created: 2025-10-06T07:32:49.579Z

- Posts: 2

- Category: 8

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Denis.Mongin (2025-10-06T07:32:49.636Z)

Dear HPC team
I am participating on several project with different PI.
Given that The use of baobab can lead to paying use, it is of importance to associate my activity to the proper PI. So two questions:
- How do I become a PI ? (I have projects myself)
- How do I switch PI ?
- How do I simply check that I indeed have changed PI ?
Thanks !
Denis


## Post 2 by @Yann.Sagon (2025-10-10T13:53:01.438Z)

Dear @Denis.Mongin[@Denis.Mongin](https://hpc-community.unige.ch/u/denis.mongin)
Denis.Mongin:
> How do I become a PI ? (I have projects myself)
What we expect from a PI is to be for example a team leader or someone who has budget for doing projects. At least a PI should be someone that will stay at UNIGE more time than a student. As you see, the definition is quite vague:) As you have projects yourself, you are a good candidate to be a PI.
ref: hpc:access_the_hpc_clusters [eResearch Doc][hpc:access_the_hpc_clusters [eResearch Doc]](https://doc.eresearch.unige.ch/hpc/access_the_hpc_clusters#account)
Denis.Mongin:
> How do I switch PI ?
To switch between PI: hpc:faq [eResearch Doc][hpc:faq [eResearch Doc]](https://doc.eresearch.unige.ch/hpc/faq#i_m_a_user_and_i_ve_noticed_th)
Denis.Mongin:
> How do I simply check that I indeed have changed PI ?
You have to explicitly choose the PI you want to account to when you submit a job (srun, salloc, sbatch). If you don’t do that, your default account (pi) is accounted. You can check who is your default PI like so:
```
sacctmgr show user mongin
```
Contact us by email if you want to be your own PI and or that we update your user account with a second PI.
