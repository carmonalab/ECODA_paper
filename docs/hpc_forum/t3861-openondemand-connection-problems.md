# OpenOnDemand connection problems

- Source: https://hpc-community.unige.ch/t/3861

- Created: 2025-03-13T07:00:01.353Z

- Posts: 4

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Kyrylo.Krasnykov (2025-03-13T07:00:01.390Z)

## Primary informations
Username: krasnyko
Cluster: baobab, bamboo
## Description
Dear HPC team,
I have encountered an issue with my OnDemand connection since yesterday (12.03). After going through SWITCH edu  log-in steps, the following error is thrown:
Error – can’t find user for krasnyk3
Run ‘nginx_stage --help’ to see a full list of available command line options.
Is it cluster account-related? The username krasnyk3 used to be linked to my student e-mail address (etu-unige.ch), and has already been expired months ago.
Thanks a lot in advance,
Kyrylo Krasnykov.


## Post 2 by @Adrien.Albert (2025-03-13T10:25:28.871Z)

It should work now, could you try again ?


## Post 3 by @Kyrylo.Krasnykov (2025-03-13T13:15:06.653Z)

It works! Thanks a lot!
May I know the source of the problem?
Thanks,
Kyrylo


## Post 4 by @Adrien.Albert (2025-03-13T14:50:07.804Z)

hi @Kyrylo.Krasnykov[@Kyrylo.Krasnykov](https://hpc-community.unige.ch/u/kyrylo.krasnykov)
Here more information:
[2025] Current issues on HPC Cluster[[2025] Current issues on HPC Cluster](https://hpc-community.unige.ch//hpc-community.unige.ch/t/2025-current-issues-on-hpc-cluster/3788/4) HPC ChangeLog[HPC ChangeLog](https://hpc-community.unige.ch/c/hpc-announce/hpc-changelog/9)
> Description
We encountered an issue with OpenOnDemand authentication. While attempting to fix outsider access in collaboration with the Authentication team, the authentication rules were affected. As a result, some users with dual identities (Collaborator/Student) may have experienced account mismatches with the HPC system. 
Resolution
The fix has been rolled back, and the issue should no longer be present. We plan to test an alternative solution to allow outsider access to OpenOnDemand. 
Status…
