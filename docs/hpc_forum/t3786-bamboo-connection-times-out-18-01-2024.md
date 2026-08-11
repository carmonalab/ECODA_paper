# Bamboo connection times out - 18.01.2024

- Source: https://hpc-community.unige.ch/t/3786

- Created: 2025-01-18T18:52:56.687Z

- Posts: 5

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Cody.Cardenas (2025-01-18T18:52:56.729Z)

Unable to ssh or scp Bamboo:
Expected to be able to log in or download files from bamboo, but instead I get the `port 22: connection timed out` error.
```
ssh cardenac@login1.bamboo.hpc.unige.ch
ssh: connect to host login1.bamboo.hpc.unige.ch port 22: Connection timed out
```


## Post 2 by @Berk.Gercek (2025-01-19T09:23:07.315Z)

I’ve also had similar issues for the past few months, consistently on weekends. It seems like bamboo in particular is unreachable on Sundays and sometimes Saturdays?


## Post 3 by @Yann.Sagon (2025-01-20T08:15:01.291Z)

The reason is this:
shopping
shopping173×262 3.13 KB
[shopping173×262 3.13 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/78f283df321b8e0bff0c8d09ace2a26ac15f2f4a.webp)
Joking of course. If you have time and date and cluster when you are unable to connect, we can check the logs.
edit: we can reproduce the issue and we can’t connect to the login node. We are investigating.


## Post 4 by @Yann.Sagon (2025-01-20T09:07:42.989Z)

Dear @Cody.Cardenas[@Cody.Cardenas](https://hpc-community.unige.ch/u/cody.cardenas) and @Berk.Gercek[@Berk.Gercek](https://hpc-community.unige.ch/u/berk.gercek)  the issue is solved for now, we’ll monitor the issue tightly to see what is going on.
[2025] Current issues on HPC Cluster[[2025] Current issues on HPC Cluster](https://hpc-community.unige.ch/t/2025-current-issues-on-hpc-cluster/3788) HPC ChangeLog[HPC ChangeLog](https://hpc-community.unige.ch/c/hpc-announce/hpc-changelog/9)
> Dear users, 
We have some problems with login1 on Bamboo, we are investigation to find a solution. 
For the time being, running jobs are continuing, but it isn’t possible to login to the login node. 
We’ll keep you informed, thank you for your understanding. 
Kind regards 
edit: problem solved, we’ll add the login node to our monitoring server. 
Status : Solved green_circle
start: 2025-01-18T10:00:00Z (UTC) 
end: 2025-01-20T08:00:00Z (UTC)
Best


## Post 5 by @Cody.Cardenas (2025-01-20T09:40:52.120Z)

@Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon) thanks so much for quickly solving the issue. Hopefully you can track down whats causing it.
