# Current issue with FileZilla

- Source: https://hpc-community.unige.ch/t/3889

- Created: 2025-03-27T10:02:45.799Z

- Tags: baobab

- Posts: 8

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @nicolas.clairis (2025-03-27T10:02:45.861Z)

I’m experiencing a weird issue right now. The access to the Baobab server via Putty is working perfectly fine and I can see my script is currently running, but when trying to connect to sftp://login1.baobab.hpc.unige.ch via FileZilla to transfer files I’m receiving the following error message:
Statut : Access denied
Erreur : Échec de l’authentification.
Erreur : Erreur critique : Impossible d’établir une connexion au serveur
I tried multiple times since at least half an hour now and keep getting the same error message over and over. Any idea of what’s going on there?


## Post 2 by @nicolas.clairis (2025-03-27T11:01:24.692Z)

Actually, I tried to disconnect and reconnect in Putty and now my Putty access also seems denied, but I asked someone else who can connect with no issue so I’m guessing something weird happened with my account?


## Post 3 by @Lilou.Dehondt (2025-03-27T11:10:30.943Z)

Hi, sorry for not helping, I just want to let the HPC team know that both me and my colleague are currently unable to connect to Baobab too, via Putty and FileZilla. I therefore think it’s a wider issue as it does not only affect a single account


## Post 4 by @Zhongwei.Zhang (2025-03-27T11:12:51.275Z)

I am having the same issue when trying to connect to Baobab via SSH. My access has been denied since this morning.


## Post 5 by @Chloe.Mayere (2025-03-27T11:38:50.603Z)

Same here with ssh… Asks for password in loop.


## Post 6 by @Mathias.Elbaz (2025-03-27T11:51:27.323Z)

Same for me. It asks for my password in loop


## Post 8 by @Yann.Sagon (2025-03-28T09:05:37.664Z)

Dear all,
related to this?
Access denied to Baobab server for multiple people[Access denied to Baobab server for multiple people](https://hpc-community.unige.ch/t/access-denied-to-baobab-server-for-multiple-people/3893/9) HPC issues[HPC issues](https://hpc-community.unige.ch/c/hpc-support/hpc-issues/14)
> Dear Users, 
As mentioned in our previous communication, we have been busy with Bamboo maintenance. Thank you for your patience and understanding. 
This morning, I sent an email requesting everyone’s help in keeping Baobab running smoothly by removing unnecessary data from scratch storage. Your cooperation is greatly appreciated. 
Additionally, we have been working to identify users who may be incorrectly using the login node instead of a compute node, which could have contributed to instability…


## Post 9 by @nicolas.clairis (2025-03-28T09:15:52.750Z)

Yes thanks for reopening the login and fixing the issue
