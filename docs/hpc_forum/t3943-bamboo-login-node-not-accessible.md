# Bamboo login node not accessible

- Source: https://hpc-community.unige.ch/t/3943

- Created: 2025-04-30T14:08:42.664Z

- Tags: bamboo

- Posts: 6

- Category: 1

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Lucille.Delisle1 (2025-04-30T14:08:42.710Z)

Hi,
I don’t manage to access bamboo login node by ssh with or without VPN but I can access the open on demand of bamboo and I can ssh from baobab front-end…
Best,
Lucille


## Post 2 by @Adrien.Albert (2025-05-01T22:38:42.211Z)

on Bamboo Logs:
```
Apr 30 16:50:07 login1 sshd[1898953]: Accepted publickey for delislel from XXX.XXX.XXX.XXX  port 59632 ssh2: RSA SHA256:En5lSqs+dQtjdP5vctO+nHekEK3E2c6WLjhhszwsj
Apr 30 16:50:07 login1 sshd[1898953]: pam_unix(sshd:session): session opened for user delislel(uid=XXXXX) by delislel(uid=0)
```
It seems that you have successfully connected to login1.bamboo


## Post 3 by @Matija.Trickovic (2025-05-02T07:47:46.083Z)

Hello, I have the same issue:
ssh: connect to host login1.bamboo.hpc.unige.ch port 22: Operation timed out
Yggdrasil and Baobab work well.
Thanks,
Matija


## Post 4 by @Lucille.Delisle1 (2025-05-02T08:18:12.267Z)

I am connected to bamboo through baobab… :wink: but I would prefer to connect to bamboo directly…


## Post 5 by @Yann.Sagon (2025-05-02T08:54:32.360Z)

Dear @Lucille.Delisle1[@Lucille.Delisle1](https://hpc-community.unige.ch/u/lucille.delisle1) and @Matija.Trickovic[@Matija.Trickovic](https://hpc-community.unige.ch/u/matija.trickovic) we are trying to fix the issue, sorry if you get disconnected in the meantime.


## Post 6 by @Yann.Sagon (2025-05-02T09:09:05.628Z)

This is solved, there was a network configuration issue.
Best regards
