# Can't log in to bamboo

- Source: https://hpc-community.unige.ch/t/3608

- Created: 2024-08-26T08:30:11.747Z

- Tags: bamboo

- Posts: 4

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @daniel.forerosanchez (2024-08-26T08:30:11.783Z)

Hi, I am having issues now logging in into bamboo, I noticed the issue has been there since yesterday but it still persists.
```
>ssh -vv dforeros@login1.bamboo.hpc.unige.ch
OpenSSH_8.9p1 Ubuntu-3ubuntu0.10, OpenSSL 3.0.2 15 Mar 2022
debug1: Reading configuration data /home/daniel/.ssh/config
debug1: /home/daniel/.ssh/config line 1: Applying options for login1.bamboo.hpc.unige.ch
debug1: Reading configuration data /etc/ssh/ssh_config
debug1: /etc/ssh/ssh_config line 19: include /etc/ssh/ssh_config.d/*.conf matched no files
debug1: /etc/ssh/ssh_config line 21: Applying options for *
debug2: resolving "login1.bamboo.hpc.unige.ch" port 22
debug1: Connecting to login1.bamboo.hpc.unige.ch [129.194.9.186] port 22.
```
and times out.
Thanks for your help


## Post 2 by @Yoann.Boget (2024-08-26T09:02:57.826Z)

Same here.
It seems down


## Post 3 by @Maura.Brunetti (2024-08-26T09:20:36.992Z)

Hello, same problem, and I confirm it started yesterday


## Post 4 by @Adrien.Albert (2024-08-26T12:04:25.816Z)

Dear All,
It should work again
