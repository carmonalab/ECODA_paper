# Login Issues on Yggdrasil

- Source: https://hpc-community.unige.ch/t/3816

- Created: 2025-02-10T08:32:56.502Z

- Posts: 5

- Category: 8

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @daniel.forerosanchez (2025-02-10T08:32:56.543Z)

Hi, since friday I have been experiencing log-in issues on Yggdrasil and occasionally Bamboo. Today I can use Bamboo but not Yggdrasil.
```
ssh -vv dforeros@login1.yggdrasil.hpc.unige.ch
OpenSSH_9.6p1 Ubuntu-3ubuntu13.5, OpenSSL 3.0.13 30 Jan 2024
debug1: Reading configuration data /home/daniel/.ssh/config
debug1: /home/daniel/.ssh/config line 5: Applying options for login1.yggdrasil.hpc.unige.ch
debug1: Reading configuration data /etc/ssh/ssh_config
debug1: /etc/ssh/ssh_config line 19: include /etc/ssh/ssh_config.d/*.conf matched no files
debug1: /etc/ssh/ssh_config line 21: Applying options for *
debug2: resolving "login1.yggdrasil.hpc.unige.ch" port 22
debug1: Connecting to login1.yggdrasil.hpc.unige.ch [129.194.64.11] port 22.
debug1: connect to address 129.194.64.11 port 22: Connection refused
ssh: connect to host login1.yggdrasil.hpc.unige.ch port 22: Connection refused
```
Thanks in advance


## Post 2 by @Adrien.Albert (2025-02-10T16:34:03.334Z)

hi @daniel.forerosanchez[@daniel.forerosanchez](https://hpc-community.unige.ch/u/daniel.forerosanchez)
Friday, You have reached the maximum number of failed authentications and have been banned for 15 minutes.
The logs show that you were able to log in today. Is the problem still there?


## Post 3 by @daniel.forerosanchez (2025-02-20T13:52:49.146Z)

Hi, It was ok for a while there but now it is super slow on logging in from the terminal and won’t do it at all from VS code. At the moment I can’t really work because even if I manage to log in, it is super slow.


## Post 4 by @Adrien.Albert (2025-02-20T14:09:01.575Z)

Hi,
The user has been contacted and login node rebooted


## Post 5 by @daniel.forerosanchez (2025-02-24T08:57:15.750Z)

Hi again, can’t login now either. Not from the terminal nor VS Code.
