# Cannot connect to Baobab

- Source: https://hpc-community.unige.ch/t/4051

- Created: 2025-08-18T09:55:53.005Z

- Tags: baobab

- Posts: 6

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Heloise.Allaman (2025-08-18T09:55:53.050Z)

Hello,
I’ve been trying to connect to Baobab via SSH fo two hours, but I get error " port 22: Connection refused". Maybe some of you have any idea of the problem ?
Thank you


## Post 2 by @Yann.Sagon (2025-08-18T12:58:17.686Z)

According to the logs, you could connect today afternoon. Is this correct or do you still have the issue?
```
Aug 18 13:58:42 login1 sshd[3212246]: Accepted keyboard-interactive/pam for allamanh from 10.20.167.xx port 51135 ssh2
Aug 18 13:58:42 login1 sshd[3212246]: pam_unix(sshd:session): session opened for user allamanh(uid=xxx) by allamanh(uid=0)
```
Best
Yann


## Post 3 by @Heloise.Allaman (2025-08-18T13:13:56.670Z)

Hello,
I could ideed conect this afternoon. However, this error happens regularly, i.e. i cannot connect for some minutes / some hours, especially when I try to connect via VsCode. Do you have any idea of the issue ?
Best,
Héloïse


## Post 4 by @Yann.Sagon (2025-08-18T13:26:59.676Z)

Heloise.Allaman:
> this error happens regularly, i.e. i cannot connect for some minutes / some hours, especially when I try to connect via VsCode. Do you have any idea of the issue ?
You probably have a bad Baobab password stored somewhere, for example in VScode. If you do too many connection attempts with a wrong password, your ip address is blacklisted for 15 minutes.


## Post 5 by @Heloise.Allaman (2025-08-19T07:09:19.527Z)

I don’t think that is the issue because I always type my password manually, and I sometimes get the error “Connection failed“ without even typing my password. Furthermore, altough a bit less frequent than with VScode, I get a similar error when trying to connect directly via the terminal using SSH


## Post 6 by @Yann.Sagon (2025-08-19T12:49:57.275Z)

I’ve checked the access logs from today with your username, there was 13 failed attempts.
