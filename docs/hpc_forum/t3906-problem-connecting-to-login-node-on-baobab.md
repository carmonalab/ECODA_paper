# Problem connecting to login node on Baobab

- Source: https://hpc-community.unige.ch/t/3906

- Created: 2025-04-03T10:03:45.953Z

- Tags: baobab

- Posts: 9

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Lauric.Ferrat (2025-04-03T10:03:45.953Z)

Hi Yann,
I don’t know if my problem is related to this too;
when trying to login to Baobab using ssh I am 90% of the time refused with either this message:
Connection closed by 129.194.9.190 port 22
or this one
ssh: connect to host login1.baobab.hpc.unige.ch port 22: Connection refused
the other 10% of the time I can login.


## Post 2 by @Yann.Sagon (2025-04-03T10:10:45.527Z)

Dear @Lauric.Ferrat[@Lauric.Ferrat](https://hpc-community.unige.ch/u/lauric.ferrat)
I tried to log using your username.
I saw this issue:
```
(baobab)-[root@admin1 ~]$ su - ferratl

Due to MODULEPATH changes, the following have been reloaded:
  1) Tcl/8.6.13     2) binutils/2.40     3) bzip2/1.0.8     4) libffi/3.4.4     5) libreadline/8.2     6) ncurses/6.4     7) zlib/1.2.13

The following have been reloaded with a version change:
  1) GCCcore/12.3.0 => GCCcore/13.2.0     2) Python/3.11.3 => Python/3.11.5     3) SQLite/3.42.0 => SQLite/3.43.1     4) XZ/5.4.2 => XZ/5.4.4
```
You should check your `.bashr` file and correct it to load only one toolchain or compiler version.
Can you fix that first and try again to see if that solves the issue?
Best
Yann


## Post 3 by @Lauric.Ferrat (2025-04-03T17:57:55.350Z)

thank you for the suggestion!
Today I have not been able to login but I will try again tomorrow.
Kind regards,
Lauric


## Post 4 by @Lauric.Ferrat (2025-04-04T08:20:02.958Z)

Yann.Sagon:
> toolchain
Dear Yann,
I have managed to login using the ondemand and with the terminal of Rstudio I have deleted .bashrc but I still cannot log in with my terminal.
I had installed a public key on baobab and so I tried to force the use of a password
‘ssh -o PreferredAuthentications=password -o PubkeyAuthentication=no ferratl@@login1.baobab.hpc.unige.ch’ but I encounter this error:
Permission denied (publickey,password,keyboard-interactive,hostbased)


## Post 5 by @Yann.Sagon (2025-04-04T08:28:15.427Z)

Lauric.Ferrat:
> ferratl@@login1.baobab.hpc.unige.ch’
There is a typo: “double @”.


## Post 6 by @Lauric.Ferrat (2025-04-04T08:53:35.352Z)

That was it… :face_exhaling:
Thank you


## Post 7 by @Lauric.Ferrat (2025-04-04T08:59:08.272Z)

In fact the problem just came back it not just a typos, I was able to connect once and then when I retried I had “connect to host login1.baobab.hpc.unige.ch port 22: Connection refused”


## Post 8 by @Yann.Sagon (2025-04-04T09:46:29.399Z)

What is your ip please? Are you trying to connect as well using filezilla for example using a wrong password? If this the case, your ip is blacklisted for 15 minutes.
See here for details: hpc:faq [eResearch Doc][hpc:faq [eResearch Doc]](https://doc.eresearch.unige.ch/hpc/faq#when_i_tried_to_connect_to_the)


## Post 9 by @Lauric.Ferrat (2025-04-04T10:15:28.070Z)

2a01:cb15:6f:da00:b544:6925:f2ed:bac0%
I am not using a filezilla or equivalent - things has been working well for the past 30 minutes - so maybe as I was fiddling with different way to connect I just got on the ip blacklist.
Thank you again for your help!
