# X11 connection rejected

- Source: https://hpc-community.unige.ch/t/3286

- Created: 2024-02-01T13:42:03.497Z

- Posts: 3

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Matthias.Kruckow (2024-02-01T13:42:03.497Z)

In the last days I more and more often have seen the X11 rejection, too.
I have seen that the configuration about this should be in `/etc/ssh/sshd_config`. Interestingly, it turned out, that I don’t have read permission for this on yggdrasil and this file was last modified a few days ago on Jan 26 (potentially coinciding the start of showing the rejection messages).
Hence, I’m wondering, whether there is something wrong in this configuration file to cause this rejection messages to pop up.


## Post 2 by @Yann.Sagon (2024-02-02T08:21:33.719Z)

Hi,
I’ve created a new post thread as the one you were answering is already 3 years old.
The files `/etc/ssh/sshd_config` are identical on both clusters. This file was modified the 26th of January due to then end of the maintenance: we disable ssh access during the maintenance. We didn’t changed something else on the file.
Please give us more details:
- from which OS are you connecting?
- which ssh client?
- what is the command line you use to connect to Yggdrasil?
- does the error pop immediately when connecting or later when you launch something graphical?
Is this an option for you to use x2go[x2go](https://doc.eresearch.unige.ch/hpc/access_the_hpc_clusters#gui_accessdesktop_with_x2go)?


## Post 3 by @Matthias.Kruckow (2024-02-02T08:43:23.532Z)

Here some more details you asked for:
- My OS is: Ubuntu 22.04.3 LTS
- my ssh version is: OpenSSH_8.9p1 Ubuntu-3ubuntu0.6, OpenSSL 3.0.2
- I usually connect via: `ssh -Y kruckow@login1.yggdrasil.hpc.unige.ch` (using the `-X` option instead of `-Y` doesn’t change it)
- The error pops up form time to time later. It is not connected to any commands I’m running. It looks to happen sometimes, when an automatic handshake is done in the background. The error message doesn’t influence look to not influence any operations. Even it pops up, I can open a new window or an open window would stay intact (here e.g. using `display`).
