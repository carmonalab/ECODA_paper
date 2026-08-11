# Baobab firewall

- Source: https://hpc-community.unige.ch/t/3978

- Created: 2025-06-10T12:56:37.349Z

- Tags: baobab

- Posts: 3

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Paul.Coppin (2025-06-10T12:56:37.401Z)

Dear HPC,
Since the upgrade last week, we can no longer access the login node from our machines at Dufour which are in the IP range: 192.33.218.XXX
This is independent of the login node restarting earlier this afternoon, i.e. issue already observed yesterday and not resolved by the reboot just now.
As an example, from gridvm10.unige.ch (192.33.218.240):
```
[coppinp@gridvm10 ~]$ ping login1.baobab.hpc.unige.ch
PING login1.baobab.hpc.unige.ch (129.194.9.190) 56(84) bytes of data.
--- login1.baobab.hpc.unige.ch ping statistics ---
4 packets transmitted, 0 received, 100% packet loss, time 2999ms

[coppinp@gridvm10 ~]$ ssh -v coppinp@login1.baobab.hpc.unige.ch
OpenSSH_7.4p1, OpenSSL 1.0.2k-fips 26 Jan 2017
debug1: Reading configuration data /etc/ssh/ssh_config
debug1: /etc/ssh/ssh_config line 58: Applying options for *
debug1: Connecting to login1.baobab.hpc.unige.ch [129.194.9.190] port 22. 
# then the connection hangs
```
Could it be possible that we are blocked by the firewall?
All the best,
Paul


## Post 2 by @Gael.Rossignol (2025-06-11T08:45:37.368Z)

Dear Paul,
You’re right some routes was not correctly deployed. Config has been updated could you please check if it’s working fine?
Best regards,


## Post 3 by @Paul.Coppin (2025-06-11T13:55:37.971Z)

Hi Gael,
Indeed it works fine now. Thank you for the fix!
All the best,
Paul
