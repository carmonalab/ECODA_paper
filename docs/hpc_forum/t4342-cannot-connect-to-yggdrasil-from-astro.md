# Cannot connect to yggdrasil from astro

- Source: https://hpc-community.unige.ch/t/4342

- Created: 2026-07-08T10:43:31.711Z

- Tags: yggdrasil

- Posts: 5

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @elia.cenci (2026-07-08T10:43:31.818Z)

Hello,
I can successfully connect to bonsai. However,  when I am trying to connect to Yggdrasil/Baobab after the SSH-key change, it fails with:
```
connect to address 129.194.9.68 port 22: Connection refused
```
I should not be an authentication issue since I do not get a permission denied report and it also fails when using:
```
ssh -F /dev/null -vvv cenci@login1.yggdrasil.hpc.unige.ch
```
Any clues why this happens?
Thank you.
Elia


## Post 2 by @Yann.Sagon (2026-07-08T14:21:53.888Z)

Dear @elia.cenci[@elia.cenci](https://hpc-community.unige.ch/u/elia.cenci) is this still blocked?
We are using fail2ban to block IPs which do too many connections with wrong credentials.
See this faq[faq](https://doc.eresearch.unige.ch/hpc/faq#when_i_tried_to_connect_to_the) for more information.


## Post 3 by @elia.cenci (2026-07-09T08:12:09.484Z)

Hello,
Yes, it is still blocked. Is there any diagnostics I shall run to find out what the problem is?
I have tried the following with corresponding output:
Host: login1.yggdrasil.hpc.unige.ch
Resolved IP: 129.194.9.68
Command: nc -G 5 -vz login1.yggdrasil.hpc.unige.ch 22
Result: Operation timed out
ping result included: Communication prohibited by filter from isdcrouter.isdc.unige.ch (129.194.168.1) Source shown: 10.195.169.154
Best,
Elia


## Post 4 by @Yann.Sagon (2026-07-09T09:09:39.301Z)

elia.cenci:
> isdcrouter.isdc.unige.ch
I see! Until the security issue is solved by Astro IT, you need to be connected either by Wi-Fi or VPN from the Astro department to reach Yggdrasil.


## Post 5 by @elia.cenci (2026-07-09T13:10:49.124Z)

Thank you, with the WiFi it works. I was connected via the cable, as we were instructed to do.
