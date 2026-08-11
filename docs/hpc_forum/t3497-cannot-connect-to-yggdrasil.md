# Cannot connect to yggdrasil

- Source: https://hpc-community.unige.ch/t/3497

- Created: 2024-06-18T06:31:55.348Z

- Posts: 19

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Denis.Mongin (2024-06-18T06:31:55.388Z)

Dear HPC team
I don’t know if it is only me, but I tried to connect to yggdrasil (from windows using putty), same as I do for baobab using login1.yggdrasil.hpc.unige.ch, but the session times out.
Am I doing something wrong ?


## Post 2 by @Adrien.Albert (2024-06-18T09:41:38.231Z)

Hi @Denis.Mongin[@Denis.Mongin](https://hpc-community.unige.ch/u/denis.mongin)
When creating a post in the HPC Support > HPC Issue section, please use the provided template.
This helps ensure your post is clear, detailed, and understandable for everyone. It also ensures you include all relevant information.
Could you give the output of ssh using `-vvv` option ?


## Post 3 by @Denis.Mongin (2024-06-18T15:12:04.867Z)

```
OpenSSH_for_Windows_8.1p1, LibreSSL 3.0.2
debug3: Failed to open file:C:/Users/dmog/.ssh/config error:2
debug3: Failed to open file:C:/ProgramData/ssh/ssh_config error:2
debug2: resolving "login1.yggdrasil.hpc.unige.ch" port 22
debug2: ssh_connect_direct
debug1: Connecting to login1.yggdrasil.hpc.unige.ch [129.194.64.11] port 22.
debug3: finish_connect - ERROR: async io completed with error: 10060, io:000002E9367D4040
debug1: connect to address 129.194.64.11 port 22: Connection timed out
ssh: connect to host login1.yggdrasil.hpc.unige.ch port 22: Connection timed out
```


## Post 4 by @Adrien.Albert (2024-06-19T13:24:44.883Z)

Hi @Denis.Mongin[@Denis.Mongin](https://hpc-community.unige.ch/u/denis.mongin)
As you are using windows, could you try to use putty to connect on cluster:
https://doc.eresearch.unige.ch/hpc/access_the_hpc_clusters#from_windows[https://doc.eresearch.unige.ch/hpc/access_the_hpc_clusters#from_windows](https://doc.eresearch.unige.ch/hpc/access_the_hpc_clusters#from_windows)


## Post 5 by @Damien.Benis (2024-06-20T15:39:57.866Z)

Hello,
I cannot connect to the Yggdrasil cluster either, from X2GO or Filezilla.
Best
Damien


## Post 6 by @Giuseppe.Chindemi (2024-06-20T17:08:11.851Z)

Same with SSH… any news?


## Post 7 by @Genevieve.Savard (2024-06-20T19:08:44.105Z)

Same here, on the command line, keeps asking for password even though the good one is entered. On X2Go, it says “Connection refused”. My other computer was still logged on yggdrasil and jobs are still running and the cluster seems to be working… Looks like the trouble is just login access?


## Post 8 by @Adrien.Albert (2024-06-21T07:17:19.817Z)

Dear All,
It has been fixed. Could you please try again and let me know what happens?


## Post 9 by @Giuseppe.Chindemi (2024-06-21T08:20:27.569Z)

Thanks @Adrien.Albert[@Adrien.Albert](https://hpc-community.unige.ch/u/Adrien.Albert) ! I can connect via SSH and my jobs are still running, all good here.


## Post 10 by @Douglas.Stumpp1 (2024-06-21T13:27:17.676Z)

Hi, the issue was indeed resolved this morning but it is now back (both on x2go, ssh and fillezila).


## Post 11 by @Raphael.Rubino (2024-06-21T13:31:31.010Z)

The user “cabrerap” is currently using Yggdrasil login node as a compute node, all CPUs, RAM and swap are full, thus the login node is unreachable (or very slow).


## Post 12 by @Maura.Brunetti (2024-06-21T13:57:50.355Z)

Bonjour, je vois que yggdrasil ne repond toujours pas… est-ce qu’on peut faire quelque chose pour contacter `cabrerap’ et liberer le login node?


## Post 13 by @Douglas.Stumpp1 (2024-06-21T14:06:33.130Z)

Hi, can someone from the hpc team would be able to kill his job as he’s neither able to connect ? Thanks in advance


## Post 14 by @Adrien.Albert (2024-06-21T14:12:51.001Z)

Dear @Raphael.Rubino[@Raphael.Rubino](https://hpc-community.unige.ch/u/raphael.rubino)
Raphael.Rubino:
> The user “cabrerap” is currently using Yggdrasil login node as a compute node, all CPUs, RAM and swap are full, thus the login node is unreachable (or very slow).
Thank you for confirming my suspicions. However, please contact us directly if you are specifically targeting someone, rather than posting on hpc-community .


## Post 15 by @Adrien.Albert (2024-06-21T14:14:20.080Z)

Hi @Douglas.Stumpp[@Douglas.Stumpp](https://hpc-community.unige.ch/u/douglas.stumpp)
Login node have been rebooted at: 2024-06-21T14:04:00Z


## Post 16 by @Baptiste.Leforestier (2024-06-21T15:10:42.898Z)

Hi,
Are there not any safeguards on the login nodes to prevent this kind of occurrence?
If not, I suggest this might be a solution to deal with this kind of problems in a more systematic and automatic way in the future:
GitHub[GitHub](https://github.com/CHPC-UofU/arbiter2)
### GitHub - CHPC-UofU/arbiter2: A daemon that uses cgroups to monitor and manage...[GitHub - CHPC-UofU/arbiter2: A daemon that uses cgroups to monitor and manage...](https://github.com/CHPC-UofU/arbiter2)
A daemon that uses cgroups to monitor and manage user behavior on login nodes - CHPC-UofU/arbiter2
It’s a tool developed specifically for this, from the University of Utah’s HPC center.
Have a great weekend all!


## Post 17 by @Adrien.Albert (2024-06-21T15:33:51.372Z)

Hi @Baptiste.Leforestier[@Baptiste.Leforestier](https://hpc-community.unige.ch/u/baptiste.leforestier),
We’re always on the lookout for new ways to improve our service. I’m putting it on the agenda for our next team meeting.
Thank you !
Best Regards,


## Post 18 by @Denis.Mongin (2024-06-26T12:44:15.537Z)

As stated in my first message, I was using putty
> I don’t know if it is only me, but I tried to connect to yggdrasil (from windows using putty), same as I do for baobab using login1.yggdrasil.hpc.unige.ch, but the session times out.


## Post 19 by @Denis.Mongin (2024-06-26T12:49:01.497Z)

I tried today, and still can not connect.
It is my first time, so maybe I am doing wrong.
using putty, from windows. It works perfectly on baobab. I am doing the same, just with login1.yggdrasil.hpc.unige.ch


## Post 20 by @Yann.Sagon (2024-12-09T15:27:50.160Z)

A post was split to a new topic: Cannot connect to login1.yggdrasil[Cannot connect to login1.yggdrasil](https://hpc-community.unige.ch/t/cannot-connect-to-login1-yggdrasil/3761)
