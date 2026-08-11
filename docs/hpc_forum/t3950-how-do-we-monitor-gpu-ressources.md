# How do we monitor GPU ressources?

- Source: https://hpc-community.unige.ch/t/3950

- Created: 2025-05-14T08:39:29.689Z

- Posts: 7

- Category: 13

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Thibaut.Chataing (2025-05-14T08:39:29.689Z)

Hi,
I am struggling to visualise the consumption of my node.
I run an AI model on a gpu node in baobab and the run can be really memory intensive. I would like to monitor it to optimize the application but I don’t understand the layout of the Grafana.
Do you have some time to explain it to me
Thanks,
Thibaut


## Post 2 by @Yann.Sagon (2025-05-14T13:55:47.060Z)

Did you checked this dashboard?
image
image887×443 33.6 KB
[image887×443 33.6 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/d1a5e004426599820f3421af5d55a41a34539a9a.png)
ps: there is a bug right now: only the Yggdrasil GPUs are listed, we’ll fix that soon.


## Post 3 by @Thibaut.Chataing (2025-05-15T12:09:26.077Z)

Okay,
The PS should explain why I was not able to find the baobab gpu I was interested in. Thanks.


## Post 4 by @remy.ressegaire (2025-05-16T07:30:01.522Z)

Hello,
The bug is now fixed, you should be able to see baobab gpus in this dashboard.
Remy


## Post 5 by @Thibaut.Chataing (2026-06-04T07:30:32.854Z)

Hi,
When I try to connect to grafana (https://monitor.hpc.unige.ch/dashboards[https://monitor.hpc.unige.ch/dashboards](https://monitor.hpc.unige.ch/dashboards))
I’ve got a timeout error (ERR_CONNECTION_TIMED_OUT).
I’m connected on the UNIGE wifi at the campus BIOTECH.
Can you help me please ?
Best,
Thibaut


## Post 6 by @Gael.Rossignol (2026-06-04T07:41:18.459Z)

Dear Thibaut,
The issue is related to the following post:
Security Update – Access to Yggdrasil[Security Update – Access to Yggdrasil](https://hpc-community.unige.ch/t/security-update-access-to-yggdrasil/4298) HPC Announce[HPC Announce](https://hpc-community.unige.ch/c/hpc-announce/6)
> Dear users, 
Following a recent security incident affecting servers in the Astronomy Department, access to the Yggdrasil cluster has been restricted. 
Access changes

login1.yggdrasil is now accessible only from the UNIGE network or via the UNIGE VPN.
External connections without VPN are no longer allowed.

External users
Users who cannot use the VPN may access Yggdrasil through the web portal: 
backhand_index_pointing_right https://openondemand.yggdrasil.hpc.unige.ch[https://openondemand.yggdrasil.hpc.unige.ch](https://openondemand.yggdrasil.hpc.unige.ch) 
Recommended actions (eff…
As the tool is hosted by the Astronomy Department, access has been restricted for some time.
Sorry for the inconvenience.
Best regards,


## Post 7 by @Thibaut.Chataing (2026-06-04T07:48:45.927Z)

Dear Gael,
Thanks for the answer.
Do you mean that I can’t access grafana ? Or is there some steps I need to follow ?
Because :
- I manage to connect to yggdrasil from the terminal through ssh.
- I’m connected to the UNIGE wifi
Best,
Thibaut
