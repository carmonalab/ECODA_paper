# [Baobab] Connection timed out

- Source: https://hpc-community.unige.ch/t/3686

- Created: 2024-10-14T13:18:50.447Z

- Tags: baobab

- Posts: 3

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Camille.Christe (2024-10-14T13:18:50.447Z)

ssh christec@login1.baobab.hpc.unige.ch


## Post 2 by @Adrien.Albert (2024-10-14T13:40:16.830Z)

@Camille.Christe[@Camille.Christe](https://hpc-community.unige.ch/u/camille.christe)
I have sent you an email to resolve this issue on Zoom.
Best Regards,


## Post 3 by @Adrien.Albert (2024-10-14T15:20:27.162Z)

## Update:
During the zoom meeting we could see the following error:
```
ssh: connect to username@login1.baobab.hpc.unige.ch port 22: Connection timed out
```
This error suggests that the login1.baobab address is resolved and accessible, but the connection has not been established.
Suggestion: the problem may be that the organisation outside unige is blocking flows to login1.baobab.
Pending changes to this organisation’s network service, Baobab is  accessible via OpenOnDemand, UNIGE VPN or from other clusters.
This issue is related to : :question: FAQ: I tried to connect without success.[FAQ: I tried to connect without success.](https://doc.eresearch.unige.ch/hpc/faq#i_tried_to_connect_without_suc)
