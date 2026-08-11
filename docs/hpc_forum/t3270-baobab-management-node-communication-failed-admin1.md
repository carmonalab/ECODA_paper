# Baobab Management node communication failed: admin1

- Source: https://hpc-community.unige.ch/t/3270

- Created: 2024-01-24T16:05:09.616Z

- Posts: 5

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Anton.Hanke (2024-01-24T16:05:09.663Z)

## Primary informations
Username: hankea + piasentin
Cluster: Baobab (ex: Yggdrasil)
## Description
High usage of Head node


## Post 2 by @Yann.Sagon (2024-01-24T16:35:05.817Z)

Dear @Anton.Hanke[@Anton.Hanke](https://hpc-community.unige.ch/u/anton.hanke) I don’t understand the meaning of the title. Users don’t have access to a server named admin1.
About the description, I guess you are talking about login2.baobab? Indeed a user was running a computation directly on the login node instead of using a compute node. We’ve killed the scripts and notified the user.
Feel free to comment if the issue was something else.
Best
Yann


## Post 3 by @Anton.Hanke (2024-01-24T16:48:47.031Z)

I am not sure actually where this message was coming from.
The usage of login2 by the scripts was so high, that labmates and me could not open new ssh sessions to the node and the one that worked for me only printed
the message:
Managment node communication failed: admin1
directly after login.
So yes killing that script removed the problems for me.
Thank you :slight_smile:


## Post 4 by @Yann.Sagon (2024-01-25T09:04:17.758Z)

Hi, thanks for the explanation, I never saw this message. But yes probably related with the high load of login2.
Great it is working now! And for anyone reading this message: remember it is forbidden to launch jobs directly on the login node (see cluster usage[cluster usage](https://www.unige.ch/eresearch/en/services/hpc/terms-use/)).


## Post 5 by @maciej.falkiewicz (2024-01-26T09:44:51.662Z)

@Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon) don’t we have the same situation again? I can ssh but bash is not launching.
Not able to rsync files to scratch as well.
Kind regards,
Maciej Falkiewicz
