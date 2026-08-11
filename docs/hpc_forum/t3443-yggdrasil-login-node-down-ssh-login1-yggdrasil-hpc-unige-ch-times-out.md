# Yggdrasil login node down [ssh login1.yggdrasil.hpc.unige.ch times out]

- Source: https://hpc-community.unige.ch/t/3443

- Created: 2024-05-10T07:10:06.152Z

- Posts: 5

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Arianna.Nigioni (2024-05-10T07:10:06.215Z)

## Primary informations
Username: nigioni
Cluster: Yggdrasil
## Description
I have been trying to log into Yggdrasil  through ssh since yesterday, 9/05, but there seems to be a problem with the connection because I do not get asked for a password and I get “Operation timed out” error.
## Steps to Reproduce
I type ssh nigioni@login1.yggdrasil.hpc.unige.ch and, after ~1 minute, I get the message ssh: connect to host login1.yggdrasil.hpc.unige.ch port 22: Operation timed out
## Expected Result
I was expecting to get asked for my password
## Actual Result
I got the error message ssh: connect to host login1.yggdrasil.hpc.unige.ch port 22: Operation timed out


## Post 2 by @William.Ceva (2024-05-10T08:43:39.808Z)

Can confirm I am having the same issue, as well as other people I have talked to.


## Post 3 by @Arianna.Nigioni (2024-05-10T09:12:08.831Z)

Now I am able to connect again!


## Post 4 by @William.Ceva (2024-05-10T09:31:54.804Z)

same for me as well now


## Post 5 by @Gael.Rossignol (2024-05-13T07:31:20.574Z)

Dear Users,
Server login1.yggdrasil has failled on friday because of heavy load. It has been rebooted and is now fully operationnal.
Bets regards,
