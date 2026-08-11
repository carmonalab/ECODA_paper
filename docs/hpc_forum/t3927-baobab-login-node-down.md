# Baobab login node down

- Source: https://hpc-community.unige.ch/t/3927

- Created: 2025-04-19T16:11:21.350Z

- Tags: baobab

- Posts: 7

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Bharathkumar.Radhakrishnan (2025-04-19T16:11:21.392Z)

Hi all,
I can ping baobab, but it does not seem to respond to ssh login. Can the admins please have a look? Thanks!


## Post 2 by @Anton.Hanke (2025-04-22T07:16:31.888Z)

It is the same for me this morning.
I can request login w/ ssh to the login node but am rejected after providing my password.
Similarly to when the cluster is in maintenance?
I can provide a ssh log if needed.
Best wishes,
Anton


## Post 3 by @Adrien.Albert (2025-04-22T08:05:52.663Z)

Hi @Bharathkumar.Radhakrishnan[@Bharathkumar.Radhakrishnan](https://hpc-community.unige.ch/u/bharathkumar.radhakrishnan)
We’re investigating the matter, but for the moment the login node has been restarted and is available again.


## Post 4 by @Adrien.Albert (2025-04-22T08:42:25.182Z)

We suspect a user attempted to open or modify a very large file, resulting in login1’s memory (512 GB RAM) being saturated. This resulted in a temporary system outage.
We have contacted the user concerned to better understand the situation and share best practices to avoid similar problems in the future.
We apologize for the inconvenience caused over the long weekend and thank you for your understanding.


## Post 5 by @Alejandro.PozasKerstjens (2025-04-22T09:05:55.923Z)

Hello. Is it normal that I can still not access Baobab? I get a timeout abort when I try to log in via ssh.


## Post 6 by @Adrien.Albert (2025-04-22T09:11:40.243Z)

Hello @Alejandro.PozasKerstjens[@Alejandro.PozasKerstjens](https://hpc-community.unige.ch/u/alejandro.pozaskerstjens)
Could you try again please


## Post 7 by @Alejandro.PozasKerstjens (2025-04-22T09:14:44.049Z)

Hello @Adrien.Albert[@Adrien.Albert](https://hpc-community.unige.ch/u/Adrien.Albert). I have just checked and now I can confirm that I can access. Thank you very much.
