# Resuming paused jobs

- Source: https://hpc-community.unige.ch/t/3256

- Created: 2024-01-18T08:23:05.928Z

- Posts: 2

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Debajyoti.Sengupta (2024-01-18T08:23:05.978Z)

Hello,
I have a relatively large job array that I launched yesterday that was still running.
However, I have another job that required a higher priority, and as such I used `scontrol hold jobid` to pause those jobs, so the new jobs could go through.
Currently, though, I am unable to resume those held jobs. Is this the expected behaviour?
For posterity, if such a situation does arise, what is the recommended procedure?
Relevant jobid that was paused: `6767423_`


## Post 2 by @Debajyoti.Sengupta (2024-01-18T08:26:07.028Z)

Sometimes you find the answer immediately after you post a question in the forum.
Thanks @John.Raine[@John.Raine](https://hpc-community.unige.ch/u/john.raine) for the help. I can resume the held jobs with `scontrol release jobid`
