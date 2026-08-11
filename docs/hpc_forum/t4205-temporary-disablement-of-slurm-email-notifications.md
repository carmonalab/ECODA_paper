# 🔔 Temporary Disablement of Slurm Email Notifications

- Source: https://hpc-community.unige.ch/t/4205

- Created: 2026-01-26T10:18:07.159Z

- Tags: slurm

- Posts: 2

- Category: 6

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yann.Sagon (2026-01-26T10:18:07.262Z)

Some users enable email notifications when submitting Slurm jobs. This has recently resulted in extremely high volumes of emails being sent from our clusters. As an example, a single user generated more than 20,000 notification emails over the last weekend.
In several cases, we have observed thousands of emails per day sent by a single user. This situation causes significant issues for the mail infrastructure and has been repeatedly reported by the team responsible for it. While we have contacted the affected users multiple times, this approach has proven to be time‑consuming and not sustainable.
We are currently working on the implementation of a per‑user daily email quota (XX emails per day, quota to be defined) in order to limit such abuse.
:backhand_index_pointing_right: In the meantime, Slurm email notifications are disabled until further notice.
Thank you for your understanding.


## Post 2 by @Yann.Sagon (2026-01-27T16:41:17.833Z)

This is now done (see linked post for details)
