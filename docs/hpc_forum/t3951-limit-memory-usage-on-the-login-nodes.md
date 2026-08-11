# Limit Memory usage on the login nodes

- Source: https://hpc-community.unige.ch/t/3951

- Created: 2025-05-15T07:53:35.317Z

- Posts: 2

- Category: 8

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Lucille.Delisle1 (2025-05-15T07:53:35.361Z)

Dear HPC,
Regularly the front end of all three clusters are slow or not accessible anymore because there are heavy processes running on the front ends.
Would it be possible to set up a limit using cgroup that would kill any process run by a standard user using more than 8GB of RAM (or any amount that seems reasonable to you)? This is what is set up on the HPC instance I was using previously.
This way people that ran heavy commands on the login node thinking they were on an interactive job will realize their mistake and will save you time to notify them by email and reboot the login nodes.
Best,
Lucille


## Post 2 by @Gael.Rossignol (2025-06-20T07:34:36.608Z)

Dear Lucille,
We recently deploy a script to limit ressources on login nodes. Thanks you so much for helping us regarding that.
Best regards,
