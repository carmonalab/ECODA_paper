# Scratch: No such file or directory

- Source: https://hpc-community.unige.ch/t/4201

- Created: 2026-01-22T14:44:14.871Z

- Tags: baobab

- Posts: 4

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Davide.Piras (2026-01-22T14:44:14.923Z)

Hello team, it seems that for the last few days the `scratch` folder is not accessible on baobab. Is it just me, or is it a more general issue?
I cannot `cd` into it, nor read files which I’m quite sure are there.
Thank you very much!
Best wishes,
Davide


## Post 2 by @Gael.Rossignol (2026-01-23T10:55:04.478Z)

Dear Davide,
Scratch is working fine on baobab no issue reported, after checking your account it seems you remove probably by mistake the scratch simlink in your home.
To complete answer you can also use the true path “/srv/beegfs/scratch/users/p/piras/” and your data are available.
If you need I can recreate simlink for you, let me know.
Best regards,


## Post 3 by @Davide.Piras (2026-01-23T11:31:35.379Z)

Gael.Rossignol:
> /srv/beegfs/scratch/users/p/piras/
Dear Gael, thank you for finding the problem! Could you please recreate the simlink for me?
Thanks again!
Davide


## Post 4 by @Gael.Rossignol (2026-01-26T07:42:41.974Z)

Hello,
Scratch simlink is now present.
Best regards,
