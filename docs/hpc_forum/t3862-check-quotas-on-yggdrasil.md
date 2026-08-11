# Check quotas on yggdrasil

- Source: https://hpc-community.unige.ch/t/3862

- Created: 2025-03-13T12:55:08.373Z

- Tags: yggdrasil

- Posts: 5

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Thibault.Garel (2025-03-13T12:55:08.412Z)

Hi all,
I want to check my storage disk quota/usage for projects on yggdrasil under /srv/astro/projects/…
How do so please ?
The usual commands (e.g. beegfs-ctl --getquota --gid…) seem to fail with this error message :
> Unrecoverable error: No connAuthFile configured. Using BeeGFS without connection authentication is considered insecure and is not recommended. If you really want or need to run BeeGFS without connection authentication, please set connDisableAuthentication to true.
Thanks.


## Post 2 by @Adrien.Albert (2025-03-13T14:51:24.445Z)

Hello @Thibault.Garel[@Thibault.Garel](https://hpc-community.unige.ch/u/thibault.garel)
You can find an aswer in our FAQ: storage[FAQ: storage](https://doc.eresearch.unige.ch/hpc/faq#i_have_a_question_about_the_st) :slight_smile:


## Post 3 by @Thibault.Garel (2025-03-13T16:31:37.290Z)

Thanks for the reply, but unless I am missing smthg the documentation only explains how to get quotas for the home and the scratch, what about the projects under /srv/astro/projects/ ?


## Post 4 by @Yann.Sagon (2026-01-27T16:10:16.281Z)

@Thibault.Garel[@Thibault.Garel](https://hpc-community.unige.ch/u/thibault.garel) sorry, we missed you post Apologies for Unanswered Posts — Notification Settings Updated[Apologies for Unanswered Posts — Notification Settings Updated](https://hpc-community.unige.ch/t/apologies-for-unanswered-posts-notification-settings-updated/4203)
The space `/srv/astro/projects` isn’t managed by us but by astro IT team. We just mount is on the cluster. I’ll check with them what is going on.
Best regards
Yann


## Post 5 by @remy.ressegaire (2026-01-28T09:20:28.141Z)

Hello @Thibault.Garel[@Thibault.Garel](https://hpc-community.unige.ch/u/thibault.garel) ,
Which project would you like to check under /srv/astro/projects/ ?
Rémy
