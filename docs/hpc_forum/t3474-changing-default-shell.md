# Changing default shell

- Source: https://hpc-community.unige.ch/t/3474

- Created: 2024-06-05T07:39:59.689Z

- Posts: 3

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Vilius.Cepaitis (2024-06-05T07:39:59.722Z)

Hello,
Whenever I try changing my default shell to zsh, it works fine:
> (baobab)-[cepaitis@login2 ~]$ chsh
> Changing shell for cepaitis.
> New shell [/bin/bash]
> /bin/zsh
> Password:
> Shell changed.
However, after a few weeks, it changes back by itself.
> (baobab)-[cepaitis@login2 ~]$ echo $SHELL
> /bin/bash
Do you know what might be going on?
Best regards,
Vilius


## Post 2 by @Vilius.Cepaitis (2024-06-26T13:45:40.401Z)

Hello,
Just to say that I am still facing this issue.
Cheers,
Vilius


## Post 3 by @Gael.Rossignol (2024-06-26T19:40:45.638Z)

Vilius.Cepaitis:
> /bin/zsh
Dear Vilius,
Unfortunatly you can’t change the shell by your own way.
We have a centralised reference for all users and we deploy users each time we reinstall a server or a full cluster.
Now I have done the change for your user, so you will use zsh as default shell.
Best regards,
