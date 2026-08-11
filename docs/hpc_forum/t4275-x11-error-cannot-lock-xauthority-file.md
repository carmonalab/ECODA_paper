# X11 error: cannot lock .Xauthority file

- Source: https://hpc-community.unige.ch/t/4275

- Created: 2026-04-15T08:40:28.718Z

- Tags: bamboo

- Posts: 4

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Ghasem.Hajianfar (2026-04-15T08:40:28.804Z)

Hello everyone,
I am encountering the following error when trying to run jobs with X11 forwarding:
```
/usr/bin/xauth: error in locking authority file /home/users/h/user/.Xauthority
```
Has anyone experienced a similar problem or knows how to resolve it?
Thank you in advance for your help.
Best regards,
Ghasem


## Post 2 by @Gael.Rossignol (2026-04-16T09:30:02.557Z)

Ghasem.Hajianfar:
> `/home/users/h/user/.Xauthority`
Dear Ghasem,
It may be possible path is not correct :
```
/home/users/h/hajianf3/.Xauthority
```
What method are you using to have this result? (sbatch, script, …)?
Best regards,


## Post 3 by @Ghasem.Hajianfar (2026-04-16T09:49:43.178Z)

Thank you for your response.
I am trying to connect to Bamboo using the VS Code Remote extension, but I get the following error:
Screenshot from 2026-04-16 11-11-40
Screenshot from 2026-04-16 11-11-40966×337 14.8 KB
[Screenshot from 2026-04-16 11-11-40966×337 14.8 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/6a358c480dacee1f7f9b131b50f53a960df5202d.png)
I also tried connecting from the terminal using:
> ssh -Y hajianf3@login1.bamboo.hpc.unige.ch
I was able to connect to the server, but the following error appeared:
> /usr/bin/xauth: error in locking authority file /home/users/h/hajianf3/.Xauthority


## Post 4 by @Gael.Rossignol (2026-04-16T14:43:19.347Z)

Hello,
As we check it was a problem of quota.
Best regards,
