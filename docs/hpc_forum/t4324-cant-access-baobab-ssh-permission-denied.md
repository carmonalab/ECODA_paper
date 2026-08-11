# Can't access baobab - SSH permission denied

- Source: https://hpc-community.unige.ch/t/4324

- Created: 2026-06-26T08:15:09.394Z

- Posts: 2

- Category: 1

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Nicolas.Chiaruttini (2026-06-26T08:15:09.457Z)

Hello,
I tried for the first time to access the cluster, but I do not manage.
Here’s what I did:
- Created a ssh key
- Put the public key on my portal
- Started ssh agent on windows (my OS) and made sure the identity was set up
- Waited an hour
- Tried to log in
Here’s how it looks:
```
PS C:\\Windows\\system32> Start-Service ssh-agent
PS C:\\Windows\\system32> ssh-add -l
The agent has no identities.
PS C:\\Windows\\system32> ssh-add C:\\Users\\chiarutt.ssh\\id_ed25519
Enter passphrase for C:\\Users\\chiarutt\\.ssh\\id_ed25519:
Identity added: C:\\Users\\chiarutt\\.ssh\\id_ed25519 (isis\\chiarutt@CD-EM01EREV)
PS C:\\Windows\\system32> ssh -l chiarutt login1.baobab.hpc.unige.ch
chiarutt@login1.baobab.hpc.unige.ch: Permission denied (publickey,hostbased).
PS C:\\Windows\\system32> ssh login1.baobab.hpc.unige.ch
isis\\chiarutt@login1.baobab.hpc.unige.ch: Permission denied (publickey,hostbased).
PS C:\\Windows\\system32>
```
Should I use putty ? Is there a chance this can solve the issue ?
( I never accessed hpc yet but I should have a valid account)
Side note: why the heck is the code block formatting not working ? What did I do wrong ?
admin edit: you had “`\`” in front of each back tick, this is why the code block formatting wasn’t working.


## Post 2 by @Nicolas.Chiaruttini (2026-06-26T08:39:52.677Z)

Ah, nice, it works now!
@Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon) did you do anything to fix that ? Or was it just time ?
Thanks!
