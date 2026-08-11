# [Solved] SSH login refused (publickey,hostbased) despite key on my-account

- Source: https://hpc-community.unige.ch/t/4322

- Created: 2026-06-24T13:34:01.629Z

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Hugues.Vinzant (2026-06-24T13:34:01.688Z)

Hi,
Since the switch to key-only authentication, I can no longer connect to Baobab
(account vinzant5, UNIGE student):
```
vinzant5@login1.baobab.hpc.unige.ch: Permission denied (publickey,hostbased).
```
What I’ve checked:
- My ed25519 key is registered on my-account.unige.ch.
- ssh -v shows the client offers the key, and the server rejects it.
It looks like my key hasn’t been propagated to authorized_keys on Baobab. Is this a known sync delay, or is there something I need to do?
Thanks!


## Post 2 by @Hugues.Vinzant (2026-06-24T13:47:48.819Z)

~1/2h after the SSH key registration on my-account I could connect so I guess it was just a matter of sync.
