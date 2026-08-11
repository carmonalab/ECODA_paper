# [tutorial] 🔐 SSH login troubleshooting (key mismatch)

- Source: https://hpc-community.unige.ch/t/4339

- Created: 2026-07-08T07:00:45.345Z

- Posts: 1

- Category: 13

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yann.Sagon (2026-07-08T07:00:45.498Z)

# :locked_with_key: SSH login troubleshooting (key mismatch)
If you cannot connect to the Baobab cluster via SSH, it is often due to a mismatch between:
- the private key on your machine
- the public key stored in LDAP
You can easily verify this by comparing their fingerprints.
---
## :white_check_mark: Step 1 — Check your local key
On your computer, run:
```
ssh-keygen -lf ~/.ssh/id_rsa
```
(Replace `id_rsa` with `id_ed25519` or your key file if different.)
This will display a fingerprint like:
```
2048 SHA256:XXXX... yourname@yourhost
```
---
## :white_check_mark: Step 2 — Check the key stored in LDAP
If you cannot SSH, connect via Open OnDemand:
:backhand_index_pointing_right: https://openondemand.baobab.hpc.unige.ch[https://openondemand.baobab.hpc.unige.ch](https://openondemand.baobab.hpc.unige.ch)
Open a shell and run:
```
/usr/bin/sss_ssh_authorizedkeys <your_username> | ssh-keygen -lf -
```
---
## :white_check_mark: Step 3 — Compare fingerprints
Compare the `SHA256:...` values from both commands:
- :white_check_mark: If they are identical → your key is correct
- :cross_mark: If they are different → you are not using the right private key. Check in my-account.unige.ch that the key listed there is correct. If this is the case, contact us.
---
## :magnifying_glass_tilted_right: Optional debugging
You can also run:
```
ssh -vvv login.baobab.hpc.unige.ch
```
---
## :light_bulb: Common issues
- Wrong key file used
- Key not loaded in `ssh-agent`
- Multiple keys in `~/.ssh/`
- Incorrect SSH config (`~/.ssh/config`)
---
If you’re still stuck, feel free to post on this forum with:
- your command output
- the error message
:backhand_index_pointing_right: It will help others (and us) diagnose the issue quickly.
