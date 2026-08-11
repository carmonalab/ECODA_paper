# 🔐 Authentication Change on Compute Nodes

- Source: https://hpc-community.unige.ch/t/4334

- Created: 2026-07-01T15:13:01.578Z

- Posts: 3

- Category: 6

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yann.Sagon (2026-07-01T15:13:01.707Z)

## Hello everyone,
We would like to inform you about an upcoming change regarding SSH authentication on compute nodes.
### :cross_mark: What is changing
- It will no longer be possible to add your SSH public key directly to the `~/.ssh/authorized_keys` file on compute nodes.
### :white_check_mark: What is not changing
- No changes if you connect:
  - from a login node to a compute node
  - using the usual mechanisms (interactive jobs, batch jobs, etc.)
### :rocket: Why this change?
This update enables:
- :white_check_mark: the use of `ProxyJump`
- :white_check_mark: direct connection to compute nodes
- :white_check_mark: the ability to establish SSH tunnels (e.g. for Jupyter, visualization, etc.)
### :key: New host key
- A new ED25519 host key has been added alongside the existing host key.
- You may be prompted to trust this new key the first time you connect after the change.
- The signature is `SHA256:R/cy4lk5x8qKwmrIq8R9tiRdneDtorBnqzEynx8OnGI` and is indicated in our doc[doc](https://doc.eresearch.unige.ch/hpc/access_the_hpc_clusters#ssh_host_keys_verification) too.
### :light_bulb: Example using ProxyJump
In your `~/.ssh/config` file:
```
Host login
    HostName login1.<cluster>.hpc.unige.ch
    User <your_username>

Host compute-node
    HostName <compute-node>
    User <your_username>
    ProxyJump login
```
Then connect directly with:
```
ssh compute-node
```
### :date: Effective date
This change will take effect on July 2nd at 10:00.
### :compass: Notes
- This change aims to simplify and standardize access while improving usability.
- If you have scripts or workflows relying on `authorized_keys` on compute nodes, they will need to be adapted.
---
:backhand_index_pointing_right: Feel free to reply to this post to share your use cases or ask questions.


## Post 2 by @Yann.Sagon (2026-07-02T09:24:29.929Z)

The change is now effective. Feel free to answer to this post if you notice an issue.
Best regards


## Post 3 by @Yann.Sagon (2026-07-02T13:29:51.032Z)

Update on authentication and data transfers between clusters
It is now once again possible to transfer data directly between the login nodes of the three clusters using scp or rsync, without having to store a private SSH key in your home directory.
As a result, users are requested to remove any SSH public and private keys that were previously copied to their home directory, as these are no longer required.
Examples:
```
rsync -avhP /path/source/ user@login1.yggdrasil:/path/destination/
```
or
```
rsync -avhP /path/source/ user@login1.yggdrasil:/path/destination/
```
Please take a moment to check your home directory and delete any SSH keys that may have been added for this purpose. Thank you for helping us maintain a secure environment.
