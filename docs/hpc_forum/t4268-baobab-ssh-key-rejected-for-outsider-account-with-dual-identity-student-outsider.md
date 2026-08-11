# [Baobab] SSH key rejected for Outsider account with dual identity (student + outsider)

- Source: https://hpc-community.unige.ch/t/4268

- Created: 2026-04-02T09:10:57.943Z

- Tags: baobab

- Posts: 3

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Stephane.Bernardini (2026-04-02T09:10:58.095Z)

Primary informations
Cluster: baobab (yggdrasil and bamboo work fine with the same key)
Description
My Outsider account authenticates successfully via SSH public key
on both Yggdrasil and Bamboo, but the same key is rejected on
Baobab login1.
On Yggdrasil, `sss_ssh_authorizedkeys <username>` correctly returns
my registered ed25519 key. On Baobab, the key is rejected (ssh -vvv
shows type 51 after offering the key), SSH falls back to
keyboard-interactive (password prompt), and the connection is
eventually closed.
I suspect this may be related to the dual identity issue mentioned
in the March 2025 changelog (OpenOnDemand authentication fix
affecting users with Collaborator/Student dual identities). My
account has both an Outsider identity and a student identity, both
visible in the sss_ssh_authorizedkeys output on Yggdrasil.
Could you check whether sss_ssh_authorizedkeys returns my key on
login1.baobab, and whether my account is correctly provisioned on
that cluster?
Thank you


## Post 2 by @Yann.Sagon (2026-04-02T09:27:38.682Z)

Dear @Stephane.Bernardini[@Stephane.Bernardini](https://hpc-community.unige.ch/u/stephane.bernardini) you aren’t registered as outsider, but as student. Do you confirm?
Anyway, it seems we had an issue during the creation of your account on Baobab. This is now solved and your ssh key should work.
Best regards


## Post 3 by @Stephane.Bernardini (2026-04-03T13:10:46.642Z)

Hi @Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon),
Thank you for your help and for resolving the account issue on Baobab. Yes, I confirm I’m registered as a student (not an outsider), and my SSH key access is now working perfectly.
Best regards,
Stéphane Bernardini
