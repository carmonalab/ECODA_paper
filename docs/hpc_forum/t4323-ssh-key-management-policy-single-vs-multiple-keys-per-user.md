# SSH Key Management Policy – Single vs Multiple Keys per User

- Source: https://hpc-community.unige.ch/t/4323

- Created: 2026-06-26T06:47:15.688Z

- Posts: 9

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yann.Sagon (2026-06-26T06:47:15.829Z)

## :locked_with_key: SSH Key Management Policy – Single vs Multiple Keys per User
Hi everyone,
We would like to open a discussion regarding our current SSH access policy on the cluster.
At the moment:
- SSH access is key-based only
- Each user can register a single public key via the web portal
- The private key remains entirely under the user’s responsibility
Some users have suggested that we should allow multiple SSH keys per user. We would like to gather broader feedback from the community before considering any changes.
---
## :light_bulb: Arguments in favor of multiple keys
A few motivations have been raised:
- Redundancy
  - Having multiple keys could prevent lockout if a private key is lost
- Isolation / compartmentalization
  - Users may want separate keys for different devices (laptop, workstation, CI/CD, etc.)
  - Or separate keys per usage context
---
## :balance_scale: Counterpoints and considerations
From our perspective, there are important trade-offs:
### 1. Redundancy is not strictly necessary
- If a private key is lost, generating a new key pair is fast and straightforward
- The user can then update their public key via the portal
- Therefore, the operational benefit of redundancy is relatively limited
### 2. Increased attack surface
- Allowing multiple keys increases the probability that:
  - One key is poorly protected
  - One key is unintentionally exposed
- Any compromised key would grant access to the cluster
### 3. Lack of fine-grained control
- Ideally, multiple keys should come with scoping mechanisms, such as:
  - Restricting a key to specific hosts
  - Restricting usage contexts
- However:
  - Our infrastructure does not currently support fine-grained key scoping
  - Cluster nodes are not managed via a centralized identity system like Active Directory
  - Implementing such controls would be complex
- Without proper scoping, multiple keys effectively behave as full-access credentials
---
## :wrench: Current exception
For completeness, note that we already allow multiple SSH keys in a specific case:
- Shared service accounts, typically used by multiple administrators
In this context:
- Multiple keys are required for operational reasons
- Key ownership is tied to individual administrators
- This use case is managed differently from standard user accounts
---
## :red_question_mark: Open questions
We would like your input on the following:
- Do you have concrete use cases that require multiple SSH keys?
- Would multiple keys significantly improve your workflows?
- How do you currently manage key security on your side?
- Would you expect additional controls (e.g. per-key restrictions)?
---
## :compass: Possible directions
Depending on feedback, we could explore:
- Keeping the current single key model
- Allowing multiple keys with limitations
- Defining best practices for key management instead
---
We look forward to hearing your thoughts!
Thanks in advance for your feedback :folded_hands:


## Post 2 by @Yann.Sagon (2026-06-29T11:24:43.548Z)

A post was merged into an existing topic: Mandatory SSH Key Authentication for Baobab HPC Access[Mandatory SSH Key Authentication for Baobab HPC Access](https://hpc-community.unige.ch/t/mandatory-ssh-key-authentication-for-baobab-hpc-access/4318/9)


## Post 3 by @Martin.Voelkle (2026-06-29T07:12:24.288Z)

Have you heard about our lord and savior Kerberos? It would sidestep the whole issue without needing a special management interface and would not suffer from the increased attack surface.


## Post 4 by @Lucille.Delisle1 (2026-06-29T09:47:28.820Z)

I can share how I manage key security on my side. I have a single key.
I have a ssh-key for my laptop with a passphrase (generated randomly). I am using KeePassXC to add it to the ssh agent and remove it when I lock my computer. This is very handy because as soon as I open KeePassXC I do not need to pay attention to my ssh key passphrase.
I am using linux and I do not know how easy it is to use this setting in windows/mac. I can give more details if people are interested.


## Post 5 by @Yann.Sagon (2026-06-29T11:24:09.749Z)

Martin.Voelkle:
> Have you heard about our lord and savior Kerberos?
Yes this would be interesting, it is similar as having a single SSH key managed by DiSTIC.


## Post 6 by @Frederic.Masclaux (2026-07-01T12:59:23.204Z)

Having a single key is perfect when using a single computer and a single connection (like a single cluster). It is probably the real topic of this discussion : what is convenient for Unige cluster is maybe not recommended overall.
With several computers, the key must be copied from one computer to another. How to do the copy pratically? USB drive? Email? If one of these computers is compromised, how to identify it? What about key revocation?
With several connections/clusters, then one key gives access to all. Even if they are completely unrelated. Even if 2 of these connections are never accessed from the same computer.
When comes the time to renew the SSH key, then it is for all connections at once. What about key rotations?


## Post 7 by @Yann.Sagon (2026-07-01T13:23:23.447Z)

Frederic.Masclaux:
> With several computers, the key must be copied from one computer to another. How to do the copy pratically? USB drive? Email?
Are you talking about the public key? If yes, the public key is read directly by ldap, we don’t copy it locally.
Frederic.Masclaux:
> When comes the time to renew the SSH key, then it is for all connections at once. What about key rotations?
If you read the public key from ldap, this isn’t an issue.


## Post 8 by @Frederic.Masclaux (2026-07-01T14:44:14.323Z)

The private key should be copied from one worskstation (laptop, desktop) to another.
The ldap is for Unige only, not for the other possible connections/clusters


## Post 9 by @Yann.Sagon (2026-07-01T14:55:05.968Z)

Frederic.Masclaux:
> The ldap is for Unige only, not for the other possible connections/clusters
Exactly. LDAP only applies to Unige-managed systems.
For access to non-Unige servers or clusters, users are free to manage their own SSH key pairs independently (i.e. outside of the central registry), according to their needs,


## Post 10 by @Frederic.Masclaux (2026-07-01T16:01:33.676Z)

One key per context is
Yann.Sagon:
> For access to non-Unige servers or clusters, users are free to manage their own SSH key pairs independently (i.e. outside of the central registry), according to their needs,
Using one specific key per context/connection is clearly a good practice, decreasing the risk for other connections/clusters than Unige. If the “specific keys” method is really applied by the user, it solves the biggest security problem. The problem remains that using a copied single key encourages the bad practices, prevent device isolation and increase exposure.
