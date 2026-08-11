# Source: https://doc.eresearch.unige.ch/hpc/access_the_hpc_clusters
# Snapshot: 2026-08-11
# Crawled: 2026-08-11T14:32:29Z

---

## Access the clusters

### Account

To access the HPC clusters (Baobab/Yggdrasil/bamboo) and submit jobs, you need a valid HPC account.

Access is reserved for:
- Members of the University of Geneva (Unige)
- Members of HES-SO GE
- External collaborators
- Guest

> 
Each HPC account in linked to a Primary Investigator (PI)(=repondant = responsable (different terms but same purpose)). A PI  is a member of the University of Geneva (Unige) who invites and takes responsibility for an user’s access to the HPC service. The PI(repondant) ensures that the access is justified, appropriate, and compliant with university policies.

#### Standard Account

This is the default account type for Unige members. To request an account please follow fill the form on DW: https://dw.unige.ch/openentry.html?tid=hpc

To connect to the HPC cluster, you must authenticate using An SSH key registered in your [ISIS account](https://my-account.unige.ch/)

#### External Account

An **external account** is intended for collaborators from other institutions who are working closely with Unige researchers. (if you only need an HPC access [Outsiders account](https://doc.eresearch.unige.ch/hpc/access_the_hpc_clusters#outsider_account) is more appropriated)

**How to request one:**
- Ask your [department administrator](https://memento.unige.ch/download_file/230/374) to register the external collaborator in the university database: [https://catalogue-si.unige.ch/en/isis](https://catalogue-si.unige.ch/en/isis)
- More info: [https://memento.unige.ch/doc/0214/](https://memento.unige.ch/doc/0214/)
- Once registered, follow the usual HPC account request procedure as described above.

#### Outsider Account

An **Outsider** is a person external to Unige who has been invited by a **PI** specifically to use the HPC service. This lightweight account is limited to HPC access and avoids standard administrative procedures. (24h waiting time to get an account)

**Key details:**
- Authentication is done via the SSH Public Key (only) provided during the invitation process  (Authentication with password is disabled)
- SSH key can be updated at any time: [Update your SSH key](https://applicant.unige.ch/main/outsider-info/update)
- SSH key updates are applied daily at 1:30 PM and 5:00 AM via the UNIGE Active Directory

Only PI can invite and manage Outsiders, visit: https://gestion-externe.unige.ch/main/outsider-requests

(**Access requires approval — please contact us with a short motivation if you'd like it enabled.**)

Once access is granted, you will be able to create an invitation, which requires:

- The email address of the future Outsider
- Selection of the appropriate service (**High-Performance Computing**)
- Setting an expiration date (maximum 1 year)
- An optional note for the guest
- Acceptance of the terms of use

Once invited, the future Outsider will receive an email with detailed instructions to finalize their registration.
(**Tip: read it carefully!**)

**To renew an expired Outsider account, a new invitation must be created.**

### Inactivity Notice and Account Deletion Policy

Accounts that have been inactive for a period of **one year** will be flagged for deletion. This is to ensure that we are not storing unused data and to avoid data being left without an owner in the event that the account is deleted from the UNIGE central user directory.

Users will receive an email notification prior to the scheduled deletion of their account, giving them the opportunity to log in and keep their account active. If no response is received within one month, the account will be deleted, along with all associated data.

**Users no more at Unige:**

In the middle of each month, an automated process is triggered to identify all accounts that are no longer listed in the University of Geneva (UNIGE) electronic directory.
An account that is no longer listed indicates that the individual has left UNIGE.

These accounts are automatically flagged for deletion in the current end-of-month deletion batch.
At the end of the month, the account and all associated data are permanently deleted.

__**Example:**__

- **15 January 2025**: Check Account
    - Dark Vador is no longer listed in the UNIGE directory. => Tagged for deletion

- **31 January 2025**: Batch deletion
    - Dark Vador’s account and associated data are removed.
### Cluster connection
Once you have received an email confirming the creation of you account, you have access to our clusters: Baobab and Yggdrasil.

You can connect to the HPC clusters only through the login nodes.
The clusters are reachable from outside UNIGE as well, without the need to use a VPN.

#### login nodes

- **For Baobab** : `login1.baobab.hpc.unige.ch`
- **For Yggdrasil** : `login1.yggdrasil.hpc.unige.ch`
- **For Bamboo** : `login1.bamboo.hpc.unige.ch`

#### Connect using SSH

You can access the clusters from anywhere through `ssh` with your [ISIs](https://catalogue-si.unige.ch/isis) account and SSH key.

> It is mandatory to access the cluster using your SSH key, it isn't allowed anymore to connect using your password.

##### ssh PublicKey

1. For security measure it isn't possible to add your ssh key to AuthorizedKeyFile in the login node. This has been disabled to prevent any non registered at UNIGE to log In.
2. You ssh public key must be added to your UNIGE account profile. (it's like if the AuthorizedKeyFile is bonded to your UNIGE identity )
3. The UNIGE central directory synchronizes the ssh key every 5 min.
4. More information on [https://hpc-community.unige.ch/t/authentication-modification-sshpublickey-managment/3011](https://hpc-community.unige.ch/t/authentication-modification-sshpublickey-managment/3011)

**1. Generate your ssh-key** (It is mandatory to create it with password/passphrase for more security) by following
[this documentation (multi-platform)](https://hpc-community.unige.ch/t/tutorial-create-and-register-an-ssh-key-to-access-hpc-clusters/4343)

**2. Upload your public key** to your Isis profile by updating "My SSH public key" on:
- [my-account](https://my-account.unige.ch) for UNIGE/external users.
- [applicant](https://applicant.unige.ch/main/outsider-info/update) for Outsider users.

**Note**: Make sure you copy the public ssh key linked to the private key you're going to use. If you have regenerated your ssh key, you'll need to put your public key in [my-account](https:*my-account.unige.ch) or [applicant](https:*applicant.unige.ch/main/outsider-info/update)

##### multiple ssh key
> Multi key isn't supported anymore, but could be in the future.
It is possible to register multiple SSH public keys on the authentication server. However, [my-account.unige.ch](https://my-account.unige.ch/) does not allow this at the moment (work in progress). In the meantime, please send your request to the LDAP team directly at dl-distic-windows-team@unige.ch.

After requesting an additional sshPublicKey, if you update it via my-account, all previous references will be overwritten.

##### SSH Host Keys Verification

> 
We are deploying new host keys in ed25519 format!

The first time you connect, your SSH client will ask you to confirm the **server fingerprint**.
The format may differ depending on your SSH client (MD5, SHA256, etc.).

###### Current fingerprints

| Cluster / Host | Key Type | Fingerprint (MD5) | Fingerprint (SHA256) |
| --- | --- | --- | --- |
| Baobab / login nodes | RSA | MD5:8f:75:c4:18:8a:75:f1:f1:19:4d:85:92:3b:b6:2a:e1 | SHA256:tKqp4nljL+EGVKl8T0VF2nS36DkHVFMpLxQOPg/gKvg |
| Baobab / login nodes | ED25519 | MD5:23:7a:4f:a9:c0:5d:41:43:1a:b3:a8:c2:7f:30:32:29 | SHA256:R/cy4lk5x8qKwmrIq8R9tiRdneDtorBnqzEynx8OnGI |

###### Important

Fingerprints are provided **for information only**.

- Users **do not need to compute them manually**.
- It is **not possible to compute the server fingerprint before connecting**, since the host key is retrieved during the first SSH connection.
- The SSH client automatically displays the fingerprint when you connect for the first time.
- You only need to **compare what is displayed with the values in the table above**.

###### How to compute fingerprints (for administrators)

MD5 format

```console
ssh-keygen -l -E md5 -f /etc/ssh/ssh_host_ed25519_key.pub
256 MD5:23:7a:4f:a9:c0:5d:41:43:1a:b3:a8:c2:7f:30:32:29 root@admin1.roseau (ED25519)

```

SHA256 format (default modern SSH)

```console
ssh-keygen -l -f /etc/ssh/ssh_host_ed25519_key.pub
256 SHA256:R/cy4lk5x8qKwmrIq8R9tiRdneDtorBnqzEynx8OnGI (ED25519)

```

###### X2Go / alternative fingerprint (SHA1)

Some clients (e.g., X2Go on Windows) display fingerprints in SHA1 format.

```console
for x in /etc/ssh/*.pub; do
  echo $x
  cut -d' ' -f2 < $x | base64 -d | openssl sha1 -c
done

```

Example:
```console
/etc/ssh/ssh_host_rsa_key.pub
(stdin)= 67:03:fc:6f:32:7c:19:9b:97:b9:e8:7b:12:1d:ad:a6:7b:c9:4c:9c

```

###### Notes

- ED25519 keys are now the **recommended and default**.
- Different clients display fingerprints in different formats (MD5, SHA256, SHA1).
- Always verify fingerprints on first connection to avoid MITM attacks.

If you type your username/password **wrong 3 times** in a row, you will be **banned for 15 minutes** before you can try again. Please read more in the [Troubleshooting](https://doc.eresearch.unige.ch/hpc/access_the_hpc_clusters#troubleshooting) section.

##### From Linux and Mac OS

Connect to the login node from the terminal:

```console
ssh youruser@clusterhostname

```

> add `-Y` if you need X forwarding
```console
ssh -Y youruser@clusterhostname

```

N.B. : replace `clusterhostname` with the [login node](https://doc.eresearch.unige.ch/hpc/access_the_hpc_clusters#login_nodes) of the cluster you are willing to use.

##### From Windows

To access Baobab or Yggdrasil in SSH from Windows, you can use :
- ssh client integrated in PowerShell. Recent versions of Window 10 offer a ssh client. You can use it the same way you would use it on Linux (see [above](https://doc.eresearch.unige.ch/hpc/access_the_hpc_clusters#from_linux_and_mac_os))
- PuTTY which has been the *de facto* solution for years.
  - You can download it here [PuTTY](http://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html).
  - You will find the needed information on the screenshots below.

N.B. : replace `clusterhostname` with the [login node](https://doc.eresearch.unige.ch/hpc/access_the_hpc_clusters#login_nodes) of the cluster you are willing to use.

PuTTY create session :

PuTTY enter username :

PuTTY save session :

PuTTY open session :

Once you open a session with PuTTY, you will be prompted for your password.
```console
| Password:

```

Unlike Windows systems, Linux and Unix systems will **not display** **`*`** (star character) or any other character to indicate that you've entered something/anything in a password field - it simply stays totally blank as you type. Cursor won't blink, move, etc.

Just type your password and press enter, everything will be fine.

#### Access to the compute nodes

In case you need to access the compute nodes for purpose of debugging your software, using `htop` or other tool directly on the node, you need to go through the login node first and connect from there using ssh.

```
ssh cpu001

```

It is important to note that you **cannot** access directly the compute nodes unless you have a **RUNNING job** on it.
As soon as your job is finished, you will be logged out from the compute node.

### GUI access / Desktop with X2Go

It's possible to access Baobab or Yggdrasil using a graphical interface. We support only `X2Go`.

You need to install [X2Go client](http://wiki.x2go.org/doku.php/download:start) on your computer and create a session following the screenshots below prior to connect to Baobab.

X2Go create session :

Once connected to the Linux desktop, you can launch lightweight (image viewer, etc.) applications directly.

If you need to launch a heavy graphical program such as Matlab for example, create a new interactive
session (see the [Interactive jobs](https://doc.eresearch.unige.ch/hpc/slurm#interactive_jobs) section in the Slurm page)
on a compute node. It is [forbidden](https://doc.eresearch.unige.ch/hpc/best_practices#rules_and_etiquette) to use the login node for this purpose.

### File transfer

If you need to transfer files from your computer to Baobab or Yggdrasil, you need to use the `sftp` protocol (`scp`).

#### From Linux

We suggest you use `scp` or a program like `rsync`

#### From Windows

We suggest you to use [FileZilla client](https://filezilla-project.org/download.php?show_all=1) to transfer your files.

do not download the sponsored version (the name sponsored appears in the download link), this installer may include bundled offers that may be recognized as virus

Filezilla create session :

### SSH tunnel and socks proxy
If you want running JupyterLab or VScodeServer you may be interested by  [OpenOnDemand](https://doc.eresearch.unige.ch/hpc/how_to_use_openondemand)

The **login** nodes have a firewall that prevent incomming connection other than ssh.

If you need to access a service from the cluster, please follow the:

**1.** On your local machine, Save old ssh key and create a new one
```
$ mkdir ~/.ssh/old
$ mv ~/.ssh/*  ~/.ssh/old
$ ssh-keygen -t ed25519

```

**2.** Copy the rsa.pub in [https://my-account.unige.ch/main/home](https://my-account.unige.ch/main/home) (for Unige Account) [https://applicant.unige.ch/](https://applicant.unige.ch/) (for Outsider Account) and wait for 5 min the synchronisation with AD is done.
the following command on login node should print your public ssh key registered in the AD:

(baobab)-[alberta@login1 ~]$ /usr/bin/sss_ssh_authorizedkeys $USER
ssh-ed25519  [...]

**3.** On your local machine configure the proxyjump:
```
[alberta@localhost ~]$ cat .ssh/config_baobab

host baobab
   HostName login1.baobab.hpc.unige.ch
   User alberta

Host cpu*
   HostName %h
   User alberta
   ProxyJump baobab

Host gpu*
   HostName %h
   User alberta
   ProxyJump baobab

```

**4.** Alloc a test job and open a new tab on your local machine and try to connect on the allocated node:

**On baobab:**
```
(baobab)-[alberta@login1 ~]$ salloc --time=00:05:00
salloc: Pending job allocation 5574654
salloc: job 5574654 queued and waiting for resources
salloc: job 5574654 has been allocated resources
salloc: Granted job allocation 5574654
salloc: Waiting for resource configuration
salloc: Nodes cpu001 are ready for job

```

At the same time On your local machine, connect to the compute with selecting the right ssh config file (For this example: Baobab):
( If you never connected to the compute node, you need to accept the host key, you can check the fingerprint above in this doc)

```
[alberta@localhost .ssh]$ ssh -F .ssh/config_baobab cpu001
The authenticity of host 'cpu001 (<no hostip for proxy command>)' can't be established.
RSA key fingerprint is SHA256:tKqp4nljL+EGVKl8T0VF2nS36DkHVFMpLxQOPg/gKvg.
RSA key fingerprint is MD5:8f:75:c4:18:8a:75:f1:f1:19:4d:85:92:3b:b6:2a:e1.
Are you sure you want to continue connecting (yes/no)? yes
Warning: Permanently added 'cpu001' (RSA) to the list of known hosts.
Last login: Tue Oct 24 10:49:29 2023
Installed: Thu Aug 17 14:40:08 CEST 2023

```

> 
More Information on HPC-community forum:

[https://hpc-community.unige.ch/t/tutorial-ssh-tunneling-and-socks-proxy/1795](https://hpc-community.unige.ch/t/tutorial-ssh-tunneling-and-socks-proxy/1795)
[https://hpc-community.unige.ch/t/proxyjump-ssh-not-working-on-baobab/3078/15](https://hpc-community.unige.ch/t/proxyjump-ssh-not-working-on-baobab/3078/15)
