# Unable to conect to node via ssh

- Source: https://hpc-community.unige.ch/t/4292

- Created: 2026-05-11T06:58:31.784Z

- Posts: 8

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Max.Briel (2026-05-11T06:58:31.856Z)

On yggdrasil, I normally proxyjump to an interactive node. This used to work fine until last Friday.
Currently, I am unable to connect to a specific node, even when I have an active job running on it. When using ssh to connect to the node from the login node, a connection fails as well.
Here’s a minimum working example for me.
- Login to login1
- Run:
```
> salloc

> ssh cpu00*
```
Then following error is returned:
```
Connection closed by 192.168.103.74 port 22
```
Thank you for your help!


## Post 2 by @pablo.strasser1 (2026-05-11T07:47:32.140Z)

This is just a guess, but I think they may have disabled this because of the big security flaw on linux.


## Post 3 by @Adrien.Albert (2026-05-11T09:49:02.216Z)

Hello,
Due to a critical security vulnerability, SSH access has been temporarily disabled on all compute nodes.
To mitigate the risk while minimizing disruption, we have:
- Drained all nodes and performed rolling reboots
- Ensured that running jobs were not interrupted, although nodes are currently inaccessible
SSH access will remain disabled until the vulnerability is fully patched.
Thank you for your understanding.
To follow issue:
[2026] Current issues on HPC Cluster[[2026] Current issues on HPC Cluster](https://hpc-community.unige.ch/t/2026-current-issues-on-hpc-cluster/4185/19) HPC ChangeLog[HPC ChangeLog](https://hpc-community.unige.ch/c/hpc-announce/hpc-changelog/9)
> Description
Cluster: All 
Dear users, 
We have been notified (again) of a critical security vulnerability. As a result, we must immediately patch our entire infrastructure, including compute and login nodes. 
The three login nodes will be rebooted shortly, and the cluster will be put into drain mode. 
The cluster will be resumed as soon as the nodes have been patched. 
Update:
For security reasons, we have disabled SSH access to compute nodes until they have all been patched. 
We will re-enable …


## Post 4 by @Max.Briel (2026-05-18T14:29:32.298Z)

Thank you for explaining this.
Is there any update on re-enabling access?
I, currently, cannot run interactive sessions for vscode/jupyterlab to which I connect from my local VSCode instance.


## Post 5 by @Vilius.Cepaitis (2026-05-19T08:44:23.715Z)

Hello,
I would also like to ask whether there is a rough timeline for when SSH access could be restored?
While there are workarounds, ssh access is crucial for rapid prototyping developments.
Best regards,
Vilius


## Post 6 by @Adrien.Albert (2026-05-19T09:06:58.855Z)

Hello everyone,
All compute nodes have now been patched, and SSH access has been fully restored.
You should once again be able to connect to the compute nodes via SSH without issue.
We apologize for the inconvenience caused and appreciate your patience during this operation.
Best regards,


## Post 7 by @Albert.Buchard (2026-06-03T10:09:20.848Z)

Do you guys have access to Mythos ? :star_struck:


## Post 8 by @Gael.Rossignol (2026-06-04T07:43:57.708Z)

Dear Albert,
Could you please clarify what you are referring to with “Mythos”? I am not familiar with it.
Best regards,
