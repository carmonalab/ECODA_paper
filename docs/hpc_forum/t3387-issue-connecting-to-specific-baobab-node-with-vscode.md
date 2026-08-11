# Issue connecting to specific baobab node with vscode

- Source: https://hpc-community.unige.ch/t/3387

- Created: 2024-03-21T13:02:05.722Z

- Posts: 4

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Stephen.Mulligan (2024-03-21T13:02:05.791Z)

## Primary informations
Username: $mulligas
Cluster: $baobab2
## Description
I am trying to connect my vscode session to a specific cpu nodes in order to be able to debug remotely. I ssh in the terminal to baobab, allocate the cpu node, and then update my .ssh config file to use the relevant node via a proxy jump. When connecting I am prompted for my ssh-key password as normal, but then afterwords I am prompted to log in to mulligas@cpuXXX , which I am not able to log into using my ssh-key password or computing account password.
## Steps to Reproduce
```
.ssh config contents
Host baobab
  HostName baobab2.hpc.unige.ch
  User mulligas
  RequestTTY yes

Host skycurtains~*
  RemoteCommand apptainer shell -B /srv,/home --nv /home/users/m/mulligas/temp_container/skycurtains.sif
  RequestTTY yes

Host baocurt skycurtains~baocurt
  HostName cpuXXX
  User mulligas
  ProxyJump baobab
  RequestTTY yes
```
salloc cpu on baobab. replace the cpuXXX with the allocated cpu.
open vs code remote explorer and connect to host baocurt. input ssky-key password.
then prompted to enter password for allocated cpu
## Expected Result
Connect to the remote host without being prompted to login to cpu
## Actual Result
Prompted to login to cpu, no password works.


## Post 2 by @Adrien.Albert (2024-03-26T17:53:09.097Z)

Hey @Stephen.Mulligan[@Stephen.Mulligan](https://hpc-community.unige.ch/u/stephen.mulligan)
I just follow the procedure from scratch and it’s working for me:
https://doc.eresearch.unige.ch/hpc/access_the_hpc_clusters?s[]=proxyjump#ssh_tunnel_and_socks_proxy[https://doc.eresearch.unige.ch/hpc/access_the_hpc_clusters?s[]=proxyjump#ssh_tunnel_and_socks_proxy](https://doc.eresearch.unige.ch/hpc/access_the_hpc_clusters?s%5B%5D=proxyjump#ssh_tunnel_and_socks_proxy)
Please also refer to the note at the end of the procedure:
ProxyJump ssh not working on Baobab[ProxyJump ssh not working on Baobab](https://hpc-community.unige.ch/t/proxyjump-ssh-not-working-on-baobab/3078/15) HPC Technical[HPC Technical](https://hpc-community.unige.ch/c/hpc-technical/5)
> Hi @Tomke.Schroeer[@Tomke.Schroeer](https://hpc-community.unige.ch/u/tomke.schroeer) 
I did the procedure again starting from the begening and it’s working. 

On your local machine, Save old ssh key and create a new one

mkdir ~/.ssh/old
mv ~/.ssh/*  ~/.ssh/old
ssh-keygen

On the cluster, make sure you have not id_rsa key (make a back up too) 

Copy the rsa.pub in https://my-account.unige.ch/main/home[https://my-account.unige.ch/main/home](https://my-account.unige.ch/main/home) end wait for 5-10 min the synchronisation with AD is done.

3 On your local machine configure the proxyjump: 
[alberta@localhost .ssh]$ cat ~/.ssh/config

host b…
Let me know if you still have the issue :slight_smile:


## Post 3 by @Stephen.Mulligan (2024-04-22T07:32:43.658Z)

Hi Adrien,
Thanks very much for the reply.
I was able to get this working by following the steps in the linked post, however now I am experiencing another issue.
“Failed to set up dynamic port forwarding connection over SSH to the VS Code Server.”
I am now unable to connect once again. Any help would be much appreciated.


## Post 4 by @Stephen.Mulligan (2024-04-22T07:34:51.861Z)

Killing the remote server fixed this issue.
