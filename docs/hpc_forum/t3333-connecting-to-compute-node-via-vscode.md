# Connecting to compute node via VSCode?

- Source: https://hpc-community.unige.ch/t/3333

- Created: 2024-02-22T10:27:21.194Z

- Posts: 4

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Imahn.Shekhzadeh1 (2024-02-22T10:27:21.252Z)

Hi,
In one of your e-mails, you encouraged us to use VSCode not on a login node, but on a computing node. I requrested a compute node via `salloc`, typed `hostname` into the compute node to find out the hostname (in my case, this resulted in `cpu026.baobab`), but I’m not quite sure now how to proceed.
In VSCode, I press `Ctrl + Shift + P`, then `Remote-SSH: Connect to host...`, but then I’m not quite sure what to enter. Something like
```
ssh shekhza2@cpu026.baobab
```
doesn’t work (“hostname not resolved”), so I guess the hostname would somehow have to be combined with `login2.baobab.hpc.unige.ch`?
Best,
Imahn


## Post 2 by @Adrien.Albert (2024-02-26T15:15:39.740Z)

Dear @Imahn.Shekhzadeh1[@Imahn.Shekhzadeh1](https://hpc-community.unige.ch/u/imahn.shekhzadeh1)
- You can use OpenOndemand to request easily a compute node and connect on VScodeServer.
Or
- you can follow these procedure to configure :
- hpc:access_the_hpc_clusters [eResearch Doc][hpc:access_the_hpc_clusters [eResearch Doc]](https://doc.eresearch.unige.ch/hpc/access_the_hpc_clusters?s%5B%5D=proxyjump#ssh_tunnel_and_socks_proxy)
More info to configure with the sshPublicKey
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
And the procedure from VscodeServer to connect on host
code.visualstudio.com[code.visualstudio.com](https://code.visualstudio.com/docs/remote/ssh)
### Developing on Remote Machines using SSH and Visual Studio Code[Developing on Remote Machines using SSH and Visual Studio Code](https://code.visualstudio.com/docs/remote/ssh)
Developing on Remote Machines or VMs using Visual Studio Code Remote Development and SSH


## Post 3 by @Imahn.Shekhzadeh1 (2024-02-29T10:12:01.287Z)

Hi @Adrien.Albert[@Adrien.Albert](https://hpc-community.unige.ch/u/Adrien.Albert),
I just tried OpenOnDemand, but unfortunately, I’m getting the following error when logging in to https://ondemand.baobab.hpc.unige.ch[https://ondemand.baobab.hpc.unige.ch](https://ondemand.baobab.hpc.unige.ch):
```
Error -- can't find user for shekhzad
Run 'nginx_stage --help' to see a full list of available command line options.
```
This is weird, as my username is `shekhza2`, not `shekhzad` :thinking:


## Post 4 by @Yann.Sagon (2024-03-04T15:00:57.317Z)

Hi  @Imahn.Shekhzadeh1[@Imahn.Shekhzadeh1](https://hpc-community.unige.ch/u/imahn.shekhzadeh1)
you have two account at UNIGE, one as student and one as employee. On  Baobab you can only one have one account. Should we migrate your account shekhza2 to shekhzad?
Best
Yann
