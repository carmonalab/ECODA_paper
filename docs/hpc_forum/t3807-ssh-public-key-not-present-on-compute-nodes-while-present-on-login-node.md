# SSH public key not present on compute nodes, while present on login node

- Source: https://hpc-community.unige.ch/t/3807

- Created: 2025-02-03T13:38:59.634Z

- Tags: baobab

- Posts: 7

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Christian.Scheulen (2025-02-03T13:38:59.683Z)

## Primary informations
Username: scheulen
Cluster: Baobab
## Description
SSH key used to connect successfully to login node fails for proxy jump to specific CPU node reserved via slurm.
I just reserved a CPU node (cpu330) for remote software development via VSCode. Upon trying to proxy jump to the compute node (either within VSCode or on console) the ssh key authentication fails and a password prompt for `$USER@cpu330`, which does not accept the usual UniGe password. On the initial session on the node opened after successful slurm allocation, no authorised key was found using the `/usr/bin/sss_ssh_authorisedkey` script as required for authentication, while my public key is found with the same script on the login node.
Therefore, it seems that the key is not accessible on the CPU node. I am doubtful that this is expected behaviour since I was able to access nodes via keys in the past and since I am not aware that there is an alternative authorisation method such as Kerberos tokens to get to the compute nodes from the login node.
For completeness, an MWE of my `.ssh/config` is reproduced below:
```
Host baobab
     HostName                    login1.baobab.hpc.unige.ch
     User                        scheulen
     Compression                 yes
     ForwardX11                  yes
     ForwardAgent                yes
     PubkeyAuthentication        yes
     IdentityFile                ~/.ssh/id_ed25519

Host cpu* gpu*
     HostName                    %h
     User                        scheulen
     UserKnownHostsFile          ~/.ssh/known_hosts.unige
     ProxyJump                   baobab
     Compression                 yes
     ForwardAgent                yes
     PubkeyAuthentication        yes
     IdentityFile                ~/.ssh/id_ed25519
     
```
## Steps to Reproduce
- Log on to baobab via `ssh baobab`
- Reserve CPU node using `salloc -n1 --partition=private-dpnc-cpu,shared-cpu --time=8:00:00 --mem=10G -c2`
- Log on to CPU node via `ssh cpuXXX` (XXX being the specific node number)
- Cry a bit
- Try out `/usr/bin/sss_ssh_authorisedkey $USER` on both the login node and the CPU node (using the console opened after successful negotiation of the slurm job)
## Expected Result
Connect directly to CPU node with ssh key, show ssh key uploaded via https://my-account.unige.ch[https://my-account.unige.ch](https://my-account.unige.ch) on both the login and the compute node.
## Actual Result
Upon trying to connect to the CPU node via proxy jump, the ssh key was not accepted (running ssh verbosely shows the key is supplied for authentication and also does work for the login node, but fails for the CPU node). Similarly, `/usr/bin/sss_ssh_authorizedkeys $USER` shows the key on the login node but not on the CPU node.


## Post 2 by @Christian.Scheulen (2025-02-03T13:49:20.283Z)

I seem to have just solved it by adding the key to `.ssh/authorised_keys` instead of relying on the uploaded key in my account. In case the key uploaded to the UniGe account should be accessible by the compute nodes as well, there might however still be something wrong.


## Post 3 by @Adrien.Albert (2025-02-05T21:38:38.200Z)

Hi @Christian.Scheulen[@Christian.Scheulen](https://hpc-community.unige.ch/u/christian.scheulen)
Is this post is related to your issue ?
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


## Post 4 by @Christian.Scheulen (2025-02-17T12:54:09.849Z)

Hi @Adrien.Albert[@Adrien.Albert](https://hpc-community.unige.ch/u/Adrien.Albert),
Indeed, that post was the guideline I used initially to troubleshoot my issues. Initially, I was under the (presumably wrong) impression that the public key associated to my UniGe account should also be available on the worker nodes.
Problematically, when trying to copy the ssh-key from my computer to the login node via `copy-ssh-id`, the process failed as the key was already available on there (I am assuming due to the fact that I had already uploaded the same public key to my UniGe account). Manually adding the public key to `.ssh/authorized_keys` did however resolve the issues for me.
In any case, since copying the key manually I have not had problems using proxy-jumps to the worker nodes via the login node.
Best,
Chris


## Post 5 by @Adrien.Albert (2025-02-17T13:02:01.050Z)

To make it work:
- To login on cluster we use  the sshPublicKey registered in Unige Account. You need to have  the same ssh key in Unige Account and your local machine to connect on login (entry point)
- Then in cluster, we use the traditionnally way. You need to have the sshPublicKey from `login1.baobab:~/.ssh/id_rsa.pub` in your `login1.baobab:~/.ssh/authorized_keys` (or something else if you have configured your ssh_config)
EDIT:
- Then on cluster, your ~/.ssh/ must be empty.
Exemple:
```
# screen1
(yggdrasil)-[alberta@login1 .ssh]$ ls
(yggdrasil)-[alberta@login1 .ssh]$ 
(yggdrasil)-[alberta@login1 .ssh]$ srun sleep 600
...
```
```
# screen2
(yggdrasil)-[alberta@login1 ~]$ sac
          JobID    JobName    Account      User        NodeList   NTasks               Start                 End      State 
--------------- ---------- ---------- --------- --------------- -------- ------------------- ------------------- ---------- 
       38555419      sleep      burgi   alberta          cpu001          2025-03-05T10:11:28             Unknown    RUNNING 
(yggdrasil)-[alberta@login1 ~]$ ssh cpu001
Last login: Wed Mar  5 10:01:02 2025 from login1.yggdrasil
Installed: Thu Jan 30 09:26:01 CET 2025
(yggdrasil)-[alberta@cpu001 ~]$
```


## Post 6 by @Adrien.Albert (2025-03-05T09:14:02.497Z)

@Christian.Scheulen[@Christian.Scheulen](https://hpc-community.unige.ch/u/christian.scheulen)
I’ve updated my previous answer. I made a mistake :pray:


## Post 7 by @Christian.Scheulen (2025-06-03T13:55:27.444Z)

Hi @Adrien.Albert[@Adrien.Albert](https://hpc-community.unige.ch/u/Adrien.Albert),
Apologies for only coming back to this now, I should really activate the forum digest to not miss out on messages. Thanks a lot for the additional information on keeping `$HOME/.ssh` empty.
Since I put some aliases to get access to some CERN machines, this one unfortunately does not work for me. However, keeping the ssh-keys in `.ssh/authorized_keys` is working without issues so far.
Best,
Chris
