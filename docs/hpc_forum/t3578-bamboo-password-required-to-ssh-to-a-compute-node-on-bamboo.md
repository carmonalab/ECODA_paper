# Bamboo: Password required to ssh to a compute node on bamboo

- Source: https://hpc-community.unige.ch/t/3578

- Created: 2024-07-30T11:58:21.015Z

- Posts: 2

- Category: 1

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Malte.Algren (2024-07-30T11:58:21.052Z)

Hi HPC,
I am trying to use the new cluster bamboo. But im having an issue with password requirement when sshing directly into a compute node (specifically an issue in vs code)
Without any issues i can do:
```
ssh algren@login1.bamboo.hpc.unige.ch
or
ssh bm
```
and afterwards
```
salloc -n1 --partition=shared-gpu, --time=7:00:00 --mem=15G --gres=gpu:1
```
To get (as an example ) gpu001.
But if i want to ssh directly into gpu001 with `ssh diff~bmgp`, it asks me for a password.
Is it suppose to be asking for password when it already has my public key? How do i avoid it asking?
```
Host diff~*
  RemoteCommand apptainer shell -B /srv,/home --nv /home/users/a/algren/singularity_images/diffusion-torch-2.2.0.sif

Host bm
  HostName login1.bamboo.hpc.unige.ch
  User algren
  RequestTTY yes

Host bmgpu diff~bmgpu
  HostName gpu001
  User algren
  ProxyJump bm
    
```


## Post 2 by @Adrien.Albert (2024-07-30T23:41:08.542Z)

Hi @Malte.Algren[@Malte.Algren](https://hpc-community.unige.ch/u/malte.algren)
Could you try this:
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
Be carefull the ssh-copy-id does not work you will need to use the -f option.
And then it should work.
Tested and approuved:
```
[alberta@localhost .ssh]$ ssh -F config_bamboo cpu002 
alberta@cpu002's password: 
Killed by signal 2.

[alberta@localhost .ssh]$ ssh-copy-id login1.bamboo.hpc.unige.ch -f
/bin/ssh-copy-id: INFO: Source of key(s) to be installed: "/home/alberta/.ssh/id_rsa.pub"

Number of key(s) added: 1

Now try logging into the machine, with:   "ssh 'login1.bamboo.hpc.unige.ch'"
and check to make sure that only the key(s) you wanted were added.

[alberta@localhost .ssh]$ ssh -F config_bamboo cpu002 
Last login: Wed Jul 31 01:39:11 2024 from login1.bamboo
Installed: Fri May 31 14:59:01 CEST 2024
(bamboo)-[alberta@cpu002 ~]$
```
