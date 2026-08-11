# Running positron "kernel" on a worker node

- Source: https://hpc-community.unige.ch/t/4091

- Created: 2025-09-17T10:19:47.175Z

- Posts: 4

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Julien.Prados (2025-09-17T10:19:47.234Z)

Dear HPC community,
Positron is a data science IDE based on VScode (https://positron.posit.co[https://positron.posit.co](https://positron.posit.co)).
One of the interesting feature of the tool is its ability to have a kernel running remotely over a ssh connection. It would be nice for many of us if we could have it works on our HPC cluster, in a containerised environment…
This .ssh/config setup is kind of working for me, as it allows me to run a python kernel remotely on a worker node:
```
Host yggdrasil
  HostName login1.yggdrasil.hpc.unige.ch
  User prados

Host yggdrasil-worker
  User prados 
  ForwardAgent yes  
  StrictHostKeyChecking no  
  UserKnownHostsFile=/dev/null
  ProxyCommand ssh -T yggdrasil salloc bash -c 'nc $SLURM_NODELIST 22'
```
However, the R kernels are not working in my setup. I suspect there are issues when loading R modules dependencies (got a “ldd” library error). Furthermore, I would like to be able to run latest R version that are not available on the HPC for the moment.
The ideal would be to have the kernel running in this official R container on the worker node:
```
apptainer exec docker://docker.io/rocker/r-ver:4.5.1 bash
```
But cannot figure out how to change .ssh/config to have everything work together…
Can someone help me on this ?
Thanks a lot
Julien


## Post 2 by @Maroun.BouSleiman (2026-05-19T13:40:21.641Z)

Julien.Prados:
> Dear HPC community,
> Positron is a data science IDE based on VScode (https://positron.posit.co[https://positron.posit.co](https://positron.posit.co)).
> One of the interesting feature of the tool is its ability to have a kernel running remotely over a ssh connection. It would be nice for many of us if we could have it works on our HPC cluster, in a containerised environment…
> This .ssh/config setup is kind of working for me, as it allows me to run a python kernel remotely on a worker node:
```
Host yggdrasil
  HostName login1.yggdrasil.hpc.unige.ch
  User prados

Host yggdrasil-worker
  User prados 
  ForwardAgent yes  
  StrictHostKeyChecking no  
  UserKnownHostsFile=/dev/null
  ProxyCommand ssh -T yggdrasil salloc bash -c 'nc $SLURM_NODELIST 22'
```
> However, the R kernels are not working in my setup. I suspect there are issues when loading R modules dependencies (got a “ldd” library error). Furthermore, I would like to be able to run latest R version that are not available on the HPC for the moment.
> The ideal would be to have the kernel running in this official R container on the worker node:
```
apptainer exec docker://docker.io/rocker/r-ver:4.5.1 bash
```
> But cannot figure out how to change .ssh/config to have everything work together…
> Can someone help me on this ?
> Thanks a lot
Hi Julien,
I had a similar need for positron (and RStudio) on HPC and ended up trying many things. As my own solution is too specific to share, here is a partly AI-generated answer that is inspired by some of my findings. Containerising R is the right way to go, and Positron should itself be in the image. Something to do with Positron’s R kernel being dynamically linked against libR.so…
The pattern that may work for you is to run an sshd inside the container on the worker (rootless apptainer on a port different from 22, since it is already taken by host). You keep the SLURM allocation alive in a tmux session on the login node, and tunnel to that sshd from your workstation/latptop via `ssh -W` through the login node. Your .ssh/config then points Positron at the tunnelled endpoint. In this scenario, the R session survives disconnects because the allocation is held by the tmux.
Example apptainer def. You build once `apptainer build ~/images/r-ssh.sif r-ssh.def`
```
Bootstrap: docker
  From: rocker/r-ver:4.5.1
  
  %post
      apt-get update && apt-get install -y --no-install-recommends openssh-server
      rm -rf /var/lib/apt/lists/*
      
  %startscript
      exec /usr/sbin/sshd -D -e -f "$HOME/sshd_in_container.conf" -o "Port=$KERNEL_SSH_PORT"
```
Bind extra paths (e.g. /srv/beegfs/…) by adding --bind to the apptainer instance start line if your
data lives outside $HOME.
One-time setup on the cluster (host key + per-user sshd config)
```
ssh-keygen -t ed25519 -f ~/.ssh/container_host_key -N ''
  
cat > ~/sshd_in_container.conf <<EOF
HostKey /home/$USER/.ssh/container_host_key
AuthorizedKeysFile /home/$USER/.ssh/authorized_keys
PasswordAuthentication no
PubkeyAuthentication yes
UsePAM no
PidFile none
StrictModes no
EOF
```
Launcher script, e.g. `~/bin/start-rkernel.sh` (run inside a `tmux new -s rkernel` session on the login node)
```
#!/usr/bin/env bash
  set -euo pipefail
  PORT=$(python3 -c 'import socket;s=socket.socket();s.bind(("",0));print(s.getsockname()[1])')
  echo "$PORT" > ~/.rkernel.port
exec salloc -n1 -c4 --time=08:00:00 --partition=shared-cpu bash -c '
  echo "$SLURM_NODELIST" > ~/.rkernel.host
  export APPTAINERENV_KERNEL_SSH_PORT='"$PORT"'
  apptainer instance start --bind "$HOME" "$HOME/images/r-ssh.sif" rkernel
  apptainer exec instance://rkernel tail -f /dev/null
'
```
Detach with Ctrl-b d. The allocation, container, and any R session inside it survive until you `tmux￼kill-session -t rkernel` or time limit.
`~/.ssh/config` on your laptop
```
Host yggdrasil-r
   User pradosForwardAgent yes
   StrictHostKeyChecking no
   UserKnownHostsFile /dev/null
   ProxyCommand bash -c ‘H=$(ssh yggdrasil cat .rkernel.host); P=$(ssh yggdrasil cat .rkernel.port);exec ssh -W “$H:$P” yggdrasil’
```
Point Positron Remote-SSH at yggdrasil-r.
Note that this leads to three separate ssh calls. You can do something about it by adding ControlMaster, ControlPath, ControlPersist to the main yggdrasil entry in `.ssh/config`.
I have not tested this specific solution, as my setup and needs were different, but I hope this is useful at the higher level.
Cheers,
Maroun


## Post 3 by @Maroun.BouSleiman (2026-05-19T13:42:49.873Z)

by the way, I thought this was a recent question but just realized it is from last year. Wondering what you ended up doing…


## Post 4 by @Gael.Rossignol (2026-06-01T08:41:04.879Z)

Hello,
I think now the issue has gone because you can launch R using openondemand application.
Best regards,
