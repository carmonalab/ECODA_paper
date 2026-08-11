# I/O errors on Baobab home directory

- Source: https://hpc-community.unige.ch/t/3777

- Created: 2024-12-28T16:27:12.680Z

- Tags: baobab

- Posts: 6

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Andrea.Serpolla (2024-12-28T16:27:12.729Z)

## Primary informations
Username: serpolla
Cluster: Baobab
## Description
There seems to be some problem with the home directory.
Many files are impossible to read due to I/O errors.


## Post 2 by @Arthur.Freeman (2024-12-29T09:44:15.009Z)

Same exact problem here, I created an issue already and called up Yann to no response, I can login but files I’ve created aren’t readable,
image
image795×568 89.8 KB
[image795×568 89.8 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/260c316ff52e037b089d26e88cfea3197210f31f.png)
Login1.baobab.hpc.unige.ch port 22: Connection refused, vscode-server authentication loop[Login1.baobab.hpc.unige.ch port 22: Connection refused, vscode-server authentication loop](https://hpc-community.unige.ch/t/login1-baobab-hpc-unige-ch-port-22-connection-refused-vscode-server-authentication-loop/3775/4) HPC issues[HPC issues](https://hpc-community.unige.ch/c/hpc-support/hpc-issues/14)
> The full error logs from vscode are the following : 
[13:35:44.556] Log Level: 2
[13:35:44.570] SSH Resolver called for "ssh-remote+login1.baobab.hpc.unige.ch", attempt 1
[13:35:44.571] "remote.SSH.useLocalServer": true
[13:35:44.571] "remote.SSH.useExecServer": true
[13:35:44.571] "remote.SSH.path": undefined
[13:35:44.571] "remote.SSH.configFile": undefined
[13:35:44.571] "remote.SSH.useFlock": true
[13:35:44.571] "remote.SSH.lockfilesInTmp": false
[13:35:44.572] "remote.SSH.localServerDownloa…
Remote I/O error likely means something is wrong with the filesystem, @Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon) au secours !


## Post 3 by @Cody.Cardenas (2024-12-29T14:11:05.531Z)

I have also seen this error on baobab.
```
$ ssh cardenac@login1.baobab.hpc.unige.ch
(cardenac@login1.baobab.hpc.unige.ch) Password: 
Last login: Sun Dec 29 14:51:54 2024 from app6.baobab
 ____              _           _
|  _ \            | |         | |
| |_) | __ _  ___ | |__   __ _| |__
|  _ < / _` |/ _ \| '_ \ / _` | '_ \
| |_) | (_| | (_) | |_) | (_| | |_) |
|____/ \__,_|\___/|_.__/ \__,_|_.__/
             _             _      __ 
            | |           (_)    /_ |
            | | ___   __ _ _ _ __ | |
            | |/ _ \ / _` | | '_ \| |
            | | (_) | (_| | | | | | |
            |_|\___/ \__, |_|_| |_|_|
                      __/ |          
                      |___/  

 Documentation: https://doc.eresearch.unige.ch/hpc/start
 Forum: https://hpc-community.unige.ch/
 OpenOndemand: https://openondemand.baobab.hpc.unige.ch/
 support: https://doc.eresearch.unige.ch/hpc/start#support_-_get_help

-bash: /home/users/c/cardenac/.bash_profile: Remote I/O error
(baobab)-[cardenac@login1 ~]$ 
```
But bamboo looks fine:
```
~$ ssh cardenac@login1.bamboo.hpc.unige.ch
(cardenac@login1.bamboo.hpc.unige.ch) Password: 
Last login: Sun Dec 29 13:50:50 2024 from 77-59-137-204.dclient.hispeed.ch
 ____                  _
|  _ \                | |
| |_) | __ _ _ __ ___ | |__   ___   ___
|  _ < / _` | '_ ` _ \| '_ \ / _ \ / _ \
| |_) | (_| | | | | | | |_) | (_) | (_) |
|____/ \__,_|_| |_| |_|_.__/ \___/ \___/
                 _             _      __ 
                | |           (_)    /_ |
                | | ___   __ _ _ _ __ | |
                | |/ _ \ / _` | | '_ \| |
                | | (_) | (_| | | | | | |
                |_|\___/ \__, |_|_| |_|_|
                          __/ |          
                         |___/  

 Documentation: https://doc.eresearch.unige.ch/hpc/start
 Forum: https://hpc-community.unige.ch/
 OpenOndemand: https://openondemand.baobab.hpc.unige.ch/
 support: https://doc.eresearch.unige.ch/hpc/start#support_-_get_help

(base) (bamboo)-[cardenac@login1 ~]$ 
```
I’m not sure if its connected, I was trying to run a quick debug script on the debug-cpus (PID 14190096 & 14190098) but the slurm file retuns an empty slurm-*.out file.
`seff` shows a `State: FAILED (exit code 15)` for both.


## Post 4 by @Adrien.Albert (2024-12-29T18:26:45.243Z)

Hi,
I’ll check and let you know.


## Post 5 by @Adrien.Albert (2024-12-29T18:53:29.820Z)

The problem has been solved, I will not investigate further due to the vacation week.
We’ll be back on January 3, 2025. :tada:


## Post 6 by @Arthur.Freeman (2024-12-30T11:09:46.260Z)

Thank you @Adrien.Albert[@Adrien.Albert](https://hpc-community.unige.ch/u/Adrien.Albert) !!
