# Lost in Baobab (did Baobab crash?)

- Source: https://hpc-community.unige.ch/t/4039

- Created: 2025-08-11T15:30:04.762Z

- Tags: baobab

- Posts: 18

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Nicolas.Clairis1 (2025-08-11T15:30:04.805Z)

## Primary informations
Username: clairis
Cluster: Baobab
## Description
My last batches on Baobab just crashed with the “FAILED” error message based on the email I received. I therefore wanted to go to check the logs in FileZilla. My previous connection was not working so I tried to reconnect but instead of arriving in my folder as always, I ended up in an unknown location with folder names that I have never seen before (acanas, acanas-fiteo, etc.; cf screenshot attached) and no trace of my actual Baobab folders, neither scratch nor regular scripts… Then I wanted to check on Putty since my connection was supposed to be still valid, but when launching even a small command like `ls` I get the following error message: `ls: reading directory '.': Communication error on send`. Do you know what is going on and how long this will take to repair? Is Baobab down? Why does my folder not appear but I can still connect to some unknown location? I hope you can help with all these questions as I’m currently lost in Baobab
image
image162×433 3.04 KB
[image162×433 3.04 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/a55a7b29b9fbf6f03da8fb358e1a01c68e57db9d.png)
image
image169×425 3.4 KB
[image169×425 3.4 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/1b505155234dcefebba9e26ab006a62c2d38981a.png)


## Post 2 by @michael.divia (2025-08-11T15:34:42.407Z)

I think there is a problem with the /home directory :
```
(baobab)-[divia@login1 /]$  ls -lah
ls: cannot access 'home': Communication error on send
total 224K
dr-xr-xr-x.   27 root root 4.0K Jun 12 14:42 .
dr-xr-xr-x.   27 root root 4.0K Jun 12 14:42 ..
-rw-r--r--     1 root root    0 Jun  3 13:15 .autorelabel
drwxr-xr-x    16 root root 2.0K Aug 19  2024 acanas
drwxr-xr-x     5 root root 2.0K Feb 25  2021 acanas-fiteo
dr-xr-xr-x.    2 root root    6 Nov  3  2024 afs
lrwxrwxrwx.    1 root root    7 Nov  3  2024 bin -> usr/bin
dr-xr-xr-x.    6 root root 4.0K Jun  3 13:09 boot
drwxr-xr-x     2 root root    0 Aug 11 10:04 cvmfs
drwxr-xr-x     3 root root   27 Jun  3 13:51 datacube
drwxr-xr-x    21 root root 3.6K Jul 31 16:36 dev
drwxr-xr-x     2 root root    0 Aug 11 09:38 dpnc
drwxr-xr-x.  162 root root 8.0K Aug 11 09:12 etc
d??????????    ? ?    ?       ?            ? home
lrwxrwxrwx.    1 root root    7 Nov  3  2024 lib -> usr/lib
lrwxrwxrwx.    1 root root    9 Nov  3  2024 lib64 -> usr/lib64
drwxr-xr-x.    2 root root    6 Nov  3  2024 media
drwxr-xr-x     2 root root    0 Jul 31 16:36 misc
drwxr-xr-x.    2 root root    6 Nov  3  2024 mnt
drwxr-xr-x     2 root root    0 Jul 31 16:36 net
drwxr-xr-x.   12 root root  173 Jul 14 18:01 opt
dr-xr-xr-x  2232 root root    0 Jul 31 16:36 proc
dr-xr-x---.    5 root root 4.0K Aug  8 13:58 root
drwxr-xr-x    47 root root 1.5K Jul 31 16:42 run
lrwxrwxrwx.    1 root root    8 Nov  3  2024 sbin -> usr/sbin
drwxr-xr-x    10 root root  130 Jun 12 14:42 srv
dr-xr-xr-x    13 root root    0 Jul 31 16:36 sys
drwxrwxrwx  1881 root root 832K Aug 11 17:29 tmp
drwxr-xr-x     2 root root    0 Aug 11 06:37 unep
drwxr-xr-x     2 root root    0 Aug  2 12:57 unige
drwxr-xr-x.   12 root root  144 Jun  3 12:29 usr
drwxr-xr-x.   22 root root 4.0K Jun  3 13:15 var
```


## Post 3 by @Nicolas.Clairis1 (2025-08-11T15:40:15.590Z)

Don’t know exactly what this means in practice but indeed. I tried to disconnect and reconnect via Putty and I get a similar error message:
```
Could not chdir to home directory /home/users/c/clairis: Communication error on send
-bash: /home/users/c/clairis/.bash_profile: Communication error on send
```
Now `ls` works in Putty but it seems that I’m again in that mysterious place as for FileZilla:
image
I hope this gets fixed soon.


## Post 4 by @Jade.Awada (2025-08-11T15:51:31.881Z)

hello!!
I guess that there is a problem, I’m having the same issue when I try to access Baobab:
Could not chdir to home directory /home/users/a/awada: Communication error on send
/usr/bin/xauth: error in locking authority file /home/users/a/awada/.Xauthority-bash: /home/users/a/awada/.bash_profile: Communication error on send


## Post 5 by @Ariel.PerezMellor (2025-08-12T06:37:38.051Z)

Dear all,
I am experiencing the same issue.
Could not chdir to home directory /home/users/p/perezmel: Communication error on send
-bash: /home/users/p/perezmel/.bash_profile: Communication error on send
I will appreciate any further information in this respect. Thanks.


## Post 6 by @nicolas.clairis (2025-08-12T07:15:37.701Z)

The issue seems solved now (at least I can connect via FileZilla and Putty from my home pc)


## Post 7 by @Marco.Froelich (2025-08-12T07:16:41.880Z)

I have a different error message, but access to Baobab login has been unsuccessful for the same timeline (still on-going).
What I tried: ssh to login1.baobab.hpc.unige.ch
What happens: Instant error message appears: <Could not establish connection to “login1.baobab.hpc.unige.ch”: Failed to create the remote server’s install directory.>


## Post 8 by @nicolas.clairis (2025-08-12T07:21:22.189Z)

I have access both via Putty and FileZilla but weirdly the folder where the log files of my last batches should have been stored has completely vanished. If I try to recreate a folder with the same name, I get the following error message: `Erreur :	mkdir /home/users/c/clairis/scripts/logs/MIST: received failure with description 'Failure`. The folder nevertheless still appears on Putty, but if I try to get inside, I get the following error message
`-bash: cd: MIST: Communication error on send` I don’t really understand what is going on here. When typing `squeue` I see that many batches of other people are still ongoing so I guess that the cluster is not completely down


## Post 9 by @nicolas.clairis (2025-08-12T07:23:01.448Z)

When I try to submit a batch, the submission works fine but then it crashes immediately… I hope we get some info as to what is going on soon


## Post 10 by @michael.divia (2025-08-12T07:36:42.943Z)

Home directorys seem to be back up but now there is the exact same problem but for all folders inside our respecive home directorys.


## Post 11 by @michael.divia (2025-08-12T07:40:52.028Z)

[2025] Current issues on HPC Cluster[[2025] Current issues on HPC Cluster](https://hpc-community.unige.ch//hpc-community.unige.ch/t/2025-current-issues-on-hpc-cluster/3788/19) HPC ChangeLog[HPC ChangeLog](https://hpc-community.unige.ch/c/hpc-announce/hpc-changelog/9)
> Description
Cluster: Baobab 
Last night, we encountered an issue affecting the home storage system, caused by an unexpected stoppage of the meta services. As a result, some jobs running at the time, as well as those triggered during this period, may have failed. 
The services have now been fully restored, and all systems are operating normally. Our team is currently investigating the root cause to prevent similar occurrences in the future. 
We appreciate your understanding. 
– 
HPC Team 
Status …
They say it is resolved but it’s not fully working on my side.


## Post 12 by @Alexander.Froch (2025-08-12T07:50:38.291Z)

Same here (our whole office has the same issue). The login node still shows the issue when we connect. So looks like the home is still broken there


## Post 13 by @michael.divia (2025-08-12T08:08:24.277Z)

@nicolas.clairis[@nicolas.clairis](https://hpc-community.unige.ch/u/nicolas.clairis) @Alexander.Froch[@Alexander.Froch](https://hpc-community.unige.ch/u/alexander.froch) Everything is back and working on my side. Go have a look if everything works on your side now.


## Post 14 by @nicolas.clairis (2025-08-12T08:11:07.111Z)

thanks for the notice! Indeed same here, the folders that had vanished appear again and I can launch batches normally again as well.


## Post 16 by @Gael.Rossignol (2025-08-12T08:24:24.516Z)

Dear all,
Some files was always inacessible, so after checking all logs I fully reboot all storage home.
Could you please test if all is now solved?
Best regards,


## Post 17 by @michael.divia (2025-08-12T08:25:28.435Z)

Hello, on my side the problem just poped back for like 2 minutes and now I have access to all my files.


## Post 18 by @nicolas.clairis (2025-08-12T08:25:54.832Z)

I see that would explain why my batches just crashed again without any error message then :sweat_smile:
looks ok now on my side at least


## Post 19 by @Ariel.PerezMellor (2025-08-12T09:18:14.321Z)

It is already solved. Thanks!
