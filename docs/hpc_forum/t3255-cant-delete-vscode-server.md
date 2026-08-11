# Can't delete .vscode-server

- Source: https://hpc-community.unige.ch/t/3255

- Created: 2024-01-17T14:57:41.872Z

- Posts: 16

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Debajyoti.Sengupta (2024-01-17T14:57:41.950Z)

Hello,
I have been having some issues with connecting vscode to a compute node through SSH, so I decided to do a fresh install. I started by trying to remove the `.vscode-server` folder on my home directory. But it doesn’t seem to let me do that at the moment. I can confirm that no vscode processes were running for at the time.
```
(baobab)-[senguptd@login2 ~]$ rm -rf ~/.vscode-server/
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/node_modules/graceful-fs': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/node_modules/socks/build/client': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/node_modules/fd-slicer': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/node_modules/rc/lib': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/node_modules/node-pty/lib': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/node_modules/tas-client-umd': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/node_modules/jschardet/dist': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/node_modules/decompress-response': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/node_modules/util-deprecate': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/node_modules/fs-constants': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/node_modules/minimist': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/node_modules/vscode-textmate': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/node_modules/yazl': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/node_modules/semver/internal': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/node_modules/semver/functions': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/node_modules/semver/classes': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/node_modules/semver/ranges': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/node_modules/mkdirp/lib': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/node_modules/lru-cache': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/node_modules/node-vsce-sign/src': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/node_modules/pend': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/node_modules/node-abi': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/extensions/ms-vscode.js-debug/resources/dark': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/extensions/ms-vscode.js-debug/src': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/extensions/gulp': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/extensions/jake': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/extensions/node_modules/typescript/lib/cs': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/extensions/node_modules/typescript/lib/it': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/extensions/node_modules/typescript/lib/tr': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/extensions/node_modules/typescript/lib/ja': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/extensions/node_modules/typescript/lib/pt-br': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/extensions/node_modules/typescript/lib/es': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/extensions/node_modules/typescript/lib/zh-cn': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/extensions/node_modules/typescript/lib/zh-tw': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/extensions/ms-vscode.vscode-js-profile-table/out': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/extensions/ms-vscode.vscode-js-profile-table/resources': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/extensions/emmet/dist/node': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/extensions/ms-vscode.js-debug-companion/resources': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/extensions/merge-conflict/dist': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/extensions/html-language-features/icons': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/extensions/html-language-features/client/dist/node': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/extensions/git/resources/icons/dark': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/extensions/git/dist': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/extensions/php-language-features/dist': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/extensions/debug-auto-launch/dist': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/extensions/debug-auto-launch/.vscode': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/extensions/grunt/images': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/extensions/grunt/dist': Device or resource busy
```
Any clues as to what’s going wrong?


## Post 2 by @Adrien.Albert (2024-01-17T15:07:47.088Z)

Hi @Debajyoti.Sengupta[@Debajyoti.Sengupta](https://hpc-community.unige.ch/u/debajyoti.sengupta)
Debajyoti.Sengupta:
```
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/node_modules/tas-client-umd': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/node_modules/jschardet/dist': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/node_modules/decompress-response': Device or resource busy
rm: cannot remove '/home/users/s/senguptd/.vscode-server/bin/0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2/node_modules/util-deprecate': Device or resource busy
rm: cannot remove '/home/users/s/sengupt
```
Did you use VScodeServer with OpenOnDemand  ?


## Post 3 by @Debajyoti.Sengupta (2024-01-17T15:12:27.908Z)

Hi Adrien,
No I am currently using the traditional SSH connection. For some reason I cannot seem to access the link (https://ondemand.baobab.hpc.unige.ch/[https://ondemand.baobab.hpc.unige.ch/](https://ondemand.baobab.hpc.unige.ch/)) when I log in with my SWITCHedu-ID.
```
Unauthorized
This server could not verify that you are authorized to access the document requested. Either you supplied the wrong credentials (e.g., bad password), or your browser doesn't understand how to supply the credentials required.
```


## Post 4 by @Adrien.Albert (2024-01-17T15:23:16.878Z)

@Debajyoti.Sengupta[@Debajyoti.Sengupta](https://hpc-community.unige.ch/u/debajyoti.sengupta)
I give you access to OpenOnDemand, Welcome to club 33 ! :european_castle:
You need to close all your running vscode server process before delete files


## Post 5 by @Debajyoti.Sengupta (2024-01-17T15:30:08.099Z)

Hi Adrien,
Thanks. I still can’t seem to remove the folders in the directory. I am pretty sure there are no vscode-server processes running for me (nor is any vscode window open on my local machine).


## Post 6 by @Adrien.Albert (2024-01-17T15:40:58.701Z)

Did you used conda or apptainer ?


## Post 7 by @Debajyoti.Sengupta (2024-01-17T15:43:48.475Z)

Hi @Adrien.Albert[@Adrien.Albert](https://hpc-community.unige.ch/u/adrien.albert)
I currently shell in using apptainer, if that is what you are asking. Here’s my ssh config.
```
Host skycurtains~*
    RemoteCommand apptainer shell --nv -B /srv,/home /home/users/s/senguptd/UniGe/astro/skycurtains/container/curtains_gaia.sif && cd /home/users/s/senguptd/UniGe/astro/skycurtains/
    RequestTTY yes

Host baobab
    HostName baobab2.hpc.unige.ch
    User senguptd
    IdentityFile path/to/file
    RequestTTY yes

Host nodebaobab  skycurtains~nodebaobab
    HostName cpu224
    User senguptd
    ProxyJump baobab
    IdentityFile path/to/file
    RequestTTY yes
```
My original problem was VSCode kept trying to connect to the compute node but got stuck on `waiting for server log`


## Post 8 by @Adrien.Albert (2024-01-17T16:04:36.444Z)

could you please disconnect from the cluster ? I will check something


## Post 9 by @Debajyoti.Sengupta (2024-01-17T16:05:10.365Z)

Sure thing. All connections terminated.


## Post 10 by @Adrien.Albert (2024-01-17T16:22:00.967Z)

There is something somewere which is writting in file at 16:59, you need to find the origins.
```
(baobab)-[senguptd@login2 .vscode-server]$ ll -la
total 4
drwxr-xr-x  3 senguptd private_dpnc   5 Jan 17 17:04 .
drwx------ 41 senguptd private_dpnc  82 Jan 17 17:00 ..
-rw-r--r--  1 senguptd private_dpnc 529 Jan 17 16:59 .0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2.log
-rw-r--r--  1 senguptd private_dpnc   8 Jan 17 16:59 .0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2.pid
-rwx------  1 senguptd private_dpnc  37 Jan 17 16:59 .0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2.token
drwxr-xr-x  3 senguptd private_dpnc   1 Jan 17 09:40 bin
-rw-r--r--  1 senguptd private_dpnc 224 Jan 17 16:59 .cli.0ee08df0cf4527e40edc9aa28f4b5bd38bbff2b2.log
```
on cpu314 you something a bit questionable, why srun 9h with nothing other process running ?
```
senguptd 1921480 1921479  0 16:47 ?        00:00:00 srun sleep 9h
senguptd 1921494 1921488  0 16:47 ?        00:00:00 /usr/bin/sleep 9h
```


## Post 11 by @Debajyoti.Sengupta (2024-01-17T16:27:18.362Z)

Troubling. Those might be from trying to mount vscode on a compute node. Is there any way to force delete those?
As for CPU node, I usually request a compute node for 9 hours where I can Then connect vscode to. Currently which is idle on request.


## Post 12 by @Adrien.Albert (2024-01-17T16:44:19.292Z)

Well there is no process on cluster, I really do not have idea right know.
```
(baobab)-[root@admin1 senguptd]$ clush -bw @all "ps -edf |grep sengu |grep -v grep"
(baobab)-[root@admin1 senguptd]$
```
Well I change the name of the directory
`mv .vscode-server/ .vscode-server-trash`
Concerning your `srun sleep 9h`. 2024 has  just started, so new resolution !! :rofl:
Could you try VScodeServer through OpenOnDemand ?


## Post 13 by @Adrien.Albert (2024-01-17T16:46:10.459Z)

PS: I canceled your job, to try deletion.


## Post 14 by @Debajyoti.Sengupta (2024-01-18T08:19:37.272Z)

Thanks a lot, it seems to work now.
I will try code-server through OoD today and report.


## Post 15 by @John.Raine (2024-01-18T08:56:20.303Z)

Hi Adrien,
the problem with OoD is that you cannot run vscode within a different apptainer container providing access to the development environment.
Without that funcitonality it is preferable to use the standard remote-development options with vscode through the desktop client.
Cheers,
Johnny


## Post 16 by @Adrien.Albert (2024-01-18T09:51:01.920Z)

Hi @John.Raine[@John.Raine](https://hpc-community.unige.ch/u/john.raine)
I don’t understand your situation. Could you please try to detail your use case, maybe there is something we can do to improve the experience with OOD.


## Post 17 by @Adrien.Albert (2024-01-23T13:41:40.310Z)
