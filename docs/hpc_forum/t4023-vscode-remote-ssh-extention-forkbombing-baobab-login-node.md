# VSCode remote-ssh Extention "Forkbombing" Baobab login node

- Source: https://hpc-community.unige.ch/t/4023

- Created: 2025-07-31T15:40:15.950Z

- Tags: baobab

- Posts: 3

- Category: 13

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Sara.Ricossa (2025-07-31T15:40:15.998Z)

Hi. A coworker and I are having this weird issue with VSCode since about 5pm today. The connection lagged for about 5 minutes, and when it came back everything was fine, except when we try to establish a connection using ssh with the `remote-ssh`[remote-ssh](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-ssh) extension, a large number of processes are spawned on the login node, making everything sluggish until we kill them.
This is the output of `ps -elf | grep ricossa` before trying to log in with vscode:
```
4 S root       82220    3356  0  80   0 -  5103 -      17:15 ?        00:00:00 sshd: ricossa [priv]
4 S ricossa    82225       1  0  80   0 -  5527 ep_pol 17:15 ?        00:00:00 /usr/lib/systemd/systemd --user
5 S ricossa    82227   82225  0  80   0 -  6532 -      17:15 ?        00:00:00 (sd-pam)
5 S ricossa    82288   82220  0  80   0 -  5103 -      17:15 ?        00:00:00 sshd: ricossa@pts/12
0 S ricossa    82289   82288  0  80   0 -  2750 do_wai 17:15 pts/12   00:00:00 -bash
0 R ricossa    84329   82289  0  80   0 -  2537 -      17:17 pts/12   00:00:00 ps -elf
0 S ricossa    84330   82289  0  80   0 -  1601 pipe_r 17:17 pts/12   00:00:00 grep --color=auto ricossa
```
This is the output while trying to log in with vscode:
```
-bash: fork: retry: Resource temporarily unavailable
-bash: fork: retry: Resource temporarily unavailable
-bash: fork: retry: Resource temporarily unavailable
-bash: fork: retry: Resource temporarily unavailable
-bash: fork: Resource temporarily unavailable
```
This is the output after shutting down vscode:
```
4 S root       82220    3356  0  80   0 -  5103 -      17:15 ?        00:00:00 sshd: ricossa [priv]
4 S ricossa    82225       1  0  80   0 -  5527 ep_pol 17:15 ?        00:00:00 /usr/lib/systemd/systemd --user
5 S ricossa    82227   82225  0  80   0 -  6532 -      17:15 ?        00:00:00 (sd-pam)
5 S ricossa    82288   82220  0  80   0 -  5103 -      17:15 ?        00:00:00 sshd: ricossa@pts/12
0 S ricossa    82289   82288  0  80   0 -  2750 do_wai 17:15 pts/12   00:00:00 -bash
0 S ricossa    88493       1  0  80   0 - 32148 futex_ 17:18 ?        00:00:00 /home/users/r/ricossa/.vscode-server/code-488a1f239235055e34e673291fb8d8c810886f81 command-shell --cli-data-dir /home/users/r/ricossa/.vscode-server/cli --parent-process-id 88472 --on-host=127.0.0.1 --on-port
0 R ricossa   108935   82289  0  80   0 -  2537 -      17:36 pts/12   00:00:00 ps -elf
0 S ricossa   108936   82289  0  80   0 -  1601 pipe_r 17:36 pts/12   00:00:00 grep --color=auto ricossa
```
And typing `pstree -p -s 88493` we get,
```
systemd(1)───code-488a1f2392(88493)─┬─{code-488a1f2392}(88500)
                                    ├─{code-488a1f2392}(88501)
                                    ├─{code-488a1f2392}(88502)
                                    ├─{code-488a1f2392}(88503)
                                    ├─{code-488a1f2392}(88504)
                                    ├─{code-488a1f2392}(88505)
                                    ├─{code-488a1f2392}(88506)
                                    ├─{code-488a1f2392}(88507)
                                    ├─{code-488a1f2392}(88508)
                                    ├─{code-488a1f2392}(88509)
                                    ├─{code-488a1f2392}(88510)
                                    ├─{code-488a1f2392}(88511)
                                    ├─{code-488a1f2392}(88512)
                                    ├─{code-488a1f2392}(88513)
                                    ├─{code-488a1f2392}(88514)
                                    ├─{code-488a1f2392}(88515)
                                    ├─{code-488a1f2392}(88516)
                                    ├─{code-488a1f2392}(88517)
                                    ├─{code-488a1f2392}(88518)
                                    ├─{code-488a1f2392}(88519)
                                    ├─{code-488a1f2392}(88520)
                                    ├─{code-488a1f2392}(88521)
                                    ├─{code-488a1f2392}(88522)
                                    ├─{code-488a1f2392}(88523)
                                    ├─{code-488a1f2392}(88524)
                                    ├─{code-488a1f2392}(88525)
                                    ├─{code-488a1f2392}(88526)
                                    ├─{code-488a1f2392}(88527)
                                    ├─{code-488a1f2392}(88528)
                                    ├─{code-488a1f2392}(88529)
                                    ├─{code-488a1f2392}(88530)
                                    ├─{code-488a1f2392}(88531)
                                    ├─{code-488a1f2392}(88532)
                                    ├─{code-488a1f2392}(88533)
                                    ├─{code-488a1f2392}(88534)
                                    ├─{code-488a1f2392}(88535)
                                    ├─{code-488a1f2392}(88536)
                                    ├─{code-488a1f2392}(88537)
                                    ├─{code-488a1f2392}(88538)
                                    ├─{code-488a1f2392}(88539)
                                    ├─{code-488a1f2392}(88540)
                                    ├─{code-488a1f2392}(88541)
                                    ├─{code-488a1f2392}(88542)
                                    ├─{code-488a1f2392}(88543)
                                    ├─{code-488a1f2392}(88544)
                                    ├─{code-488a1f2392}(88545)
                                    ├─{code-488a1f2392}(88546)
                                    ├─{code-488a1f2392}(88547)
                                    └─{code-488a1f2392}(88548)
```
This process makes the shell very sluggish, with constant messages of `fork: retry: Resource temporarily unavailable` until it is killed.
On the client side, I get these log messages:
```
[18:01:49.309] Log Level: 2
[18:01:49.324] VS Code version: 1.102.3
[18:01:49.324] Remote-SSH version: remote-ssh@0.120.0
[18:01:49.324] linux x64
[18:01:49.325] SSH Resolver called for "ssh-remote+baobab", attempt 1
[18:01:49.327] remote.SSH.useLocalServer = true
[18:01:49.327] remote.SSH.useExecServer = true
[18:01:49.327] remote.SSH.bindHost = {}
[18:01:49.327] remote.SSH.path = undefined
[18:01:49.327] remote.SSH.configFile = undefined
[18:01:49.327] remote.SSH.useFlock = true
[18:01:49.327] remote.SSH.lockfilesInTmp = false
[18:01:49.328] remote.SSH.localServerDownload = auto
[18:01:49.328] remote.SSH.remoteServerListenOnSocket = false
[18:01:49.328] remote.SSH.defaultExtensions = []
[18:01:49.328] remote.SSH.defaultExtensionsIfInstalledLocally = []
[18:01:49.328] remote.SSH.loglevel = 2
[18:01:49.328] remote.SSH.enableDynamicForwarding = true
[18:01:49.329] remote.SSH.enableRemoteCommand = false
[18:01:49.329] remote.SSH.serverPickPortsFromRange = {}
[18:01:49.330] remote.SSH.serverInstallPath = {}
[18:01:49.330] remote.SSH.permitPtyAllocation = false
[18:01:49.330] remote.SSH.preferredLocalPortRange = undefined
[18:01:49.330] remote.SSH.useCurlAndWgetConfigurationFiles = false
[18:01:49.330] remote.SSH.experimental.chat = false
[18:01:49.330] remote.SSH.experimental.enhancedSessionLogs = false
[18:01:49.330] remote.SSH.httpProxy = {}
[18:01:49.330] remote.SSH.httpsProxy = {}
[18:01:49.338] SSH Resolver called for host: baobab
[18:01:49.338] Setting up SSH remote "baobab"
[18:01:49.342] Acquiring local install lock: /tmp/vscode-remote-ssh-e1897406-install.lock
[18:01:49.343] Looking for existing server data file at /home/sbr/.config/Code/User/globalStorage/ms-vscode-remote.remote-ssh/vscode-ssh-host-e1897406-488a1f239235055e34e673291fb8d8c810886f81-0.120.0-es/data.json
[18:01:49.343] No existing data file
[18:01:49.343] Using commit id "488a1f239235055e34e673291fb8d8c810886f81" and quality "stable" for server
[18:01:49.343] Extensions to install: 
[18:01:49.347] Install and start server if needed
[18:01:49.351] PATH: /home/sbr/.local/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/opt/android-sdk/platform-tools:/opt/cuda/bin:/opt/cuda/nsight_compute:/opt/cuda/nsight_systems/bin:/usr/lib/jvm/default/bin:/usr/bin/site_perl:/usr/bin/vendor_perl:/usr/bin/core_perl:/usr/lib/rustup/bin
[18:01:49.351] Checking ssh with "ssh -V"
[18:01:49.356] > OpenSSH_10.0p2, OpenSSL 3.5.1 1 Jul 2025

[18:01:49.358] askpass server listening on /run/user/1000/vscode-ssh-askpass-d05cd3982f43cf60963a9752b82662dd9c2cef1b.sock
[18:01:49.359] Spawning local server with {"serverId":1,"ipcHandlePath":"/run/user/1000/vscode-ssh-askpass-a45199f3bc7dcebace2d69a1bd7571c0c3939db0.sock","sshCommand":"ssh","sshArgs":["-v","-T","-D","38235","-o","ConnectTimeout=15","baobab"],"serverDataFolderName":".vscode-server","dataFilePath":"/home/sbr/.config/Code/User/globalStorage/ms-vscode-remote.remote-ssh/vscode-ssh-host-e1897406-488a1f239235055e34e673291fb8d8c810886f81-0.120.0-es/data.json"}
[18:01:49.359] Local server env: {"SSH_AUTH_SOCK":"/run/user/1000/gcr/ssh","SHELL":"/bin/bash","DISPLAY":":0","ELECTRON_RUN_AS_NODE":"1","SSH_ASKPASS":"/home/sbr/.vscode/extensions/ms-vscode-remote.remote-ssh-0.120.0/out/local-server/askpass.sh","VSCODE_SSH_ASKPASS_NODE":"/opt/visual-studio-code/code","VSCODE_SSH_ASKPASS_EXTRA_ARGS":"","VSCODE_SSH_ASKPASS_MAIN":"/home/sbr/.vscode/extensions/ms-vscode-remote.remote-ssh-0.120.0/out/askpass-main.js","VSCODE_SSH_ASKPASS_HANDLE":"/run/user/1000/vscode-ssh-askpass-d05cd3982f43cf60963a9752b82662dd9c2cef1b.sock"}
[18:01:49.361] Spawned 78124
[18:01:49.361] Using connect timeout of 17 seconds
[18:01:49.424] > local-server-1> Running ssh connection command: ssh -v -T -D 38235 -o ConnectTimeout=15 baobab
[18:01:49.426] > local-server-1> Spawned ssh, pid=78132
[18:01:49.448] stderr> debug1: Server host key: ssh-rsa SHA256:tKqp4nljL+EGVKl8T0VF2nS36DkHVFMpLxQOPg/gKvg
[18:01:49.531] stderr> Authenticated to login1.baobab.hpc.unige.ch ([129.194.9.190]:22) using "publickey".
[18:01:50.100] > ready: c1b2bffd4128
[18:01:50.103] > Linux 5.14.0-503.40.1.el9_5.x86_64 #1 SMP PREEMPT_DYNAMIC Wed Apr 30 17:38:54 UTC 2025
[18:01:50.104] Platform: linux
[18:01:50.105] > /bin/bash
[18:01:50.105] Parent Shell: bash
[18:01:50.105] Parent Shell pid: 78124
[18:01:50.106] Waiting for subshell to start
[18:01:50.109] > 141980
[18:01:50.109] stdout -> '141980'
[18:01:50.109] sub-process detected
[18:01:50.117] > c1b2bffd4128: running
> Script executing under PID: 141980
[18:01:50.133] > Found existing installation at /home/users/r/ricossa/.vscode-server...
> Starting VS Code CLI...
[18:01:50.136] > Removing old logfile at /home/users/r/ricossa/.vscode-server/.cli.488a1f239235055e34e673291fb8d8c810886f81.log
[18:01:50.142] > Spawned remote CLI: 141998
[18:01:50.146] > Waiting for server log...
[18:01:50.178] stderr> main: fork: retry: Resource temporarily unavailable
[18:01:51.180] stderr> main: fork: retry: Resource temporarily unavailable
[18:01:53.179] stderr> main: fork: retry: Resource temporarily unavailable
[18:01:57.181] stderr> main: fork: retry: Resource temporarily unavailable
[18:02:05.181] stderr> main: fork: Resource temporarily unavailable
[18:02:05.182] > Waiting for server log...
[18:02:05.214] stderr> main: fork: retry: Resource temporarily unavailable
[18:02:06.215] stderr> main: fork: retry: Resource temporarily unavailable
.......
```
I hope someone can help me understand why this is happening. Nothing changed on my laptop. No updates. no nothing.
Thank you very much!
*Edited to include more detailed logs
*This probably happened after the reboot[reboot](https://hpc-community.unige.ch//hpc-community.unige.ch/t/2025-current-issues-on-hpc-cluster/3788/17).


## Post 2 by @Sara.Ricossa (2025-08-01T23:45:21.752Z)

It looks like the problem has been fixed maybe in some roundabout way. I’m still gonna leave this open in case it happens again.


## Post 3 by @Adrien.Albert (2025-08-03T17:15:29.437Z)

Hi  @Sara.Ricossa[@Sara.Ricossa](https://hpc-community.unige.ch/u/sara.ricossa)
this is a duplicate issue:
Cannot connect to Yggdrasil after the maintenance[Cannot connect to Yggdrasil after the maintenance](https://hpc-community.unige.ch/t/cannot-connect-to-yggdrasil-after-the-maintenance/4019/11) HPC issues[HPC issues](https://hpc-community.unige.ch/c/hpc-support/hpc-issues/14)
> Hello All, 
Yggdrasil
@Lena.Parc[@Lena.Parc](https://hpc-community.unige.ch/u/lena.parc) : no computing tasks on the login node, no compiling, building, or anything like that. wink 
I wanted to restrict you all to just basic commands like vim, emacs, nano, ls, rsync, and Slurm on the login node, but my kind colleagues keep stopping me from being too harsh,  so consider yourself lucky! grinning_face_with_smiling_eyes (joke… or not I keep the mystery) 
@Max.Briel[@Max.Briel](https://hpc-community.unige.ch/u/max.briel) You have vscode-server processes running on the login node and you source Conda in you…


## Post 4 by @Adrien.Albert (2025-08-03T17:15:36.246Z)
