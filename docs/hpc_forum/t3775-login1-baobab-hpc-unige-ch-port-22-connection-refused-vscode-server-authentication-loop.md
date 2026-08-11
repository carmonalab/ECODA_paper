# Login1.baobab.hpc.unige.ch port 22: Connection refused, vscode-server authentication loop

- Source: https://hpc-community.unige.ch/t/3775

- Created: 2024-12-27T11:00:27.649Z

- Tags: baobab

- Posts: 5

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Arthur.Freeman (2024-12-27T11:00:27.693Z)

## Primary informations
Username: freemana
Cluster: baobab (ex: Yggdrasil)
## Description
Hello, this morning when connecting to login1.baobab.unige.ch via vscode remote ssh extension, the application was loop asking me my password. I was however able to connect via a seperate terminal. After looking up the vscode server logs, I followed instructions here :
stackoverflow.com[stackoverflow.com](https://stackoverflow.com/questions/55979701/how-do-i-resolve-failed-to-parse-remote-port-from-server)
Sean Zamora
#### How do I resolve "Failed to parse remote port from server"?[How do I resolve "Failed to parse remote port from server"?](https://stackoverflow.com/questions/55979701/how-do-i-resolve-failed-to-parse-remote-port-from-server)
vscode-remote
asked by
  
  
    Sean Zamora
  [Sean Zamora](https://stackoverflow.com/users/11450767/sean-zamora)
  on 05:19AM - 04 May 19 UTC[05:19AM - 04 May 19 UTC](https://stackoverflow.com/questions/55979701/how-do-i-resolve-failed-to-parse-remote-port-from-server)
That recommended removing the folder `~/.vscode-server` on the remote server.
I did that, and now it doesn’t even prompt a password anymore, it just refuses the connection :
ssh: connect to host login1.baobab.hpc.unige.ch port 22: Connection refused
And I can’t login via terminal anymore either, same error. I suspect I’ve been banned.
Regards,
Arthur


## Post 2 by @Arthur.Freeman (2024-12-27T11:04:45.084Z)

Probably related, now baobab is claiming my password is wrong when I’m copy pasting from a password manager and the same password works on yggdrasil,
image
image799×740 56 KB
[image799×740 56 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/be858793bb20a086b19c0a11040d411713316737.png)
What is going on ?


## Post 3 by @Arthur.Freeman (2024-12-27T12:32:48.747Z)

Ok, now able to reconnect via CLI, however getting Remote I/O error when reading files,
image
image828×578 50 KB
[image828×578 50 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/f6e0d83e7a91d61e2ef528664c44c8d37e793797.png)


## Post 4 by @Arthur.Freeman (2024-12-27T12:38:34.498Z)

The full error logs from vscode are the following :
```
[13:35:44.556] Log Level: 2
[13:35:44.570] SSH Resolver called for "ssh-remote+login1.baobab.hpc.unige.ch", attempt 1
[13:35:44.571] "remote.SSH.useLocalServer": true
[13:35:44.571] "remote.SSH.useExecServer": true
[13:35:44.571] "remote.SSH.path": undefined
[13:35:44.571] "remote.SSH.configFile": undefined
[13:35:44.571] "remote.SSH.useFlock": true
[13:35:44.571] "remote.SSH.lockfilesInTmp": false
[13:35:44.572] "remote.SSH.localServerDownload": auto
[13:35:44.572] "remote.SSH.remoteServerListenOnSocket": false
[13:35:44.572] "remote.SSH.showLoginTerminal": false
[13:35:44.572] "remote.SSH.defaultExtensions": []
[13:35:44.572] "remote.SSH.loglevel": 2
[13:35:44.572] "remote.SSH.enableDynamicForwarding": true
[13:35:44.572] "remote.SSH.enableRemoteCommand": false
[13:35:44.572] "remote.SSH.serverPickPortsFromRange": {}
[13:35:44.572] "remote.SSH.serverInstallPath": {}
[13:35:44.572] "remote.SSH.permitPtyAllocation": false
[13:35:44.572] "remote.SSH.preferredLocalPortRange: undefined
[13:35:44.572] "remote.SSH.useCurlAndWgetConfigurationFiles: false
[13:35:44.580] VS Code version: 1.95.3
[13:35:44.580] Remote-SSH version: remote-ssh@0.115.1
[13:35:44.580] linux x64
[13:35:44.582] SSH Resolver called for host: login1.baobab.hpc.unige.ch
[13:35:44.582] Setting up SSH remote "login1.baobab.hpc.unige.ch"
[13:35:44.587] Acquiring local install lock: /tmp/vscode-remote-ssh-e2b59fc2-install.lock
[13:35:44.588] Looking for existing server data file at /home/gordon/.config/Code/User/globalStorage/ms-vscode-remote.remote-ssh/vscode-ssh-host-e2b59fc2-f1a4fb101478ce6ec82fe9627c43efbf9e98c813-0.115.1-es/data.json
[13:35:44.588] No existing data file
[13:35:44.589] Using commit id "f1a4fb101478ce6ec82fe9627c43efbf9e98c813" and quality "stable" for server
[13:35:44.590] Script variables:
 {
  "InstallExitCode.24": "AlreadyInProgress",
  "InstallExitCode.25": "ServerDownloadFailed",
  "InstallExitCode.26": "NoDownloaderAvailable",
  "InstallExitCode.27": "UnsupportedArch",
  "InstallExitCode.28": "StatusCheckFailed",
  "InstallExitCode.29": "NeedInsidersArch",
  "InstallExitCode.30": "NoDownloaderAvailableForStatusCheck",
  "InstallExitCode.31": "ServerTransferFailed",
  "InstallExitCode.32": "ServerFailedToStart",
  "InstallExitCode.33": "NeedInsidersWindows",
  "InstallExitCode.34": "CreateInstallDirFailed",
  "InstallExitCode.35": "UnsupportedPlatform",
  "InstallExitCode.36": "ServerTerminatedCVE20201416",
  "InstallExitCode.37": "UnpackFailed",
  "InstallExitCode.38": "ChangeDirFailed",
  "InstallExitCode.AlreadyInProgress": "24",
  "InstallExitCode.ServerDownloadFailed": "25",
  "InstallExitCode.NoDownloaderAvailable": "26",
  "InstallExitCode.NoDownloaderAvailableForStatusCheck": "30",
  "InstallExitCode.UnsupportedArch": "27",
  "InstallExitCode.StatusCheckFailed": "28",
  "InstallExitCode.NeedInsidersArch": "29",
  "InstallExitCode.ServerTransferFailed": "31",
  "InstallExitCode.ServerFailedToStart": "32",
  "InstallExitCode.NeedInsidersWindows": "33",
  "InstallExitCode.CreateInstallDirFailed": "34",
  "InstallExitCode.UnsupportedPlatform": "35",
  "InstallExitCode.ServerTerminatedCVE20201416": "36",
  "InstallExitCode.UnpackFailed": "37",
  "InstallExitCode.ChangeDirFailed": "38",
  "InstallUnpackCode.Success": "success",
  "InstallUnpackCode.Error": "error",
  "InstallUnpackCode.MissingFiles": "missingFiles",
  "uuid": "c079a0eac13c",
  "startMarker": "c079a0eac13c: running",
  "commitId": "f1a4fb101478ce6ec82fe9627c43efbf9e98c813",
  "quality": "stable",
  "token": "1aa111aa-a11a-1a1a-11a1-aaa1a11111aa",
  "vscodeAgentFolder": "$HOME/.vscode-server",
  "allowClientDownload": "1",
  "forceClientDownload": "0",
  "cliNameInArchive": "code",
  "ignoreWgetConfigFlag": " --no-config ",
  "ignoreCurlConfigFlag": " --disable ",
  "wgetTriesSegment": "--tries=1",
  "listenArgs": "--on-host=127.0.0.1 --on-port",
  "getDownloadServerStartTrigger": "c079a0eac13c:trigger_server_download",
  "getDownloadServerEndTrigger": "c079a0eac13c:trigger_server_download_end",
  "getProgressDownloading": "c079a0eac13c%%1%%",
  "getProgressInstalling": "c079a0eac13c%%2%%"
}
[13:35:44.592] Install and start server if needed
[13:35:44.599] PATH: /home/gordon/.local/bin:/home/gordon/miniconda3/bin:/home/gordon/miniconda3/condabin:/home/gordon/miniconda3/bin:/home/gordon/.local/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin:/snap/bin
[13:35:44.599] Checking ssh with "ssh -V"
[13:35:44.610] > OpenSSH_9.6p1 Ubuntu-3ubuntu13.5, OpenSSL 3.0.13 30 Jan 2024

[13:35:44.613] askpass server listening on /run/user/1000/vscode-ssh-askpass-e128601ef0a46e1512fc2f7952f5c259fabf248e.sock
[13:35:44.614] Spawning local server with {"serverId":1,"ipcHandlePath":"/run/user/1000/vscode-ssh-askpass-24f09e755c9d8d20f1aedeec756e8a377c261805.sock","sshCommand":"ssh","sshArgs":["-v","-T","-D","38531","-o","ConnectTimeout=15","login1.baobab.hpc.unige.ch"],"serverDataFolderName":".vscode-server","dataFilePath":"/home/gordon/.config/Code/User/globalStorage/ms-vscode-remote.remote-ssh/vscode-ssh-host-e2b59fc2-f1a4fb101478ce6ec82fe9627c43efbf9e98c813-0.115.1-es/data.json"}
[13:35:44.614] Local server env: {"SSH_AUTH_SOCK":"/run/user/1000/keyring/ssh","SHELL":"/bin/bash","DISPLAY":":0","ELECTRON_RUN_AS_NODE":"1","SSH_ASKPASS":"/home/gordon/.vscode/extensions/ms-vscode-remote.remote-ssh-0.115.1/out/local-server/askpass.sh","VSCODE_SSH_ASKPASS_NODE":"/snap/code/176/usr/share/code/code","VSCODE_SSH_ASKPASS_EXTRA_ARGS":"","VSCODE_SSH_ASKPASS_MAIN":"/home/gordon/.vscode/extensions/ms-vscode-remote.remote-ssh-0.115.1/out/askpass-main.js","VSCODE_SSH_ASKPASS_HANDLE":"/run/user/1000/vscode-ssh-askpass-e128601ef0a46e1512fc2f7952f5c259fabf248e.sock"}
[13:35:44.621] Spawned 17471
[13:35:44.622] Using connect timeout of 17 seconds
[13:35:44.707] > local-server-1> Running ssh connection command: ssh -v -T -D 38531 -o ConnectTimeout=15 login1.baobab.hpc.unige.ch
[13:35:44.711] > local-server-1> Spawned ssh, pid=17480
[13:35:44.715] stderr> OpenSSH_9.6p1 Ubuntu-3ubuntu13.5, OpenSSL 3.0.13 30 Jan 2024
[13:35:44.810] stderr> debug1: Server host key: ssh-rsa SHA256:tKqp4nljL+EGVKl8T0VF2nS36DkHVFMpLxQOPg/gKvg
[13:35:44.987] Got askpass request: {"request":" freemana@login1.baobab.hpc.unige.ch's password: "}
[13:35:44.987] Showing password prompt
[13:35:44.988] Listening for interwindow password on /run/user/1000/vscode-ssh-askpass-baaca9d08ff88e98a99b1fd56774346afcb4ce60.sock
[13:35:44.988] Writing password prompt to globalState
[13:35:46.920] Got password response
[13:35:46.921] Interactor gave response: ******************************
[13:35:46.921] Cleaning up other-window auth server
[13:35:46.923] Using connect timeout of 17 seconds
[13:35:47.186] stderr> Authenticated to login1.baobab.hpc.unige.ch ([129.194.9.190]:22) using "password".
[13:35:47.423] stderr> /home/users/f/freemana/.bash_profile: line 5: /home/users/f/freemana/.bashrc: Remote I/O error
[13:35:47.435] > ready: c079a0eac13c
[13:35:47.454] > Linux 4.18.0-553.el8_10.x86_64 #1 SMP Fri May 24 13:05:10 UTC 2024
[13:35:47.454] Platform: linux
[13:35:47.466] > /bin/bash
[13:35:47.467] Parent Shell: bash
[13:35:47.467] Parent Shell pid: 17471
[13:35:47.482] > 3206021
[13:35:47.482] Waiting for pid of spawned 'sh' subshell: '3206021'...
[13:35:47.506] > c079a0eac13c: running
> Script executing under PID: 3206021
[13:35:47.520] > Installing to /home/users/f/freemana/.vscode-server...
[13:35:47.525] > c079a0eac13c%%1%%
> Downloading with wget
> wget is from busybox: no
[13:35:47.533] stderr> Program 'wget' appears to support flag '--no-config'
[13:35:48.086] > Download complete
> c079a0eac13c%%2%%
> tar --version: tar (GNU tar) 1.30
> Copyright (C) 2017 Free Software Foundation, Inc.
> License GPLv3+: GNU GPL version 3 or later <https://gnu.org/licenses/gpl.html>.
> This is free software: you are free to change and redistribute it.
> There is NO WARRANTY, to the extent permitted by law.
> 
> Written by John Gilmore and Jay Fenlason.
[13:35:48.086] stderr> tar: code: Cannot write: Remote I/O error
[13:35:48.193] stderr> tar: code: Cannot utime: Remote I/O error
[13:35:48.193] stderr> tar: code: Cannot close: Remote I/O error
[13:35:48.194] stderr> tar: Exiting with failure status due to previous errors
[13:35:48.199] > ERROR: tar exited with a non-zero exit code: 2
> Trigger local server download
> c079a0eac13c:trigger_server_download
> artifact==cli-alpine-x64==
> destFolder==/home/users/f/freemana/.vscode-server==
> destFolder2==/vscode-cli-f1a4fb101478ce6ec82fe9627c43efbf9e98c813.tar.gz==
> c079a0eac13c:trigger_server_download_end
> Waiting for client to transfer server archive...
> Waiting for /home/users/f/freemana/.vscode-server/vscode-cli-f1a4fb101478ce6ec82fe9627c43efbf9e98c813.tar.gz.done and vscode-server.tar.gz to exist
>  
[13:35:48.200] Got request to download on client for {"artifact":"cli-alpine-x64","destPath":"/home/users/f/freemana/.vscode-server/vscode-cli-f1a4fb101478ce6ec82fe9627c43efbf9e98c813.tar.gz"}
[13:35:48.200] server download URL: https://update.code.visualstudio.com/commit:f1a4fb101478ce6ec82fe9627c43efbf9e98c813/cli-alpine-x64/stable
[13:35:48.201] Downloading VS Code server locally...
[13:35:48.313] Downloaded VS Code server to /tmp/a29053c9-2fd6-45ee-8c90-f953248d9f88
[13:35:48.314] Renamed VS Code server to /tmp/vscode_server_1735302948313/vscode-cli-f1a4fb101478ce6ec82fe9627c43efbf9e98c813.tar.gz
[13:35:48.314] Preparing to scp to host login1.baobab.hpc.unige.ch
[13:35:48.315] PATH: /home/gordon/.local/bin:/home/gordon/miniconda3/bin:/home/gordon/miniconda3/condabin:/home/gordon/miniconda3/bin:/home/gordon/.local/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin:/snap/bin
[13:35:48.315] Checking ssh with "ssh -V"
[13:35:48.324] > OpenSSH_9.6p1 Ubuntu-3ubuntu13.5, OpenSSL 3.0.13 30 Jan 2024

[13:35:48.325] Testing scp with "scp"
[13:35:48.331] scp exited with code: 1
[13:35:48.331] Got stderr from scp: usage: scp [-346ABCOpqRrsTv] [-c cipher] [-D sftp_server_path] [-F ssh_config]
           [-i identity_file] [-J destination] [-l limit] [-o ssh_option]
           [-P port] [-S program] [-X sftp_option] source ... target
[13:35:48.332] Copying file to remote with scp -o ConnectTimeout=15 'vscode-cli-f1a4fb101478ce6ec82fe9627c43efbf9e98c813.tar.gz' 'vscode-cli-f1a4fb101478ce6ec82fe9627c43efbf9e98c813.tar.gz.done' 'login1.baobab.hpc.unige.ch':'/home/users/f/freemana/.vscode-server'
[13:35:48.332] Using cwd: file:///tmp/vscode_server_1735302948313
[13:35:48.660] > freemana@login1.baobab.hpc.unige.ch's password: 
[13:35:48.660] Showing password prompt
[13:35:51.361] >  
[13:35:54.332] >  
[13:35:57.215] >  
[13:35:57.800] Got password response
[13:35:57.801] "Copy server to host" wrote data to terminal: "******************************"
[13:35:57.820] > 
[13:35:58.140] > vscode-cli-f1a4fb101478ce6ec82fe9627c43efbf9e   0%    0     0.0KB/s   --:-- ETA
[13:35:59.139] > vscode-cli-f1a4fb101478ce6ec82fe9627c43efbf9e  82% 7264KB   7.1MB/s   00:00 ETA
[13:35:59.343] > vscode-cli-f1a4fb101478ce6ec82fe9627c43efbf9e 100% 8835KB   7.2MB/s   00:01    
[13:35:59.397] > vscode-cli-f1a4fb101478ce6ec82fe9627c43efbf9e   0%    0     0.0KB/s   --:-- ETA
[13:35:59.409] > vscode-cli-f1a4fb101478ce6ec82fe9627c43efbf9e 100%    9     0.7KB/s   00:00    
[13:35:59.424] > scp: remote fsetstat: Failure
[13:35:59.692] "Copy server to host" terminal command done
[13:36:00.374] > Found flag and server on host
> c079a0eac13c%%2%%
> tar --version: tar (GNU tar) 1.30
> Copyright (C) 2017 Free Software Foundation, Inc.
> License GPLv3+: GNU GPL version 3 or later <https://gnu.org/licenses/gpl.html>.
> This is free software: you are free to change and redistribute it.
> There is NO WARRANTY, to the extent permitted by law.
> 
> Written by John Gilmore and Jay Fenlason.
[13:36:00.375] stderr> tar: code: Cannot write: Remote I/O error
[13:36:00.436] stderr> tar: code: Cannot utime: Remote I/O error
[13:36:00.436] stderr> tar: code: Cannot close: Remote I/O error
[13:36:00.436] stderr> tar: Exiting with failure status due to previous errors
[13:36:00.442] > ERROR: tar exited with a non-zero exit code: 2
> Already attempted local download, failing
> c079a0eac13c: start
> exitCode==37==
> listeningOn====
> osReleaseId==rocky==
> arch==x86_64==
> vscodeArch==x64==
> bitness==64==
> tmpDir==/tmp==
> platform==linux==
> unpackResult==error==
> didLocalDownload==1==
> downloadTime==447==
> installTime==223==
> serverStartTime====
> execServerToken==1aa111aa-a11a-1a1a-11a1-aaa1a11111aa==
[13:36:00.452] > platformDownloadPath==cli-alpine-x64==
> c079a0eac13c: end
[13:36:00.452] Received install output: 
exitCode==37==
listeningOn====
osReleaseId==rocky==
arch==x86_64==
vscodeArch==x64==
bitness==64==
tmpDir==/tmp==
platform==linux==
unpackResult==error==
didLocalDownload==1==
downloadTime==447==
installTime==223==
serverStartTime====
execServerToken==1aa111aa-a11a-1a1a-11a1-aaa1a11111aa==aaaaaaaaAaaaaaaaAaaa==aaa-aaaaaa-a11==

[13:36:00.453] Terminating local server
[13:36:00.455] Resolver error: Error: Failed to install the VS Code Server
	at v.ServerInstallError (/home/gordon/.vscode/extensions/ms-vscode-remote.remote-ssh-0.115.1/out/extension.js:2:493365)
	at h (/home/gordon/.vscode/extensions/ms-vscode-remote.remote-ssh-0.115.1/out/extension.js:2:486874)
	at t.handleInstallOutput (/home/gordon/.vscode/extensions/ms-vscode-remote.remote-ssh-0.115.1/out/extension.js:2:488844)
	at e (/home/gordon/.vscode/extensions/ms-vscode-remote.remote-ssh-0.115.1/out/extension.js:2:545236)
	at process.processTicksAndRejections (node:internal/process/task_queues:95:5)
	at async /home/gordon/.vscode/extensions/ms-vscode-remote.remote-ssh-0.115.1/out/extension.js:2:567561
	at async t.withShowDetailsEvent (/home/gordon/.vscode/extensions/ms-vscode-remote.remote-ssh-0.115.1/out/extension.js:2:571256)
	at async /home/gordon/.vscode/extensions/ms-vscode-remote.remote-ssh-0.115.1/out/extension.js:2:541941
	at async T (/home/gordon/.vscode/extensions/ms-vscode-remote.remote-ssh-0.115.1/out/extension.js:2:539992)
	at async t.resolveWithLocalServer (/home/gordon/.vscode/extensions/ms-vscode-remote.remote-ssh-0.115.1/out/extension.js:2:541481)
	at async P (/home/gordon/.vscode/extensions/ms-vscode-remote.remote-ssh-0.115.1/out/extension.js:2:564693)
	at async t.resolve (/home/gordon/.vscode/extensions/ms-vscode-remote.remote-ssh-0.115.1/out/extension.js:2:568667)
	at async /home/gordon/.vscode/extensions/ms-vscode-remote.remote-ssh-0.115.1/out/extension.js:2:839059
[13:36:00.459] ------

[13:36:00.462] Local server exit: 15
```
So it seems there are `Remote I/O errors` preventing vscode-server from installing.
```
...
[13:35:47.423] stderr> /home/users/f/freemana/.bash_profile: line 5: /home/users/f/freemana/.bashrc: Remote I/O error
...
[13:35:48.086] stderr> tar: code: Cannot write: Remote I/O error
[13:35:48.193] stderr> tar: code: Cannot utime: Remote I/O error
[13:35:48.193] stderr> tar: code: Cannot close: Remote I/O error
```
Is there perhaps some issue with the data partitions ?


## Post 5 by @Arthur.Freeman (2024-12-30T11:10:17.656Z)

Solved following I/O errors on Baobab home directory - #6 by Arthur.Freeman[I/O errors on Baobab home directory - #6 by Arthur.Freeman](https://hpc-community.unige.ch/t/i-o-errors-on-baobab-home-directory/3777/6) investigation
