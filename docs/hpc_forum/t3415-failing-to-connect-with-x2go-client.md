# Failing to connect with X2Go client

- Source: https://hpc-community.unige.ch/t/3415

- Created: 2024-04-15T13:20:24.443Z

- Posts: 3

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Pierre.Megevand (2024-04-15T13:20:24.500Z)

## Primary informations
Username: megevanp
Cluster: baobab
## Description
Dear all, I cannot connect to Baobab using the X2Go client any more. I can retrieve files using FileZilla, and launch jobs using Putty, but X2Go throws an error. Below is the error log from X2Go:
```
NXPROXY - Version 3.5.0

Copyright (C) 2001, 2010 NoMachine.
See http://www.nomachine.com/ for more information.

Info: Proxy running in client mode with pid '8700'.
Session: Starting session at 'Mon Apr 15 15:15:29 2024'.
Info: Connecting to remote host 'localhost:58556'.
Info: Connection to remote proxy 'localhost:58556' established.
Info: Connection with remote proxy completed.
Warning: Unrecognized session type 'unix-kde-depth_32'. Assuming agent session.
Warning: Failed to read data from the X auth command.
Warning: Generated a fake cookie for X authentication.
Info: Using WAN link parameters 768/24/1/0.
Info: Using cache parameters 4/4096KB/8192KB/8192KB.
Info: Using pack method '16m-jpeg-9' with session 'unix-kde-depth_32'.
Info: Using ZLIB data compression 1/1/32.
Info: Using ZLIB stream compression 1/1.
Info: Using cache file '/cygdrive/C/Users/megevanp/X2GO~1/cache-unix-kde-depth_32/S-FD8A465E10D163255E9C16DFA06E698A'.
Info: Forwarding X11 connections to display 'localhost:0'.
Session: Session started at 'Mon Apr 15 15:15:29 2024'.
Warning: Cookie mismatch in the X authentication data.
Session: Terminating session at 'Mon Apr 15 15:15:34 2024'.
Info: Your session was closed before reaching a usable state.
Info: This can be due to the local X server refusing access to the client.
Info: Please check authorization provided by the remote X application.
Session: Session terminated at 'Mon Apr 15 15:15:34 2024'.
Connection timeout, aborting
NXPROXY - Version 3.5.0

Copyright (C) 2001, 2010 NoMachine.
See http://www.nomachine.com/ for more information.

Info: Proxy running in client mode with pid '18144'.
Session: Starting session at 'Mon Apr 15 15:16:19 2024'.
Info: Connecting to remote host 'localhost:58556'.
Info: Connection to remote proxy 'localhost:58556' established.
Info: Connection with remote proxy completed.
Warning: Unrecognized session type 'unix-kde-depth_32'. Assuming agent session.
Warning: Failed to read data from the X auth command.
Warning: Generated a fake cookie for X authentication.
Info: Using WAN link parameters 768/24/1/0.
Info: Using cache parameters 4/4096KB/8192KB/8192KB.
Info: Using pack method '16m-jpeg-9' with session 'unix-kde-depth_32'.
Info: Using ZLIB data compression 1/1/32.
Info: Using ZLIB stream compression 1/1.
Info: Using cache file '/cygdrive/C/Users/megevanp/X2GO~1/cache-unix-kde-depth_32/S-FD8A465E10D163255E9C16DFA06E698A'.
Info: Forwarding X11 connections to display 'localhost:0'.
Session: Session started at 'Mon Apr 15 15:16:19 2024'.
Warning: Cookie mismatch in the X authentication data.
Session: Terminating session at 'Mon Apr 15 15:16:23 2024'.
Info: Your session was closed before reaching a usable state.
Info: This can be due to the local X server refusing access to the client.
Info: Please check authorization provided by the remote X application.
Session: Session terminated at 'Mon Apr 15 15:16:23 2024'.
Connection timeout, aborting
```
Thank you in advance for any advice,
Pierre


## Post 2 by @Adrien.Albert (2024-04-16T10:31:54.755Z)

Dear @Pierre.Megevand[@Pierre.Megevand](https://hpc-community.unige.ch/u/pierre.megevand)
Have you tried The Desktop session on https://ondemand.baobab.hpc.unige.ch/[https://ondemand.baobab.hpc.unige.ch/](https://ondemand.baobab.hpc.unige.ch/)?
Since OpenOnDemand offers an easy-to-use desktop, we’re recommending it against X2go.


## Post 3 by @Yann.Sagon (2024-06-14T12:23:55.970Z)

Dear @Pierre.Megevand[@Pierre.Megevand](https://hpc-community.unige.ch/u/pierre.megevand)
we were able to reproduce the issue to connect to x2go.We fixed the issue, please try again.
