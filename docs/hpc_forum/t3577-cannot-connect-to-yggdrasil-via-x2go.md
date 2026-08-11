# Cannot connect to yggdrasil via X2Go

- Source: https://hpc-community.unige.ch/t/3577

- Created: 2024-07-29T22:54:32.728Z

- Tags: yggdrasil

- Posts: 14

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @William.Ceva (2024-07-29T22:54:32.772Z)

Hello,
I cannot connect to yggdrasil using X2Go.  Whenever I try to connect, the X2Go status changes to “running”, the window showing the yggdrasil screen opens briefly (but remains black, nothing loads in the window), and then closes again, with X2Go’s status changing to “Finished”, before I am prompted to log in again.
I have confirmed this issue is specifically with my X2Go connection to yggdrasil, since I can successfully connect to yggdrasil using ssh in terminal, and I can connect to other HPC clusters (for example, Lesta) using X2Go.  The issue is therefore neither solely related to my yggdrasil log in, nor solely due to my X2Go installation.
Any help would be greatly appreciated,
Will Ceva


## Post 2 by @Adrien.Albert (2024-07-30T07:23:43.317Z)

Hi @William.Ceva[@William.Ceva](https://hpc-community.unige.ch/u/william.ceva)
I have commented out the conda initializer in your `.bashrc`, could you try again?
On the other hand, does the ssh connection work fine?


## Post 3 by @William.Ceva (2024-07-30T13:38:23.563Z)

Hello,
Unfortunately, the issue still persists.  SSH still works though.
To confirm: I do not need to be in Switzerland to connect, correct?  (I am currently abroad).
I can connect to baobab using X2Go without issues, so I assume being abroad is not the issue.
Just in case this is relevant: until now, I do not believe I had tried to connect to yggdrasil since the maintenance ~2 weeks ago.
I am interested in trying the solution mentioned in this thread: X2go session closing immediately after connection - #3 by Pierre.Kuenzli[X2go session closing immediately after connection - #3 by Pierre.Kuenzli](https://hpc-community.unige.ch/t/x2go-session-closing-immediately-after-connection/886/3)
Can you confirm that it is safe for me to delete my ~/.x2go/ directory on yggdrasil via ssh?
Thank you,
Will Ceva


## Post 4 by @William.Ceva (2024-07-30T20:22:34.674Z)

I went ahead and tried the solution included in this link: X2go session closing immediately after connection - #3 by Pierre.Kuenzli[X2go session closing immediately after connection - #3 by Pierre.Kuenzli](https://hpc-community.unige.ch/t/x2go-session-closing-immediately-after-connection/886/3)
It did not solve the issue.


## Post 5 by @Adrien.Albert (2024-07-30T23:12:20.716Z)

Hi @William.Ceva[@William.Ceva](https://hpc-community.unige.ch/u/william.ceva)
According to the .xsession-x2go logs and the .ICEauthorithy, I suspect you are using ICEVM instead XFCE.
```
(yggdrasil)-[root@login1 ceva3]$  ll -ltra /home/users/c/ceva3
-rw-------   1 ceva3 hpc_users    1278 Jul 30 22:16 .xsession-x2go-login1.yggdrasil-errors.old
-rw-------   1 ceva3 hpc_users    1278 Jul 30 22:22 .xsession-x2go-login1.yggdrasil-errors
-rw-------   1 ceva3 hpc_users   14410 Jul 30 22:22 .ICEauthority
drwx------   3 ceva3 hpc_users       5 Jul 30 22:22 .gnupg
-rw-r--r--   1 ceva3 hpc_users     909 Jul 30 22:36 .bashrc
-rw-------   1 ceva3 hpc_users   12493 Jul 30 22:46 .bash_history
drwx------  18 ceva3 hpc_users      37 Jul 31 01:04 .
```
Could you give your x2go configuration for the yggdrasil session ?
for me:
image
image1379×959 52.4 KB
[image1379×959 52.4 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/2b8627d09c4221fee8fed258dc247fa9a76492d0.png)


## Post 6 by @William.Ceva (2024-07-31T02:11:18.505Z)

Hi,
I have always had it set to XFCE.
I changed it to ICEVM, then tried connecting, (which showed nothing at all), then changed back to XFCE and tried connecting again, and got the same issue as before (window opens briefly, then closes, then “Finished”, etc).
Screenshot from 2024-07-31 04-07-47
Screenshot from 2024-07-31 04-07-47807×744 74.3 KB
[Screenshot from 2024-07-31 04-07-47807×744 74.3 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/7b31707125ec384e9be6982a2145574a1c710b0c.png)


## Post 7 by @Adrien.Albert (2024-07-31T09:17:03.125Z)

Hi @William.Ceva[@William.Ceva](https://hpc-community.unige.ch/u/william.ceva)
I made some change, Could you try again please ?


## Post 8 by @William.Ceva (2024-07-31T12:42:47.264Z)

Hello,
I just tried again, same result as before but slower.
Is there any other way I can access Yggdrasil with GUI?  This would be 3 days of not being able to work due to this issue.
Sorry for not replying earlier in the day, but I am in a different time zone (6hrs behind Geneva).


## Post 9 by @William.Ceva (2024-07-31T13:03:08.051Z)

I also made a new session on x2go for connecting to yggdrasil (in case the issue was due to some specific session setting) but I still get the same result.
Also no change if I use a different Wi-Fi network  - I do not think it is an issue on my end, since I can connect to baobab and LESTA with x2go without issue.


## Post 10 by @William.Ceva (2024-07-31T13:26:20.696Z)

In case it helps, here are the full details from x2go during my last connection attempt:
```
NXPROXY - Version 3.5.99.26

Copyright (c) 2001, 2011 NoMachine (http://www.nomachine.com)
Copyright (c) 2008-2014 Oleksandr Shneyder <o.shneyder@phoca-gmbh.de>
Copyright (c) 2014-2016 Ulrich Sibiller <uli42@gmx.de>
Copyright (c) 2014-2016 Mihai Moldovan <ionic@ionic.de>
Copyright (c) 2011-2016 Mike Gabriel <mike.gabriel@das-netzwerkteam.de>
Copyright (c) 2015-2016 Qindel Group (http://www.qindel.com)

NXCOMP, NX protocol compression and NX extensions to this software
are copyright of the aforementioned persons and companies.

Redistribution and use of the present software is allowed according
to terms specified in the file LICENSE.nxcomp which comes in the
source distribution.

All rights reserved.

NOTE: This software has received contributions from various other
contributors, only the core maintainers and supporters are listed as
copyright holders. Please contact us, if you feel you should be listed
as copyright holder, as well.

NX protocol compression is derived from DXPC project.

Copyright (c) 1995,1996 Brian Pane
Copyright (c) 1996,1997 Zachary Vonler and Brian Pane
Copyright (c) 1999 Kevin Vigor and Brian Pane
Copyright (c) 2000,2003 Gian Filippo Pinzari and Brian Pane

All rights reserved.

See https://github.com/ArcticaProject/nx-libs for more information.

Info: Proxy running in server mode with pid '10981'.
Session: Starting session at 'Wed Jul 31 15:23:14 2024'.
Info: Using errors file '/home/will/.x2go/S-ceva3-54-1722432191_stDXFCE_dp24/sessions'.
Info: Using stats file '/home/will/.x2go/S-54/stats'.
Loop: WARNING! Overriding auxiliary X11 port with new value '1'.
Warning: Overriding auxiliary X11 port with new value '1'.
Info: Using abstract X11 socket in kernel namespace for accessing DISPLAY=:0.
Info: Connecting to remote host 'localhost:32601'.
Info: Connected to remote proxy on FD#5.
Info: Connection with remote proxy completed.
Loop: WARNING! Unrecognized session type 'unix-kde-depth_24'. Assuming agent session.
Warning: Unrecognized session type 'unix-kde-depth_24'. Assuming agent session.
Info: Using ADSL link parameters 1408/24/1/0.
Info: Using cache parameters 4/4096KB/8192KB/8192KB.
Info: Using pack method 'adaptive-9' with session 'unix-kde-depth_24'.
Info: Using ZLIB data compression 1/1/32.
Info: Using ZLIB stream compression 4/4.
Info: No suitable cache file found.
Info: Forwarding X11 connections to display ':0'.
Info: Forwarding auxiliary X11 connections to display ':0'.
Session: Session started at 'Wed Jul 31 15:23:30 2024'.
Info: Established X server connection.
Info: Using shared memory parameters 0/0K.
```


## Post 11 by @Adrien.Albert (2024-07-31T14:10:42.769Z)

Hi @William.Ceva[@William.Ceva](https://hpc-community.unige.ch/u/william.ceva)
I sent you a zoom invitation by email, can you join ?


## Post 12 by @William.Ceva (2024-07-31T18:02:14.127Z)

For those who may encounter this issue in the future, the issue appeared to be caused by the setting I added for the `PS1` parameter in my `.bashrc` file (which customizes the appearance of bash shells).
I had the following line written in my `.bashrc` file:
`PS1="\u@\h:\w\$ "`
By commenting out this line in my `.bashrc` file, I was able to connect yggdrasil via x2go without issues.


## Post 13 by @William.Ceva (2024-08-16T09:57:17.104Z)

I should also add: even after applying the above solution, I would still sometimes end up having the same original issue reappear randomly.
The solution is to use pure ssh connection to move your .bashrc file out of your ~/ directory, into some temporary directory, and then try x2go again.  This fixes the original issue, in the 2 random instances for me where it reappeared again.
I was then able to move my .bashrc back into ~/ and the x2go connection would continue to work.


## Post 14 by @Adrien.Albert (2024-08-16T13:41:08.274Z)

Hi @William.Ceva[@William.Ceva](https://hpc-community.unige.ch/u/william.ceva)
I’ve just published the new FAQ and included the step for troubleshooting x2go :
https://doc.eresearch.unige.ch/hpc/faq#x2go-desktop[https://doc.eresearch.unige.ch/hpc/faq#x2go-desktop](https://doc.eresearch.unige.ch/hpc/faq#x2go-desktop)
