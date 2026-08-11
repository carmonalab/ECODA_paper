# Cannot connect to yggdrasil

- Source: https://hpc-community.unige.ch/t/3471

- Created: 2024-06-04T13:27:06.072Z

- Tags: yggdrasil

- Posts: 26

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @William.Ceva (2024-06-04T13:27:06.114Z)

## Primary informations
Username: ceva3
Cluster: Yggdrasil
## Description
I cannot connect to yggdrasil right now.  As of ~15 minutes ago, whenever I try to connect to yggdrasil, I get a “Timeout connecting” error.  I am using X2Go to try and connect, but I assume X2Go is not the issue here, since I have confirmed with other users (who connect using different methods) that yggdrasil appears to be unavailable.


## Post 2 by @Anton.Chudaykin (2024-06-04T13:57:25.332Z)

I experience similar issue. The server is extremely slow.


## Post 3 by @Adrien.Albert (2024-06-04T14:06:47.099Z)

Hi
Here the related issue. You should have access again :slight_smile:
[2024] Current issues on HPC Cluster[[2024] Current issues on HPC Cluster](https://hpc-community.unige.ch//hpc-community.unige.ch/t/2024-current-issues-on-hpc-cluster/3245/12) HPC ChangeLog[HPC ChangeLog](https://hpc-community.unige.ch/c/hpc-announce/hpc-changelog/9)
> Yggdrasil : Login node stuck
Dear users, 
The login node (login1.yggdrasil) is blocked due to a computation process. 
The server has been restarted and is available again. 
Thank you for your understanding. 
Best regards, 

Status :  Resolved white_check_mark
Start: 2024-06-04T13:40:00Z (UTC) 
End:2024-06-04T14:05:00Z (UTC)
Best Regards,


## Post 4 by @William.Ceva (2024-06-11T07:17:02.371Z)

The connection to yggdrasil appears to be down again, and I confirmed this with other users.


## Post 5 by @Lorenzo.Bini (2024-06-11T08:55:20.865Z)

Me as well, same kind of error


## Post 6 by @Giuseppe.Chindemi (2024-06-11T11:46:03.012Z)

Any news? I’m also experiencing the same issue.


## Post 7 by @Anton.Chudaykin (2024-06-11T11:46:07.203Z)

I can not establish connection to yggdrasil as well.


## Post 8 by @William.Ceva (2024-06-11T12:18:56.903Z)

I am able to connect now


## Post 9 by @Adrien.Albert (2024-06-11T12:23:55.087Z)

[2024] Current issues on HPC Cluster[[2024] Current issues on HPC Cluster](https://hpc-community.unige.ch//hpc-community.unige.ch/t/2024-current-issues-on-hpc-cluster/3245/13) HPC ChangeLog[HPC ChangeLog](https://hpc-community.unige.ch/c/hpc-announce/hpc-changelog/9)
> Yggdrasil : Login node stuck
Dear users, 
The login node (login1.yggdrasil) is blocked due to a computation process. 
The server has been restarted and is available again. 
Thank you for your understanding. 
Best regards, 

Status :  Resolved white_check_mark
Start: 2024-06-11T07:30:00Z (UTC) 
End:2024-06-11T12:15:00Z (UTC)


## Post 10 by @Giuseppe.Chindemi (2024-07-08T14:05:15.495Z)

I can’t connect to the login node anymore…


## Post 11 by @Jonathan.Mutal (2024-07-08T14:47:16.916Z)

Same here, I can ping the login node but I cannot connect.


## Post 12 by @Adrien.Albert (2024-07-08T15:17:06.715Z)

Hi @Jonathan.Mutal[@Jonathan.Mutal](https://hpc-community.unige.ch/u/jonathan.mutal)
It should work now. :slight_smile:


## Post 13 by @Elliot.Jiwani-Brown (2024-07-08T15:32:34.076Z)

Hi, I am also struggling to connect to the yggdrasil login node. Either with X2Go and ssh, both give me a timeout error and fail to login


## Post 14 by @Giuseppe.Chindemi (2024-07-08T15:33:56.129Z)

Mh, it seems still unavailable on my end.


## Post 15 by @Giuseppe.Chindemi (2024-07-11T21:03:31.277Z)

Is Yggdrasil unreachable again?


## Post 16 by @Jonathan.Mutal (2024-07-12T04:58:25.621Z)

Yes, It looks like it is unreachable. We can ping it, but when we are asked for the password, it does not log in.
Thanks in advance to the support team for working on the issue.


## Post 17 by @Adrien.Albert (2024-07-12T08:19:58.949Z)

Hi All,
The Linux authentication service was on failure and have been fixed. It should work again.


## Post 18 by @Arianna.Nigioni (2024-07-23T07:18:04.550Z)

The connection to Yggdrasil seems down again for both me and other users.


## Post 19 by @Adrien.Albert (2024-07-23T08:26:42.853Z)

Dear All,
A user process failed a core dump, which resulted in a cascade of errors, causing the connection node to freeze.
The user in question has been contacted and reminded of the impact on the community of running jobs on the login node. I also forwarded the three forum posts concerning this issue to the user.
The login node was restarted at approximately 9:55 AM and is now available again.
I apologize for any inconvenience caused.
Best regards


## Post 20 by @Matthias.Kruckow (2024-08-02T08:16:14.565Z)

My open connections to yggdrasil freeze and I can’t reconnect. Even a ping to the login node fails currently.
