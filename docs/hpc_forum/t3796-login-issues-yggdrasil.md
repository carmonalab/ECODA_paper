# Login Issues Yggdrasil

- Source: https://hpc-community.unige.ch/t/3796

- Created: 2025-01-23T13:26:50.436Z

- Tags: yggdrasil

- Posts: 20

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @felix.vecchi (2025-01-23T13:26:50.482Z)

## Primary informations
Username: fvecchi
Cluster: Yggdrasil
## Description
I cant login on Yggdrasil and neither can my colleagues. Is there something going on?
Thanks Felix


## Post 2 by @Siddharth.Bhatnagar (2025-01-23T13:59:51.512Z)

I confirm the same (at least since 12:45pm)


## Post 3 by @Siddharth.Bhatnagar (2025-01-23T15:24:38.750Z)

Update: I’m able to login. However, now I get a “Remote I/O error.”


## Post 4 by @daniel.forerosanchez (2025-01-23T15:39:24.225Z)

I confirm I can now log in too but still get Remote IO error as before the login issue.


## Post 5 by @Adrien.Albert (2025-01-23T15:49:13.061Z)

I did some witchcraft, do you still have problems?


## Post 6 by @Siddharth.Bhatnagar (2025-01-23T15:50:39.900Z)

The witchcraft helped, thank you! :wink:


## Post 7 by @daniel.forerosanchez (2025-01-23T15:56:58.795Z)

I am still getting IO error with `ls` command


## Post 8 by @Davide.Pietrobon (2025-01-23T16:22:28.400Z)

I am also experiencing the same Remote I/O error when attempting to access my directory on the Yggdrasil cluster. Specifically, the error occurs when I try to list the contents of my home directory (`/home/users/p/pietrobo`) using the `ls` command.
## Steps to reproduce
- Log in to the cluster:
```
ssh pietrobo@login1
```
- Navigate to my home directory:
```
cd /home/users/p/pietrobo
```
- Attempt to list files:
```
ls
```
- Error output:
```
ls: reading directory '.': Remote I/O error
```
Thank you for your help.


## Post 9 by @maciej.falkiewicz (2025-01-23T16:25:49.020Z)

+1 me too (adding more characters)


## Post 10 by @Adrien.Albert (2025-01-24T09:05:18.280Z)

The situation seems to have stabilized, but is anyone still experiencing the issue?


## Post 11 by @daniel.forerosanchez (2025-01-24T13:04:22.668Z)

Problem seems to be back as of 14h04.


## Post 12 by @Adrien.Albert (2025-01-24T13:15:37.077Z)

I am investigating  I keep you inform


## Post 13 by @Adrien.Albert (2025-01-24T13:26:30.462Z)

A good candidate for the root cause has been found and contacted. I’ll keep an eye on him and login1.
Let me know if you still have an issue.


## Post 14 by @daniel.forerosanchez (2025-01-24T13:30:48.339Z)

Still an issue but maybe it will take a bit longer?


## Post 15 by @Adrien.Albert (2025-01-24T13:53:44.971Z)

Okay guys,  I need volunteers @daniel.forerosanchez[@daniel.forerosanchez](https://hpc-community.unige.ch/u/daniel.forerosanchez) if you still have the issue please join the following meeting:
Zoom Video[Zoom Video](https://unige.zoom.us/j/63790873024)
### Join our Cloud HD Video Meeting[Join our Cloud HD Video Meeting](https://unige.zoom.us/j/63790873024)
Zoom is the leader in modern enterprise video communications, with an easy, reliable cloud platform for video and audio conferencing, chat, and webinars across mobile, desktop, and room systems. Zoom Rooms is the original software-based conference...


## Post 16 by @daniel.forerosanchez (2025-01-27T09:27:53.063Z)

Hi, just an update, Yggdrasil hanging for me when trying command.


## Post 17 by @Nicolas.Cuny (2025-01-27T09:54:55.410Z)

I also have troubles to connect or copy files from yggdrasil to my computer. It’s either filing or taking a lot of time.


## Post 18 by @Siddharth.Bhatnagar (2025-01-27T09:56:09.542Z)

Same problem for me, unable to login or extract files.


## Post 20 by @Jaime.RomanGarza (2025-01-27T13:36:31.181Z)

Hi, I can’t access Yggdrasil or download data from it, this is after today’s morning fix to the login


## Post 21 by @Yann.Sagon (2025-01-28T08:55:34.369Z)

As we are doing a cluster maintenance today, I’m closing this topic that will be obsolete. If you still have an issue after the maintenance, feel free to create a new post.
