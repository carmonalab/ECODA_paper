# Cannot connect to baobab via X2Go

- Source: https://hpc-community.unige.ch/t/4084

- Created: 2025-09-12T08:47:29.114Z

- Tags: baobab

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @William.Ceva (2025-09-12T08:47:29.179Z)

Hello,
Following the changes + reboot on baobab (i.e. as of yesterday, when baobab became available again), I am unable to connect to baobab via X2Go.  I can connect via normal ssh, but not with X2Go.
I confirm that the X2Go connection settings I am using to attempt to connect to baobab are identical to those I am using to connect to yggdrasil and bamboo (to both of which I can connect to using X2Go without issue).
I also confirm that I have already used an ssh connection to move my .bashrc file into my ~/tmp/ directory, then try to connect to baobab via X2Go, and yet this still does not solve the issue.  Nor does moving my .bashrc back into my home directory and then trying to connect with X2Go again; I confirm that my .bashrc file is not the issue either.
Below is a screenshot with the error message I receive when I attempt to connect to baobab via X2Go.
Any help solving this issue is appreciated,
Will Ceva
Screenshot 2025-09-12 104507
Screenshot 2025-09-12 104507565×668 47.5 KB
[Screenshot 2025-09-12 104507565×668 47.5 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/8da4c6b7b3309395c7de5a62ab0a8895b61795b9.png)


## Post 2 by @Yann.Sagon (2025-09-15T07:16:08.612Z)

Dear @William.Ceva[@William.Ceva](https://hpc-community.unige.ch/u/william.ceva) this is fixed, thanks for the feedback.
By the way, we’re not sure whether to continue supporting x2go, given that it’s now possible to open a remote graphical session through OpenOnDemand. Would using OOD be a solution for you?
Best regards
Yann
