# SSH Connection to Baobab

- Source: https://hpc-community.unige.ch/t/4083

- Created: 2025-09-11T16:36:56.594Z

- Tags: baobab

- Posts: 3

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Marco.Froelich (2025-09-11T16:36:56.685Z)

Hi,
I am having issues connecting via SSH from VSCode to Boabab login node since the maintenance (last successful login was the morning of the maintenance start). What is attempted:
SSH to: login1.baobab.hpc.unige.ch, as  at Port 22 (tried 80 also, to align with potential maintenance update).
There is long attempt at connection with no error appearing. After 5min or so:
```
Error: CodeError(ServerOriginTimeout)
```
Connection through FileZilla works fine. Any ideas what could be wrong?
Thanks for your work,
Marco
Update: I was able to retrieve another error while trying to connect:
```
Failed to connect to the remote extension host server (Error: WrappedError(WrappedError { message: “error creating temp download dir”, original: “Permission denied (os error 13) at path "/tmp/.tmpQkZ2gS"” }))
```


## Post 2 by @Marco.Froelich (2025-09-15T07:10:35.354Z)

Update: My terminal is showing errors more precisely now. It stems from Downloading the VSCode Server:
Error installing server: error creating temp download dir: Permission denied (os error 13) at path "/tmp/.tmpSvLYFX”
Any help is appreciated :slight_smile:


## Post 3 by @Yann.Sagon (2025-09-15T07:18:31.719Z)

Hi @Marco.Froelich[@Marco.Froelich](https://hpc-community.unige.ch/u/marco.froelich) this is now fixed, this was the same issue as Cannot connect to baobab via X2Go[Cannot connect to baobab via X2Go](https://hpc-community.unige.ch/t/cannot-connect-to-baobab-via-x2go/4084)
Best
Yann
