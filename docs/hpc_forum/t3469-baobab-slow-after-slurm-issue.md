# Baobab slow after Slurm issue

- Source: https://hpc-community.unige.ch/t/3469

- Created: 2024-06-03T14:19:41.544Z

- Posts: 6

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Malte.Algren (2024-06-03T14:19:41.579Z)

Hi, since the slurm crashed I have experienced a significant decrease in the speed of baobab as well as a very slow internet connection.
Is it possible parts of the cluster are not fully up and running yet?
Sorry for the vague description but haven’t changed anything in my scripts since Friday and the cluster has been unresponsive today.
Best,
Malte


## Post 2 by @Imahn.Shekhzadeh1 (2024-06-03T16:02:36.876Z)

I was also experiencing a very slow internet connection (at least 7x slower than usual).


## Post 3 by @Samuel.Klein (2024-06-04T06:15:59.098Z)

I’m having a similar issue, when I try to debug in VSCode it always times out.


## Post 4 by @Imahn.Shekhzadeh1 (2024-06-11T13:36:39.108Z)

Dear @support[@support](https://hpc-community.unige.ch/groups/support), any update on the slow internet speed on `baobab`?


## Post 5 by @Yann.Sagon (2024-06-11T15:03:22.929Z)

Hello,
Maybe for this reason? a lot of vscode running/sleeping, using all the ram.
If a vscode specialist can let us know what’s happening, this would be great!
image
image1890×1000 193 KB
[image1890×1000 193 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/8ca9c15ea89227cf95ec6c155a4fe923d58161d3.png)
Maybe worth to try to upgrade vscode ? Number of processes increases every time remote window is reloaded · Issue #211462 · microsoft/vscode · GitHub[Number of processes increases every time remote window is reloaded · Issue #211462 · microsoft/vscode · GitHub](https://github.com/microsoft/vscode/issues/211462)


## Post 6 by @Ludovic.Dumoulin (2024-06-11T15:26:49.681Z)

Hello,
Is it the login node ?
If it is, it might be possible that some people are using the loggin node as vs code server. But I don’t know much about it. When I was running a vscode server on my workstation I saw kind of similar “.vscode-server/cli/servers/…” lines. I don’t remember exactly.
> The Visual Studio Code Server is a service you can run on a remote development machine, like your desktop PC or a virtual machine (VM). It allows you to securely connect to that remote machine from anywhere through a local VS Code client, without the requirement of SSH.
