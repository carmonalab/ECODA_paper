# Can't connect on bamboo with ssh key nor password

- Source: https://hpc-community.unige.ch/t/4317

- Created: 2026-06-18T11:21:18.687Z

- Tags: bamboo

- Posts: 5

- Category: 1

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Marc.Gay-Balmaz (2026-06-18T11:21:18.748Z)

I didn’t connect for some time (maybe a month), and when I tried back today, it says me that there are too many authentication failures. I didn’t change my ssh key, but probably my password. I also tried to connect forcing a password use, but the new password doesn’t work either. Do I ask the question at the right place ?


## Post 2 by @Yann.Sagon (2026-06-26T08:00:16.771Z)

Dear Marc, sorry we didn’t answered before.
In the meantime, it is now mandatory to connect using an ssh-key:)
I’ve checked: your account is still valid (at UNIGE and HPC).
Your keys are:
- ssh-ed25519 AAAAC3NzaC1lZD…
- ssh-rsa AAAAB3NzaC1yc2EAAAADAQABAAABgQ…
Please try again and let us know if not working.
- I tried on cluster : `<cluster>`
- I’m using ssh client: `<name> <version>`
- My OS is: `<name>` `<version>`
- The command line I’m using to connect (if any) `<command line>`


## Post 3 by @Martin.Millon (2026-06-26T08:16:49.228Z)

I am having a similar issue. Is the login node still login1.baobab.hpc.unige.ch. I found also baobab2.hpc.unige.ch in the documentation. Which one is correct ?


## Post 4 by @Yann.Sagon (2026-06-26T08:25:53.298Z)

Martin.Millon:
> I found also baobab2.hpc.unige.ch in the documentation. Which one is correct ?
Only login1 is working. Where did you found this information in the documentation please so we can correct it.
If you need help please give us information:
- I tried on cluster : `<cluster>`
- I’m using ssh client: `<name>` `<version>`
- My OS is: `<name>` `<version>`
- The command line I’m using to connect (if any) `<command line>`


## Post 5 by @Martin.Millon (2026-06-26T08:33:38.173Z)

Martin.Millon:
> baobab2.hpc.unige.ch
Here : Le cluster Baobab - Scientific and Parallel Computing Group - UNIGE[Le cluster Baobab - Scientific and Parallel Computing Group - UNIGE](https://spc.unige.ch/en/teaching/courses/parallelisme/le-cluster-baobab/)
- I tried on cluster : `login1.baobab.hpc.unige.ch`
- I’m using ssh client: OpenSSH_10.2p1, LibreSSL 3.3.6
- My OS is: `MacOS Tahoe26.5`
- The command line I’m using to connect (if any) `:`ssh -i ~/.ssh/id_ed25519_baobab millon@login1.baobab.hpc.unige.ch
