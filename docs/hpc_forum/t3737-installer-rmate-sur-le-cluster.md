# Installer rmate sur le cluster

- Source: https://hpc-community.unige.ch/t/3737

- Created: 2024-11-21T09:16:51.665Z

- Posts: 3

- Category: 10

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Theo.Desbordes (2024-11-21T09:16:51.711Z)

Bonjour,
Je souhaiterai éditer mes scripts sur le cluster depuis mon laptop avec sublime text, comme il est décri ici : Editing files remotely via SSH on SublimeText 3[Editing files remotely via SSH on SublimeText 3](https://wrgms.com/editing-files-remotely-via-ssh-on-sublimetext-3/)
Comme indiqué sur le tuto, il faut installer un script dans`/usr/local/bin/`
Cela serait-il posible ?
Merci d’avance !
Théo


## Post 2 by @Yann.Sagon (2024-11-21T09:23:46.918Z)

Dear @Theo.Desbordes[@Theo.Desbordes](https://hpc-community.unige.ch/u/theo.desbordes)
What is important is that the script you are talking about is present in your `$PATH`.
On our clusters, you can create a directory named `bin` in your home directory, it will be part of your `$PATH`
Here is an example to test:
Create the `bin` directory
```
(baobab)-[desborde@login1 ~]$ mkdir bin
```
Creat a dummy script in the `bin` directory
```
(baobab)-[desborde@login1 ~]$ cat  << EOF > bin/toto.sh
> #!/bin/sh
> echo toto
> EOF
```
Do not forget to make the script executable
```
(baobab)-[desborde@login1 ~]$ chmod a+x bin/toto.sh
```
And test it
```
(baobab)-[desborde@login1 ~]$ toto.sh
toto
```


## Post 3 by @Theo.Desbordes (2024-11-21T15:47:26.901Z)

Thanks very much @Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon)!
