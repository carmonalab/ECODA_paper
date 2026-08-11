# Ftp does not work on baobab

- Source: https://hpc-community.unige.ch/t/3668

- Created: 2024-10-02T12:43:10.275Z

- Posts: 6

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Denis.Mongin (2024-10-02T12:43:10.310Z)

## Primary informations
Username: mongin
Cluster: baobab
## Description
Dear baobab team
FTP file transfer does not work anymore for me.
I use winscp, which works usually fine. I gives me an error when trying to send a file to `/home/mongin/participant_number` (I was able to do 5 min ago):
```
Code d'erreur : 4
Message d'erreur du serveur  : Failure

Les raisons courante pour l'erreur code 4 sont:
- Changement de nom de fichier avec un déjà existant.
- Création d'un dossier déjà existant.
- Déplacement de fichier distant vers un disque fichiers système (disque dur) différent.
- Envoyer un fichier vers un disque fichier système (disque dur) plein.
- Dépassement du quota disque utilisateur.
```
I tried to remove some files to have less space used (removes at least 1 Go), but it still does not work, it seems that disk use is not the problem.
This happend with all the files.
I tried closing and reopening Winscp, did not work either. I updated winscp, did not change anything.
I tried to install filezilla, but I cannot on my HUG laptop (it is blocked by the HUG).
Please help, I am totally blocked


## Post 2 by @Denis.Mongin (2024-10-02T13:06:04.598Z)

Ok, it looks like it is back again, now it is working.
Not sure what caused this though.


## Post 3 by @Yann.Sagon (2024-10-04T13:07:20.868Z)

Dear @Denis.Mongin[@Denis.Mongin](https://hpc-community.unige.ch/u/denis.mongin)
first of all, a technical detail: we aren’t using ftp, but sftp, two very distinct protocols.
I checked your quota on Baobab: you are near the limit, maybe this was the reason?
Best
Yann


## Post 4 by @Denis.Mongin (2024-10-10T07:18:43.313Z)

Thanks
How can I check the quota detail, to know where I should remove things ?


## Post 5 by @Denis.Mongin (2024-10-10T07:25:43.794Z)

Hi @Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon)
It is happeniing again this morning, and my work is blocked.
I removed large folders using `rm`, but the sftp still does not work.
Could you help ?


## Post 6 by @Yann.Sagon (2024-10-10T10:02:14.981Z)

Denis.Mongin:
> How can I check the quota detail, to know where I should remove things ?
You can get your quota usage like this:
```
(baobab)-[mongin@login1 ~]$ beegfs-get-quota-home-scratch.sh
home dir: /home/mongin
scratch dir: /srv/beegfs/scratch/users/m/mongin

          user/group                 ||           size          ||    chunk files
  storage     |   name        |  id  ||    used    |    hard    ||  used   |  hard
  ----------------------------|------||------------|------------||---------|---------
home        |         mongin|284413||  557.20 GiB| 1024.00 GiB||   163705|unlimited
scratch     |         mongin|284413||      0 Byte|   unlimited||        0| 1000000
```
So you did a lot of cleanup and this isn’t anymore the issue. I’ll contact you to do a zoom to debug the issue.
