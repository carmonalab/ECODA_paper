# Problème d'accès et de copie de fichier sur mes partitions baobab

- Source: https://hpc-community.unige.ch/t/4330

- Created: 2026-06-29T11:23:10.875Z

- Tags: baobab

- Posts: 3

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Bari.Dina (2026-06-29T11:23:10.939Z)

Bonjour,
Je vous écrit concernant un soucis que j’ai rencontré récemment avec l’utilisation de baobab. Depuis le changement de méthode d’accès avec la clé SSH je n’arrive pas à copier des fichiers dans les partitions auxquelles j’ai accès. Comme vous le verrez ci-joint, j’arrive à me connecter au baobab, mais me voit refuser la permission de copier des fichier sur des dossiers auxquels j’ai accès (ce n’était pas un problème avant la connection avec la clé SSH)
Capture d’écran 2026-06-29 à 13.16.54
Capture d’écran 2026-06-29 à 13.16.542116×1488 409 KB
[Capture d’écran 2026-06-29 à 13.16.542116×1488 409 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/9db1484e4fa96a4729e068ea37b5b45301844f20.jpeg)


## Post 2 by @Gael.Rossignol (2026-06-29T12:57:15.287Z)

Bonjour,
D’après la copie d’écran vous essayez de copier des fichiers depuis votre poste de travail vers le cluster HPC.
De ce fait il vous fait faire la commande “scp” et mettre votre clef privée dans l’option “-i”.
Depuis votre poste :
```
scp -r -i /*** chemin vers la clefs privée ***/ /Users/baridina/PREM/Videos/7958C dina2@login1.baobab-hpc-unige.ch:/srv/beegfs/scratch/shares/schaerm/schaer_sensit/CIPA_VIDEOS_2/
```
Dans tous les cas les droits sur les répertoires ne sont pas inpactés par la nouvelle procédure de connection.
Cordialement,


## Post 3 by @Bari.Dina (2026-06-29T13:45:45.229Z)

Bonjour Gaël,
Merci beaucoup ! J’ai réessayé depuis mon terminal local en spécifiant le chemin de ma clé SSH et le transfert a parfaitement fonctionné.
Belle journée à vous,
Bari
