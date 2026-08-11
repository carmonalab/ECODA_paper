# How to transfer files from the Baobab cluster to local storage?

- Source: https://hpc-community.unige.ch/t/3742

- Created: 2024-11-25T09:07:32.009Z

- Tags: baobab

- Posts: 2

- Category: 8

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Michael.Silvestre (2024-11-25T09:07:32.052Z)

I’ve been using OpenOnDemand to run code with a graphical interface but recently have been encountering a lot of random crashes that make working on it impossible. I’ve used Baobab to deal with a very large dataset that I’ve since trimmed down by a lot.
How can I take my data from the cluster to my local storage ?
I tried downloading it from FileZilla but encountered an error in the process.
"
Erreur :
Could not open ‘C:\Users\MasterMergeRDD.dta’ for writing: Permission denied.
Erreur :
local: unable to open C:\Users\MasterMergeRDD.dta
Erreur :
Erreur critique lors du transfert du fichier
"
Thanks for the help,
Michaël


## Post 2 by @Adrien.Albert (2024-11-26T13:30:44.815Z)

Hi @Michael.Silvestre[@Michael.Silvestre](https://hpc-community.unige.ch/u/michael.silvestre)
You get a permission denied error, which means that you are not authorized to write to the destination.
PS: we resolve it during our zoom meeting.
