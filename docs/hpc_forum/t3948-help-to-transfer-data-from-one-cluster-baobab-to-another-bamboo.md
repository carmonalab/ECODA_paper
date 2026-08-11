# Help to transfer data from one cluster (Baobab) to another (Bamboo)?

- Source: https://hpc-community.unige.ch/t/3948

- Created: 2025-05-07T15:38:55.342Z

- Tags: baobab, bamboo

- Posts: 3

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Nicolas.Clairis1 (2025-05-07T15:38:55.390Z)

Due to insufficient RAM to run my jobs, I was recommended to switch from Baobab to Bamboo (conversation here[conversation here](https://hpc-community.unige.ch/)) which has more RAM available, but I don’t know how to do such a transfer and would like to avoid to retransfer the data locally in nasac as it would take very long to do so given the size of the data + I may not have enough space on nasac to store all the temporary files that are necessary to run my analysis.
Can you guide me through the best way to do such a transfer?
I just tried to use rsync following hpc:best_practices [eResearch Doc][hpc:best_practices [eResearch Doc]](https://doc.eresearch.unige.ch/hpc/best_practices?s%5B%5D=rsync#rsync) but it doesn’t seem to work (maybe because I try to transfer a folder from scratch?). Here is the command that I launched:
```
(baobab)-[clairis@login1 scratch]$ DST=$/home/users/c/clairis/scratch
(baobab)-[clairis@login1 scratch]$ DIR=fMRI_analysis
(baobab)-[clairis@login1 scratch]$ BAMBOO=login1.bamboo
(baobab)-[clairis@login1 scratch]$ rsync -aviuzPrg ${DIR} ${BAMBOO}:${DST}
(clairis@login1.bamboo) Password:
sending incremental file list
rsync: [Receiver] mkdir "/home/users/c/clairis/$/home/users/c/clairis/scratch" failed: No such file or directory (2)
rsync error: error in file IO (code 11) at main.c(783) [Receiver=3.2.3]
```
Thanks in advance for any help!


## Post 3 by @Adrien.Albert (2025-05-09T14:00:47.890Z)

There is a typo in DST definition:
```
rsync: [Receiver] mkdir "/home/users/c/clairis/$/home/users/c/clairis/scratch" failed: No such file or directory (2)
# ---------------------------------------------^ Here
```
```
(baobab)-[clairis@login1 scratch]$ DST=$/home/users/c/clairis/scratch
# -------------------------------------^ Remove the '$'
```
let me know it’s working :slight_smile:


## Post 4 by @nicolas.clairis (2025-05-09T17:07:17.024Z)

Yes it seems to work now!!! Thanks a lot :folded_hands:
