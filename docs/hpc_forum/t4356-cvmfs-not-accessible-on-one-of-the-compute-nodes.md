# /cvmfs not accessible on one of the compute nodes

- Source: https://hpc-community.unige.ch/t/4356

- Created: 2026-07-21T06:57:27.980Z

- Tags: baobab

- Posts: 1

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Paul.Coppin (2026-07-21T06:57:28.049Z)

## Primary informations
Username: coppinp
Cluster: Baobab
## Description
Jobs fail on cpu350 with the message:
`cd: /cvmfs/dampe.cern.ch/centos7/opt//releases/DmpSoftware-6-0-18_geant4_10_06p03/: No such file or directory`
## Steps to Reproduce
```
#!/bin/bash
REPO=/cvmfs/dampe.cern.ch
NODE=$(hostname -s)
if [ ! -d "$REPO" ]; then
    echo "$NODE NO_DIR"
elif timeout 20 ls "$REPO" >/dev/null 2>&1; then
    echo "$NODE OK"
else
    echo "$NODE FAIL_OR_TIMEOUT"
fi
```
## Expected Result
I checked a few other nodes in the private-dpnc-cpu queue, which all seem fine. Just cpu350 return the output: cpu350 NO_DIR
