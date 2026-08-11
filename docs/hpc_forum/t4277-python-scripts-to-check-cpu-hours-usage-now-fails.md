# Python scripts to check CPU hours usage now fails

- Source: https://hpc-community.unige.ch/t/4277

- Created: 2026-04-16T08:48:26.243Z

- Tags: slurm

- Posts: 6

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @David.Vicente (2026-04-16T08:48:26.325Z)

I was tasked to check the credit usage of the HPC clusters for my group. For this I was using a script provided to me (i assume it originated from the HPC team) to compute used CPU hours and how many credits we have remaining:
```
#!/bin/bash

module load GCCcore/13.3.0 ClusterShell/1.9.3 tabulate2/1.10.0

# total credits

echo 'Total credits'

echo "==========================================================="

ug_getNodeCharacteristicsSummary.py --partitions private-teodoro-gpu private-teodoro-cpu --cluster baobab --summary

echo "-----------------------------------------------------------"

# used credits

echo 'Used credits'

echo "==========================================================="

ug_slurm_usage_per_user.py --group private_teodoro --start 2026-01-01 --report-type account

echo "-----------------------------------------------------------"

echo 'Total DRIM credits'

# total DRIM credits

echo "==========================================================="

ug_getNodeCharacteristicsSummary.py --partitions private-drim-gpu --cluster bamboo --summary

echo "-----------------------------------------------------------"

# used credits

echo "DRIM used credits"

echo "==========================================================="

ug_slurm_usage_per_user.py --group private_drim --start 2026-01-01 --report-type account

echo "-----------------------------------------------------------
```
It used to work fine, but I tried to run it this week to report it, but now it fails. Here is the output:
```
Total credits
===========================================================
/bin/env: ‘$UG_PYTHON’: No such file or directory
-----------------------------------------------------------
Used credits
===========================================================
/usr/bin/env: ‘$UG_PYTHON’: No such file or directory
-----------------------------------------------------------
Total DRIM credits
===========================================================
/bin/env: ‘$UG_PYTHON’: No such file or directory
-----------------------------------------------------------
DRIM used credits
===========================================================
/usr/bin/env: ‘$UG_PYTHON’: No such file or directory
-----------------------------------------------------------
```
I guess the python path is not in my ENV anymore, I tried to load other modules, but it was not successful. Did something change in those scrips ?
have a nice day !
David Vicente


## Post 2 by @Gael.Rossignol (2026-04-16T14:42:31.571Z)

Hello,
The scripts have been improved, and using `module load` is no longer required.
Could you please try removing that line? It should work without it.
Best regards,


## Post 3 by @David.Vicente (2026-04-16T15:09:45.988Z)

Hi,
Thank for your reply. I tried to comment out the load module line or running the script directly, but it still fails. here is the command history:
```
(baobab)-[vincente@login1 ~]$ cat credits_check.sh 
#!/bin/bash

# module load GCCcore/13.3.0 ClusterShell/1.9.3 tabulate2/1.10.0

# total credits
echo 'Total credits'
echo "==========================================================="
ug_getNodeCharacteristicsSummary.py --partitions private-teodoro-gpu private-teodoro-cpu --cluster baobab --summary
echo "-----------------------------------------------------------"

# used credits
echo 'Used credits'
echo "==========================================================="
ug_slurm_usage_per_user.py --group private_teodoro --start 2026-01-01 --report-type account
echo "-----------------------------------------------------------"

echo 'Total DRIM credits'
# total DRIM credits
echo "==========================================================="
ug_getNodeCharacteristicsSummary.py --partitions private-drim-gpu --cluster bamboo --summary
echo "-----------------------------------------------------------"

# used credits 
echo "DRIM used credits"
echo "==========================================================="
ug_slurm_usage_per_user.py --group private_drim --start 2026-01-01 --report-type account
echo "-----------------------------------------------------------"
(baobab)-[vincente@login1 ~]$ ./credits_check.sh 
Total credits
===========================================================
/bin/env: ‘$UG_PYTHON’: No such file or directory
-----------------------------------------------------------
Used credits
===========================================================
/usr/bin/env: ‘$UG_PYTHON’: No such file or directory
-----------------------------------------------------------
Total DRIM credits
===========================================================
/bin/env: ‘$UG_PYTHON’: No such file or directory
-----------------------------------------------------------
DRIM used credits
===========================================================
/usr/bin/env: ‘$UG_PYTHON’: No such file or directory
-----------------------------------------------------------
(baobab)-[vincente@login1 ~]$ ug_getNodeCharacteristicsSummary.py --partitions private-teodoro-gpu private-teodoro-cpu --cluster baobab --summary
/bin/env: ‘$UG_PYTHON’: No such file or directory
(baobab)-[vincente@login1 ~]$ 
```
Have a nice afternoon !


## Post 4 by @Gael.Rossignol (2026-04-17T07:47:10.358Z)

Dear David,
You are right — you now need to use the wrapper without the `.py` extension.
Here is a working example:
```
(baobab)-[rossigng@login1 ~]$ ug_slurm_usage_per_user
--------------------------------------------------------------------------------

Cluster/User/Account Utilization 2026-04-01T00:00:00 - 2026-04-17T09:59:59 (1418400 secs)

Usage reported in TRES Hours

--------------------------------------------------------------------------------

Cluster    Login     Proper Name     Account    TRES Name      Used
---------  --------  --------------  ---------  -----------  ------
bamboo     rossigng  Rossignol Gaël  burgi      billing           3
yggdrasil  rossigng  Rossignol Gaël  burgi      billing           0
Total usage: 0.00M
```
Best regards,


## Post 5 by @David.Vicente (2026-04-17T08:03:39.680Z)

Hi,
The script now works ! Thank you for your time and support.
best


## Post 6 by @Marco.Froelich (2026-06-02T10:58:30.778Z)

Gael.Rossignol:
> `baobab)-[rossigng@login1 ~]$ ug_slurm_usage_per_user`
Note that the command needs to be run on login nodes.
