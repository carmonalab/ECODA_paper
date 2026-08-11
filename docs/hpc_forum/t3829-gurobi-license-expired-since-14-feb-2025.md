# Gurobi license expired since 14 Feb 2025

- Source: https://hpc-community.unige.ch/t/3829

- Created: 2025-02-17T09:09:52.498Z

- Tags: baobab

- Posts: 3

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Xin.Wen (2025-02-17T09:09:52.553Z)

Dear HPC team,
Hope you are doing well!
Please find below an issue regarding the Gurobi license.
Thank you very much in advance for your help!
Best,
Xin
## Primary informations
Username: $wenx
Cluster: $Baobab
## Description
It seems like Gurobi license has expired on Baobab since 14 Feb 2025.
e.g. job submitted: 15144577, 15144581
## Actual Result
Gurobi shell based on Python 3.9.6 can be launched with command `gurobi.sh`
Gurobi Python Interface can be loaded in Python 3.9.6 with ‘import gurobipy’
Running DEU on node cpu001.baobab
…
File “/home/users/w/wenx/baobab_python_env/lib/python3.9/site-packages/pyomo/solvers/plugins/solvers/gurobi_direct.py”, line 232, in available
raise ApplicationError(
pyomo.common.errors.ApplicationError: Could not create Model for <class ‘pyomo.solvers.plugins.solvers.gurobi_direct.GurobiDirect’> solver plugin - gurobi message=Could not create Model - gurobi message=Request denied: license expired
Set parameter GURO_PAR_SPECIAL
Set parameter TokenServer to value “admin1.baobab.hpc.unige.ch”
…


## Post 2 by @Adrien.Albert (2025-02-24T16:24:45.817Z)

Hello,
Gurobi is now up and running again. We apologize for any inconvenience caused during the disruption.
Please let us know if everything is working as expected.


## Post 3 by @Xin.Wen (2025-02-27T18:48:07.005Z)

Perfect, I tried on Bamboo and it works! Thank you very much :slight_smile:
