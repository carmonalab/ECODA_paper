# Gurobi 9.5.0 license not available

- Source: https://hpc-community.unige.ch/t/4217

- Created: 2026-02-09T08:48:57.052Z

- Posts: 3

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Cecilia.Mannik (2026-02-09T08:48:57.100Z)

## Primary informations
Username: mannikc
Cluster: Yggdrasil
## Description
Hello, I am trying to run an optimization model that makes use of Gurobi 9.5.0, however the associated license seems to not be available, including when launching Gurobi as standalone.
## Steps to Reproduce
Run following commands:
```
module load GCCcore/11.2.0

module load Python/3.9.6

module load Gurobi/9.5.0

gurobi.sh
```
## Expected Result
The license is loaded and a Gurobi shell is launched.
## Actual Result
When running `gurobi.sh:`
```
Python 3.9.6 (default, Dec 17 2021, 16:47:27) 
[GCC 11.2.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
Set parameter WLSAccessID
Set parameter WLSSecret
Set parameter LicenseID
Set parameter LogFile to value "gurobi.log"

Web license service only available for container environments
```
When running `gurobi_cl`:
```
Set parameter WLSAccessID
Set parameter WLSSecret
Set parameter LicenseID
Set parameter LogFile to value "gurobi.log"

Failed to set up a license

Error 10024: Web license service only available for container environments
```


## Post 2 by @Adrien.Albert (2026-02-10T09:57:04.151Z)

Hello @Cecilia.Mannik[@Cecilia.Mannik](https://hpc-community.unige.ch/u/cecilia.mannik)
Is there a specific reason for using an older version of Gurobi?
With the new gurobi :
```
(baobab)-[alberta@login1 ~]$ ml GCCcore/12.3.0 Gurobi
Gurobi shell based on Python 3.11.3 can be launched with command `gurobi.sh`
Gurobi Python Interface can be loaded in Python 3.11.3 with 'import gurobipy'

(baobab)-[alberta@login1 ~]$ gurobi.sh
Python 3.11.3 (main, Nov  4 2024, 14:00:24) [GCC 12.3.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
Set parameter LicenseID to value 218116
Set parameter GURO_PAR_SPECIAL
Set parameter TokenServer to value "admin1.baobab"
Set parameter LogFile to value "gurobi.log"

Gurobi Interactive Shell (linux64), Version 11.0.2
Copyright (c) 2024, Gurobi Optimization, LLC
Type "help()" for help

gurobi> 

(baobab)-[alberta@login1 ~]$ gurobi_cl --tokens

Checking status of Gurobi token server 'admin1.baobab'...

Token server functioning normally.
Maximum allowed uses: 4096, current: 0

Found 1 active servers
```
In the meantime it has been fixed but   we will review the relevance of keeping older versions of this software.


## Post 3 by @Cecilia.Mannik (2026-02-10T14:21:40.181Z)

Hello,
I can confirm it is working now, many thanks!
The reason for using an older version of Gurobi is that newer ones do not seem to be compatible with the Python setup that is known to work for running the model (in particular an older version of Pyomo is required).
