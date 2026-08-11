# Gurobi license expired since 25 Feb 2026

- Source: https://hpc-community.unige.ch/t/4232

- Created: 2026-02-25T06:39:49.548Z

- Tags: yggdrasil

- Posts: 4

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Xin.Wen (2026-02-25T06:39:49.610Z)

Dear HPC team,
Hope you are doing well!
Please find below an issue regarding the Gurobi license.
## Primary informations
Username: $wenx
Cluster: $Yggdrasil
## Description
It seems like Gurobi license has expired on Yggdrasil since 25 Feb 2026.
e.g. job submitted: 44003044, 44003045
You could also find below the information regarding the missing license by a quick check using gurobi.sh.
## Steps to Reproduce
$ ml GCCcore/12.3.0
$ ml Gurobi/11.0.2
$ gurobi.sh
## Expected Result
license found
## Actual Result
Python 3.11.3 (main, Nov  4 2024, 14:00:24) [GCC 12.3.0] on linux
Type “help”, “copyright”, “credits” or “license” for more information.
Set parameter LicenseID to value 218116
Set parameter GURO_PAR_SPECIAL
Set parameter TokenServer to value “admin1.baobab”
Set parameter LogFile to value “gurobi.log”
Request denied: license expired
Thank you very much for your help in advance !
Best regards,
Xin


## Post 2 by @Adrien.Albert (2026-02-25T10:00:03.025Z)

Hello;
Sorry for the inconvenience, the licence extension has been requested and should be updated soon.


## Post 3 by @Adrien.Albert (2026-02-25T11:06:48.300Z)

The Gurobi license has been updated. I have also added a monitoring check to notify me when the license is about to expire, in order to avoid any interruption.
I apologize once again for the inconvenience caused.


## Post 4 by @Xin.Wen (2026-02-25T14:26:02.532Z)

Perfect, thank you so much for your prompt help :slight_smile: !
