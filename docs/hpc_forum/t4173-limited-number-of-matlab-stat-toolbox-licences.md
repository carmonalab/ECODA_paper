# Limited number of matlab stat toolbox licences?

- Source: https://hpc-community.unige.ch/t/4173

- Created: 2025-12-17T10:40:59.491Z

- Tags: bamboo

- Posts: 2

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @nicolas.clairis (2025-12-17T10:40:59.535Z)

Many of my batches crashed yesterday with the following message:
```
License checkout failed.
License Manager Error -4
Maximum number of users for Statistics_Toolbox reached.
Try again later.

Troubleshoot this issue by visiting:
https://www.mathworks.com/support/lme/R2022a/4

Diagnostic Information:
Feature: Statistics_Toolbox
License path:

/home/users/c/clairis/.matlab/R2022a_licenses:/opt/ebsofts/MATLAB/2022a/licenses/license.dat:/opt/ebsofts/MATLAB/2022a/licenses/network.lic
Licensing error: -4,132."
```
I guess this comes from the fact that matlab/unige only provides with a limited number of pcs which an use the toolbox in parallel, but I was wondering if there isn’t a way to bypass that issue as this clearly limits the interest of using the cluster?


## Post 2 by @Gael.Rossignol (2025-12-18T13:34:07.065Z)

Dear Nicolas,
Regarding Matlab this is needed to reserve licence in the sbatch command. All is documented here :
doc.eresearch.unige.ch[doc.eresearch.unige.ch](https://doc.eresearch.unige.ch/hpc/applications_and_libraries?s%5B%5D=matlab#matlab)
### hpc:applications_and_libraries [eResearch Doc][hpc:applications_and_libraries [eResearch Doc]](https://doc.eresearch.unige.ch/hpc/applications_and_libraries?s%5B%5D=matlab#matlab)
But if you need to launch multiple instance of you script you mush compile matlab
doc.eresearch.unige.ch[doc.eresearch.unige.ch](https://doc.eresearch.unige.ch/hpc/applications_and_libraries?s%5B%5D=matlab#compile_your_matlab_code)
### hpc:applications_and_libraries [eResearch Doc][hpc:applications_and_libraries [eResearch Doc]](https://doc.eresearch.unige.ch/hpc/applications_and_libraries?s%5B%5D=matlab#compile_your_matlab_code)
After that you are not limited with the number of licence.
Best regards,
