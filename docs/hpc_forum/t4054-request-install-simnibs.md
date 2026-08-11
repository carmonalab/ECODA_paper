# Request - install SimNIBS

- Source: https://hpc-community.unige.ch/t/4054

- Created: 2025-08-20T12:02:27.792Z

- Posts: 11

- Category: 10

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Myriam.Trigilia (2025-08-20T12:02:27.850Z)

## Primary informations
Username: trigilia
Cluster: Baobab
## Description
Can you please assist me in installing SimNIBS ( SimNIBS 4.5 — SimNIBS 4.5.0 documentation[SimNIBS 4.5 — SimNIBS 4.5.0 documentation](https://simnibs.github.io/simnibs/build/html/index.html) ) in my hpc account? Or can you please show me if there are any steps that I can follow to install it myself?
Thank you


## Post 2 by @Yann.Sagon (2025-08-21T15:02:14.052Z)

Dear @Myriam.Trigilia[@Myriam.Trigilia](https://hpc-community.unige.ch/u/myriam.trigilia) this is done but version 4.0.1 as newer version isn’t yet available in our software repository. Is this working for you?
New software installed: SimNIBS version 4.0.1[New software installed: SimNIBS version 4.0.1](https://hpc-community.unige.ch/t/new-software-installed-simnibs-version-4-0-1/4055) HPC ChangeLog[HPC ChangeLog](https://hpc-community.unige.ch/c/hpc-announce/hpc-changelog/9)
> Dear users, we have installed a new software: SimNIBS 4.0.1: 
flatpak: symbol lookup error: /lib64/libldap.so.2: undefined symbol: EVP_md2, version OPENSSL_3.0.0

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  SimNIBS: SimNIBS/4.0.1
------------------------------------------------------------------------------…


## Post 3 by @Myriam.Trigilia (2025-08-22T13:21:38.413Z)

Dear @Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon),
Thank you for your help!
I tested SimNIBS/4.0.1 on my account. The main command (simnibs --version, simnibs --help) works and gmsh is available:
(baobab)-[trigilia@login1 ~]$ module load GCC/12.3.0
(baobab)-[trigilia@login1 ~]$ module load OpenMPI/4.1.5
(baobab)-[trigilia@login1 ~]$ module load SimNIBS/4.0.1
(baobab)-[trigilia@login1 ~]$ simnibs --version
4.0.1
(baobab)-[trigilia@login1 ~]$ gmsh --version
4.9.5
Is there something else that I need to check or is this fine?
Best regards,
Myriam


## Post 4 by @Yann.Sagon (2025-08-25T07:25:17.377Z)

Hello,
the next step is that to try to use it I guess:)
Best regards
Yann


## Post 5 by @Myriam.Trigilia (2025-08-26T10:29:24.999Z)

Hi Yann,
When testing the charm command in SimNIBS 4.0.1, the segmentation pipeline crashes at the post-processing step due to SimNIBS still using np.bool, which was removed in NumPy 1.24. Here’s the error:
Capture d’écran 2025-08-26 à 12.26.16
Capture d’écran 2025-08-26 à 12.26.161824×348 122 KB
[Capture d’écran 2025-08-26 à 12.26.161824×348 122 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/c18ffdcc63c8203770a95e4a7219074202a9a28d.png)
Would it be possible for you to downgrade NumPy to ≤1.23 (e.g. 1.22.0) so that SimNIBS runs without crashing?
Thank you so much!
Myriam


## Post 6 by @Yann.Sagon (2025-08-27T19:30:07.943Z)

Hi, thanks to Pavel’s commit[Pavel’s commit](https://github.com/easybuilders/easybuild-easyconfigs/pull/23748/commits/1bd5392972cdbb5fa8950909f8401405f9f1bf87), this should be fixed. Please give a try


## Post 7 by @Myriam.Trigilia (2025-08-28T11:17:22.743Z)

Dear Yann,
Thank you so much. It is working successfully, except that in the end Simnibs tried to save a report in the WEBP format:
[ simnibs ]INFO: Creating report
[ simnibs ]CRITICAL: Uncaught exception
KeyError: ‘WEBP’
(baobab)-[trigilia@login1 jobs]$
For completeness, can you please let me know if you are able to also add the WEBP plugin?
Thank you very very much!
Myriam


## Post 8 by @Yann.Sagon (2025-08-28T11:29:57.865Z)

Hi, can you try to simply load libwebp before running SimNIBS?
```
module load GCC/12.3.0 OpenMPI/4.1.5 SimNIBS/4.0.1 libwebp/1.3.1
```
Best


## Post 10 by @Myriam.Trigilia (2025-08-28T20:09:04.166Z)

Hi Yann,
I loaded the module as request:
Capture d’écran 2025-08-28 à 22.07.37
Capture d’écran 2025-08-28 à 22.07.37668×486 98.1 KB
[Capture d’écran 2025-08-28 à 22.07.37668×486 98.1 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/45510bbb2080cb6f100d9e16154341f35e710e6e.png)
But, I still have the same error:
[ simnibs ]CRITICAL: Uncaught exception
KeyError: ‘WEBP’
(baobab)-[trigilia@login1 jobs]$
Can you please recommend me another option?
Thank you a lot!
Myriam


## Post 11 by @Yann.Sagon (2025-08-29T06:38:30.548Z)

Hi, this should be ok now. The Pillow library was missing support for webp format. You can remove `libwebp/1.3.1` from the module load line and try again.


## Post 12 by @Myriam.Trigilia (2025-08-29T10:44:02.048Z)

Dear Yann,
Thank you so much! It worked!
Have a nice day!
Myriam
