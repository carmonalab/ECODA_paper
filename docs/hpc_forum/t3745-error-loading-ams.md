# Error loading AMS

- Source: https://hpc-community.unige.ch/t/3745

- Created: 2024-11-27T09:34:10.815Z

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Anais.Herren (2024-11-27T09:34:10.849Z)

Dear HPC team,
Username: herrena
Cluster: Bamboo
I am not able to load AMS when I type:
```
module load AMS/2023.105-intelmpi
```
I get this error:
```
These module(s) or extension(s) exist but cannot be loaded as requested: "AMS/2023.105-intelmpi"
```
Thank you in advance for your help,
Regards,
Anais


## Post 2 by @Adrien.Albert (2024-11-27T12:07:28.395Z)

Hi Anais,
Anais.Herren:
> I get this error:
```
These module(s) or extension(s) exist but cannot be loaded as requested: "AMS/20
```
If it’s not quite right, you get this error:
```
(bamboo)-[alberta@login1 ~]$ module load AMS/2023.105-intelmpi
Lmod has detected the following error:  These module(s) or extension(s) exist but cannot be loaded as requested: "AMS/2023.105-intelmpi"
   Try: "module spider AMS/2023.105-intelmpi" to see how to load the module(s).
```
The most important clue is this sentence:
```
Try: "module spider AMS/2023.105-intelmpi" to see how to load the module(s).
```
So let cluedo it:
```
(bamboo)-[alberta@login1 ~]$ module spider AMS/2023.105-intelmpi

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  AMS: AMS/2023.105-intelmpi
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Description:
      The Amsterdam Modeling Suite (AMS) provides a comprehensive set of modules for computational chemistry and materials science, from quantum mechanics to fluid thermodynamics. 

    You will need to load all module(s) on any one of the lines below before the "AMS/2023.105-intelmpi" module is available to load.

      intel-compilers/2022.2.1  impi/2021.7.1
 
    Help:
      Description
      ===========
      The Amsterdam Modeling Suite (AMS) provides a comprehensive set of modules for
      computational chemistry and materials science, from quantum mechanics to fluid
      thermodynamics.
      
      
      More information
      ================
       - Homepage: https://www.scm.com/amsterdam-modeling-suite/
```
here the next clue:
```
    You will need to load all module(s) on any one of the lines below before the "AMS/2023.105-intelmpi" module is available to load.

      intel-compilers/2022.2.1  impi/2021.7.1
```
So let cluedot it:
```
(bamboo)-[alberta@login1 ~]$ ml intel-compilers/2022.2.1  impi/2021.7.1 !$
ml intel-compilers/2022.2.1  impi/2021.7.1 AMS/2023.105-intelmpi
These environment variables need to be defined before using AMS:
  * $SCMLICENSE: path to AMS license file
  * $SCM_TMPDIR: path to user scratch directory
```
I recommend you read the entire output, it sometimes helps :wink:
Sorry, I couldn’t help myself. Oups I did it again… :notes:
