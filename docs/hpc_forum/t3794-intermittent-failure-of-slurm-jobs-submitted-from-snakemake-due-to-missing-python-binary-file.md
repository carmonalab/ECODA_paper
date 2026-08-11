# Intermittent failure of slurm jobs submitted from snakemake due to "missing" python binary file

- Source: https://hpc-community.unige.ch/t/3794

- Created: 2025-01-22T15:41:26.984Z

- Tags: baobab

- Posts: 7

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Stephen.Mulligan (2025-01-22T15:41:27.038Z)

If you are asking for help, try to provide information that can help us solve your issue, such as :
I have a snakemake workflow which in the last ~month has started to result in failed jobs due to the following error
```
/var/spool/slurmd/job14377756/slurm_script: line 4: 
/home/users/m/mulligas/venvs/snakemake/bin/python: No such file or directory
```
# what did you try:
verified that the file does in fact exist. Re-ran the workflow, the problem often “fixes itself” and jobs stop failing due to this
# what didn’t work:
seemingly random slurm job failures as part of a larger workflow
# what was the expected result:
successful submission of slurm jobs and code execution
# what was the error message:
```
/var/spool/slurmd/job14377756/slurm_script: line 4: 
/home/users/m/mulligas/venvs/snakemake/bin/python: No such file or directory
```
# path to the relevant files (logs, sbatch script, etc):
```
/home/users/m/mulligas/calo_working/energy_study/modules/.snakemake/slurm_logs/
rule_gather_results/14373124.log
rule_gather_results/14373882.log
rule_gather_results/14374047.log
rule_gather_results/14374159.log
rule_gather_results/14374570.log
rule_gather_results/14377754.log
rule_gather_results/14377755.log
rule_gather_results/14377756.log
rule_gather_results/14377757.log
rule_gather_results/14377758.log
rule_gather_results/14377869.log
rule_gather_results/14377870.log
rule_gather_results/14377871.log
rule_gather_results/14377872.log
rule_gather_results/14377873.log
```


## Post 2 by @Adrien.Albert (2025-01-22T16:49:52.525Z)

Hi @Stephen.Mulligan[@Stephen.Mulligan](https://hpc-community.unige.ch/u/stephen.mulligan)
For the future request please post on HPC issues[HPC issues](https://hpc-community.unige.ch/c/hpc-support/hpc-issues/14),use the provided template and format  your text to make it readable by everyone. This help us a lot :pray:
Could you please provide your sbatch ?
Are you sure `/home/users/m/mulligas/venvs/snakemake/bin/python` is not deleted/edited/moved during your execution ?


## Post 3 by @Stephen.Mulligan (2025-01-23T09:41:54.451Z)

Hi @Adrien.Albert[@Adrien.Albert](https://hpc-community.unige.ch/u/Adrien.Albert) , thanks very much for the response. I’ll be sure to make future posts on Hpc issuse, sorry for that.
As I’m submitting jobs with snakemake the sbatch is not readily available to me. I looked very briefly as the snakemake documentation but couldn’t find a discussion of how this is handled. I assume snakemake generates some sbatch at some point to handle job submissions, maybe you are more familiar with this?
As for the deletion/editing of `bin/python`, I’m not sure that this is not happening and so it’s definitely possible, but there’s no real pattern I can see as to why this happens for some jobs and not others.


## Post 4 by @Adrien.Albert (2025-01-24T09:27:12.525Z)

I am asking you because you are using your own snakemake:
```
/home/users/m/mulligas/venvs/snakemake/bin/python
```
But it is available via module:
```
(baobab)-[alberta@login1 ~]$ ml spider snakemake

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  snakemake: snakemake/6.6.1
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Description:
      The Snakemake workflow management system is a tool to create reproducible and scalable data analyses.

    You will need to load all module(s) on any one of the lines below before the "snakemake/6.6.1" module is available to load.

      GCC/10.3.0  OpenMPI/4.1.1
 
    Help:
      Description
      ===========
      The Snakemake workflow management system is a tool to create reproducible and scalable data analyses.
      
      
      More information
      ================
       - Homepage: https://snakemake.readthedocs.io
      
      
      Included extensions
      ===================
      amply-0.1.4, ConfigArgParse-1.5.1, connection_pool-0.0.3, datrie-0.8.2,
      PuLP-2.4, ratelimiter-1.2.0.post0, smart_open-4.2.0, snakemake-6.6.1,
      stopit-1.1.2, toposort-1.6, wrapt-1.12.1
```
Is there a reason to use your own ?


## Post 5 by @Stephen.Mulligan (2025-01-24T11:10:36.999Z)

Sorry wasn’t aware this was available. I’ve just tested this out using this instead of my own venv as I was encountering the same issue again and I’m getting the following error:
```
snakemake: error: unrecognized arguments: --executor=slurm
```
which suggests to me that this is missing the snakemake-executor-plugin-slurm needed for use with slurm, and when i search for this with module spider I can’t find anything


## Post 7 by @Yann.Sagon (2025-01-28T11:45:11.478Z)

Dear @Stephen.Mulligan[@Stephen.Mulligan](https://hpc-community.unige.ch/u/stephen.mulligan)
We’ve installed a more recent snakemake version[snakemake version](https://hpc-community.unige.ch/t/new-software-installed-snakemake-version-8-4-2/3801).
This version has `snakemake-executor-plugin-slurm-0.2.1` included.
Please give a try.
Best
Yann


## Post 8 by @Stephen.Mulligan (2025-01-30T15:04:38.106Z)

Hi Yann, thanks for the response.
Using now and seems to be working well. Have not encountered the original issue so far, but this had also stopped with my own venv so remains to be seen if it crops up again.
Regardless, thanks for the help
