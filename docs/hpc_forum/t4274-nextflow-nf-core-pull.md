# Nextflow - nf-core - pull

- Source: https://hpc-community.unige.ch/t/4274

- Created: 2026-04-13T13:18:11.926Z

- Posts: 4

- Category: 13

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Camille.Christe (2026-04-13T13:18:11.984Z)

I try to use Nextflow program and I would like to pull a script from nf-core repository.
It worked, I can see it is present typing nextflow list, real path points to the working directory I would like the script to be in but I cannot see it.
> nextflow pull nf-core/demo
> Checking nf-core/demo …
> downloaded from GitHub - nf-core/demo: nf-core/demo is a simple nf-core style bioinformatics pipeline for workshops and demos. · GitHub[GitHub - nf-core/demo: nf-core/demo is a simple nf-core style bioinformatics pipeline for workshops and demos. · GitHub](https://github.com/nf-core/demo.git) - revision: 45904cb9d1 [master]
> realpath demo.git
> /home/christec/Nextflow/demo.git
How to access a script pulled from nf-core repository?
Best regards,
Camille


## Post 2 by @Camille.Christe (2026-04-13T15:28:49.532Z)

Related to that, the demo script for this pipeline is using containers to set up environment. There is a profile for Apptainer but I was not able to set it up to make it work.
Is there anything to know when using containers in program such as nextflow or snakemake?
Best,
Camille


## Post 4 by @Gael.Rossignol (2026-04-14T12:20:37.669Z)

Dear Camille,
After reading the documentation I follow the procedure and it seems to be working fine :
```
(base) (bamboo)-[rossigng@login1 ~]$ salloc --cpus-per-task=6
salloc: Pending job allocation 3850673
salloc: job 3850673 queued and waiting for resources
salloc: job 3850673 has been allocated resources
salloc: Granted job allocation 3850673
salloc: Nodes cpu001 are ready for job
(base) (bamboo)-[rossigng@cpu001 ~]$ nextflow run nf-core/demo -profile singularity --input samplesheet.csv --outdir .
bash: nextflow: command not found
(base) (bamboo)-[rossigng@cpu001 ~]$ ml Nextflow/25.10.4
(base) (bamboo)-[rossigng@cpu001 ~]$ nextflow run nf-core/demo -profile singularity --input samplesheet.csv --outdir .

 N E X T F L O W   ~  version 25.10.4

Launching `https://github.com/nf-core/demo` [spontaneous_curry] DSL2 - revision: 45904cb9d1 [master]

------------------------------------------------------
                                        ,--./,-.
        ___     __   __   __   ___     /,-._.--~'
  |\ | |__  __ /  ` /  \ |__) |__         }  {
  | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                        `._,._,'
  nf-core/demo 1.1.0
------------------------------------------------------

Input/output options
  input              : samplesheet.csv
  outdir             : .

Generic options
  trace_report_suffix: 2026-04-14_14-20-25

Core Nextflow options
  revision           : master
  runName            : spontaneous_curry
  containerEngine    : singularity
  launchDir          : /home/users/r/rossigng
  workDir            : /home/users/r/rossigng/work
  projectDir         : /home/users/r/rossigng/.nextflow/assets/nf-core/demo
  userName           : rossigng
  profile            : singularity
  configFiles        : /home/users/r/rossigng/.nextflow/assets/nf-core/demo/nextflow.config

!! Only displaying parameters that differ from the pipeline defaults !!
------------------------------------------------------

* The pipeline
    https://doi.org/10.5281/zenodo.12192442

* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/nf-core/demo/blob/master/CITATIONS.md

executor >  local (6)
[16/d32ac7] NFCORE_DEMO:DEMO:FASTQC (SAMPLE2_PE)     | 3 of 3 ✔
[95/141ce2] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE2_PE) | 3 of 3 ✔
[-        ] NFCORE_DEMO:DEMO:MULTIQC                 -
Pulling Singularity image https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/34/34e733a9aexecutor >  local (7)
[16/d32ac7] NFCORE_DEMO:DEMO:FASTQC (SAMPLE2_PE)     | 3 of 3 ✔
[95/141ce2] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE2_PE) | 3 of 3 ✔
[a0/19d0ad] NFCORE_DEMO:DEMO:MULTIQC                 | 0 of 1
Pulling Singularity image https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/34/34e733a9aexecutor >  local (7)
[16/d32ac7] NFCORE_DEMO:DEMO:FASTQC (SAMPLE2_PE)     | 3 of 3 ✔
[95/141ce2] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE2_PE) | 3 of 3 ✔
[a0/19d0ad] NFCORE_DEMO:DEMO:MULTIQC                 | 0 of 1
Pulling Singularity image https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/34/34e733a9aexecutor >  local (7)
[16/d32ac7] NFCORE_DEMO:DEMO:FASTQC (SAMPLE2_PE)     | 3 of 3 ✔
[95/141ce2] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE2_PE) | 3 of 3 ✔
[a0/19d0ad] NFCORE_DEMO:DEMO:MULTIQC                 | 1 of 1 ✔
Pulling Singularity image https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/34/34e733a9aexecutor >  local (7)
[16/d32ac7] NFCORE_DEMO:DEMO:FASTQC (SAMPLE2_PE)     | 3 of 3 ✔
[95/141ce2] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE2_PE) | 3 of 3 ✔
[a0/19d0ad] NFCORE_DEMO:DEMO:MULTIQC                 | 1 of 1 ✔
Pulling Singularity image https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/34/34e733a9ae16a27e80fe00f863ea1479c96416017f24a907996126283e7ecd4d/data [cache /home/users/r/rossigng/work/singularity/community-cr-prod.seqera.io-docker-registry-v2-blobs-sha256-34-34e733a9ae16a27e80fe00f863ea1479c96416017f24a907996126283e7ecd4d-data.img]
-[nf-core/demo] Pipeline completed successfully-
WARN: Singularity cache directory has not been defined -- Remote image will be stored in the path: /home/users/r/rossigng/work/singularity -- Use the environment variable NXF_SINGULARITY_CACHEDIR to specify a different location
```
Best regards,


## Post 5 by @Camille.Christe (2026-04-15T11:16:45.753Z)

Thank you ! It is indeed working.
And I see how I can access the script.
