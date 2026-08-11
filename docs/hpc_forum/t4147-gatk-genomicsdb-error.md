# GATK-GenomicsDB error

- Source: https://hpc-community.unige.ch/t/4147

- Created: 2025-11-17T09:43:49.539Z

- Tags: yggdrasil

- Posts: 4

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Elisa.Denis (2025-11-17T09:43:49.609Z)

## Primary informations
Username: $denisel0
Cluster: Yggdrasil
I am using GATK to perform variant calling on my whole genome sequences. Everything is very heavy in terms of time and resources when using the function CombineGVCFs or GenotypeGVCFs so I wanted to switch to GenomicsDBImport which is supposedly more efficient.
However, I always get the same error:
```
Using GATK jar /home/users/d/denisel0/.conda/envs/popgen_gatk/share/gatk4-4.6.2.0-1/gatk-package-4.6.2.0-local.jar
Running:
java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -Xmx32G -Djava.io.tmpdir=/home/users/d/denisel0/scratch/tmp_43348987_9 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -jar /home/users/d/denisel0/.conda/envs/popgen_gatk/share/gatk4-4.6.2.0-1/gatk-package-4.6.2.0-local.jar GenomicsDBImport --genomicsdb-workspace-path /home/users/d/denisel0/scratch/tmp_genomicsdb/9 --sample-name-map /home/users/d/denisel0/scratch/sparganium/bam/rg_fixed/gvcf/gvcf_sample_map.txt -L OZ243046.1:40000001-45000000 --batch-size 50 --reader-threads 8 --overwrite-existing-genomicsdb-workspace true
10:39:02.104 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/home/users/d/denisel0/.conda/envs/popgen_gatk/share/gatk4-4.6.2.0-1/gatk-package-4.6.2.0-local.jar!/com/intel/gkl/native/libgkl_compression.so
10:39:02.309 INFO  GenomicsDBImport - ------------------------------------------------------------
10:39:02.312 INFO  GenomicsDBImport - The Genome Analysis Toolkit (GATK) v4.6.2.0
10:39:02.312 INFO  GenomicsDBImport - For support and documentation go to https://software.broadinstitute.org/gatk/
10:39:02.312 INFO  GenomicsDBImport - Executing as denisel0@cpu077.yggdrasil on Linux v5.14.0-503.40.1.el9_5.x86_64 amd64
10:39:02.312 INFO  GenomicsDBImport - Java runtime: OpenJDK 64-Bit Server VM v17.0.11-internal+0-adhoc..src
10:39:02.312 INFO  GenomicsDBImport - Start Date/Time: November 17, 2025 at 10:39:01 AM CET
10:39:02.312 INFO  GenomicsDBImport - ------------------------------------------------------------
10:39:02.313 INFO  GenomicsDBImport - ------------------------------------------------------------
10:39:02.314 INFO  GenomicsDBImport - HTSJDK Version: 4.2.0
10:39:02.314 INFO  GenomicsDBImport - Picard Version: 3.4.0
10:39:02.314 INFO  GenomicsDBImport - Built for Spark Version: 3.5.0
10:39:02.316 INFO  GenomicsDBImport - HTSJDK Defaults.COMPRESSION_LEVEL : 2
10:39:02.317 INFO  GenomicsDBImport - HTSJDK Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS : false
10:39:02.317 INFO  GenomicsDBImport - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS : true
10:39:02.317 INFO  GenomicsDBImport - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE : false
10:39:02.317 INFO  GenomicsDBImport - Deflater: IntelDeflater
10:39:02.317 INFO  GenomicsDBImport - Inflater: IntelInflater
10:39:02.317 INFO  GenomicsDBImport - GCS max retries/reopens: 20
10:39:02.317 INFO  GenomicsDBImport - Requester pays: disabled
10:39:02.317 INFO  GenomicsDBImport - Initializing engine
10:39:02.519 INFO  IntervalArgumentCollection - Processing 5000000 bp from intervals
10:39:02.521 INFO  GenomicsDBImport - Done initializing engine
10:39:02.746 INFO  GenomicsDBLibLoader - GenomicsDB native library version : 1.5.5-addd795
10:39:02.751 INFO  GenomicsDBImport - Shutting down engine
\[November 17, 2025 at 10:39:02 AM CET\] org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBImport done. Elapsed time: 0.01 minutes.
Runtime.totalMemory()=167772160

---

**A USER ERROR has occurred: Error creating GenomicsDB workspace: /home/users/d/denisel0/scratch/tmp_genomicsdb/9**

---

org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBImport$UnableToCreateGenomicsDBWorkspace: Error creating GenomicsDB workspace: /home/users/d/denisel0/scratch/tmp_genomicsdb/9
at org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBImport.overwriteCreateOrCheckWorkspace(GenomicsDBImport.java:1049)
at org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBImport.onTraversalStart(GenomicsDBImport.java:680)
at org.broadinstitute.hellbender.engine.GATKTool.doWork(GATKTool.java:1117)
at org.broadinstitute.hellbender.cmdline.CommandLineProgram.runTool(CommandLineProgram.java:150)
at org.broadinstitute.hellbender.cmdline.CommandLineProgram.instanceMainPostParseArgs(CommandLineProgram.java:203)
at org.broadinstitute.hellbender.cmdline.CommandLineProgram.instanceMain(CommandLineProgram.java:222)
at org.broadinstitute.hellbender.Main.runCommandLineProgram(Main.java:166)
at org.broadinstitute.hellbender.Main.mainEntry(Main.java:209)
at org.broadinstitute.hellbender.Main.main(Main.java:306)
```
I wanted to know if it could potentially be an issue from yggdrasil which is not able to support this.
Thanks in advance
Kind regards
Elisa Denis


## Post 2 by @Yann.Sagon (2025-11-17T13:10:56.955Z)

Dear @Elisa.Denis[@Elisa.Denis](https://hpc-community.unige.ch/u/elisa.denis),
What I understand: you need a space to host temporary files and you are using scratch for for purpose.
Can you please test with local scratch, it is way faster and is a real POSIX filesystem.
doc.eresearch.unige.ch[doc.eresearch.unige.ch](https://doc.eresearch.unige.ch/hpc/storage_on_hpc#local_storage)
### hpc:storage_on_hpc [eResearch Doc][hpc:storage_on_hpc [eResearch Doc]](https://doc.eresearch.unige.ch/hpc/storage_on_hpc#local_storage)


## Post 3 by @Elisa.Denis (2025-11-17T13:55:37.443Z)

Thanks Yann,
I don’t have access to the local scratch, do I need my PI to write to you directly?
Thanks
Elisa


## Post 4 by @Yann.Sagon (2025-11-18T09:36:36.626Z)

Hi @Elisa.Denis[@Elisa.Denis](https://hpc-community.unige.ch/u/elisa.denis)
the local scratch is only available on a compute node during the job execution, thus it is normal that you don’t see it from the login node.
ps: the forum is working again, thanks for the notification.
Best regards
