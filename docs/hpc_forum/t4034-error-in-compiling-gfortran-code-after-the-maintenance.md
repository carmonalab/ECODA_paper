# Error in compiling gfortran code after the maintenance

- Source: https://hpc-community.unige.ch/t/4034

- Created: 2025-08-06T13:00:44.039Z

- Tags: yggdrasil

- Posts: 5

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Arianna.Nigioni (2025-08-06T13:00:44.089Z)

## Primary informations
Username: nigioni
Cluster: Yggdrasil
## Description
I went in the source folder of my code (path: /home/users/n/nigioni/arianna_vJuly2024/sources). In this folder there is a Makefile which is responsible for the compiling of all the fortran files within the same folder. I typed in the terminal: make clean, then module load GCC and finally make -j
The make -j command, instead of generating the executable, gave this error:
---
gfortran: fatal error: cannot execute ‘/opt/ebsofts/GCCcore/14.2.0/libexec/gcc/x86_64-pc-linux-gnu/14.2.0/f951’: posix_spawn: Operation not permitted
compilation terminated.
gfortran: fatal error: cannot execute ‘/opt/ebsofts/GCCcore/14.2.0/libexec/gcc/x86_64-pc-linux-gnu/14.2.0/f951’: posix_spawn: Operation not permitted
gfortran: fatal error: cannot execute ‘/opt/ebsofts/GCCcore/14.2.0/libexec/gcc/x86_64-pc-linux-gnu/14.2.0/f951’: posix_spawn: Operation not permitted
compilation terminated.
compilation terminated.
make: *** [Makefile:211: structure_attached.o] Error 1
make: *** Waiting for unfinished jobs…
make: *** [Makefile:211: RK.o] Error 1
make: *** [Makefile:211: structure.o] Error 1
gfortran: fatal error: cannot execute ‘/opt/ebsofts/GCCcore/14.2.0/libexec/gcc/x86_64-pc-linux-gnu/14.2.0/f951’: posix_spawn: Operation not permitted
compilation terminated.
make: *** [Makefile:211: structure_detached.o] Error 1
---
I never had this problem before when compiling and I did not modify any file going from before to after the last maintenance. If it is needed, I can provide the contents of the Makefile.
Thank you for your help,
Arianna


## Post 2 by @Adrien.Albert (2025-08-07T13:04:50.556Z)

Hi @Arianna.Nigioni[@Arianna.Nigioni](https://hpc-community.unige.ch/u/arianna.nigioni)
I have copied your source in my home
```
(yggdrasil)-[alberta@login1 ~]$ salloc --partition=shared-cpu,shared-bigmem --time=00:10:00 -c 16 --mem=32G
salloc: Pending job allocation 41935333
salloc: job 41935333 queued and waiting for resources
salloc: job 41935333 has been allocated resources
salloc: Granted job allocation 41935333
salloc: Waiting for resource configuration
salloc: Nodes cpu156 are ready for job
(yggdrasil)-[alberta@cpu156 ~]$ pwd
/home/users/a/alberta
(yggdrasil)-[alberta@cpu156 ~]$ cd arianna_vJuly2024/sources/
(yggdrasil)-[alberta@cpu156 sources]$ ml GCCcore/14.1.0 
(yggdrasil)-[alberta@cpu156 sources]$ make clean; make -j
[...]
(yggdrasil)-[alberta@cpu156 sources]$ echo $?
0
```
I have no errors, could you tell me on which node you try to compile your software ?


## Post 3 by @Adrien.Albert (2025-08-07T13:24:15.961Z)

I have tested as a evil user :alien_monster: on the login node and encountered the same error. After checking the user limitations, it appears that you have exceeded the task limit (150), which caused the compilation to fail, this is expected behavior.
```
(yggdrasil)-[alberta@login1 ~]$ systemctl status user-401775.slice 
Failed to fork: Resource temporarily unavailable
● user-401775.slice - User Slice of UID 401775
     Loaded: loaded
    Drop-In: /usr/lib/systemd/system/user-.slice.d
             └─10-defaults.conf
             /etc/systemd/system.control/user-401775.slice.d
             └─50-BlockIOWeight.conf, 50-CPUQuota.conf, 50-MemoryMax.conf, 50-TasksMax.conf
     Active: active since Thu 2025-08-07 14:57:37 CEST; 8min ago
      Until: Thu 2025-08-07 14:57:37 CEST; 8min ago
       Docs: man:user@.service(5)
      Tasks: 150 (limit: 150)
     Memory: 450.6M (max: 8.0G available: 7.5G)
        CPU: 1min 346ms
     CGroup: /user.slice/user-401775.slice
             ├─session-33818.scope
             │ ├─2418513 "sshd: alberta [priv]"
             │ ├─2418586 "sshd: alberta@pts/68"
             │ ├─2418587 -bash
             │ ├─2449773 make -j
             │ ├─2449792 gfortran -ffixed-line-length-132 -fimplicit-none -fbacktrace -fdefault-real-8 -fdefault-double-8 -O3 -g -funroll-loops -ftree-vectorize -ftree-loop-optimize -mavx -Wall -Wno-tabs -Wno-u…
             │ ├─2449794 /opt/ebsofts/GCCcore/14.1.0/libexec/gcc/x86_64-pc-linux-gnu/14.1.0/f951 mercury/mercury6_embed.f -ffixed-form -I . -quiet -dumpdir mercury/ -dumpbase mercury6_embed.f -dumpbase-ext .f -…
             │ ├─2449861 gfortran -ffixed-line-length-132 -fimplicit-none -fbacktrace -fdefault-real-8 -fdefault-double-8 -O3 -g -funroll-loops -ftree-vectorize -ftree-loop-optimize -mavx -Wall -Wno-tabs -Wno-u…
             │ ├─2449862 gfortran -ffixed-line-length-132 -fimplicit-none -fbacktrace -fdefault-real-8 -fdefault-double-8 -O3 -g -funroll-loops -ftree-vectorize -ftree-loop-optimize -mavx -Wall -Wno-tabs -Wno-u…
             │ ├─2449863 gfortran -ffree-line-length-none -fimplicit-none -fbacktrace -fdefault-real-8 -fdefault-double-8 -O3 -g -funroll-loops -ftree-vectorize -ftree-loop-optimize -mavx -Wall -Wno-tabs -Wno-u…
             │ ├─2449864 gfortran -ffree-line-length-none -fimplicit-none -fbacktrace -fdefault-real-8 -fdefault-double-8 -O3 -g -funroll-loops -ftree-vectorize -ftree-loop-optimize -mavx -Wall -Wno-tabs -Wno-u…
```
The compilation process attempts to use all available resources, which is why it’s recommended (or rather, required) to use a compute node for such tasks.


## Post 4 by @Arianna.Nigioni (2025-08-07T14:22:01.366Z)

Ok I understand, now it works, thank you !


## Post 5 by @Adrien.Albert (2025-08-12T07:40:58.629Z)

Extra: Best Practice
When compiling, it’s strongly recommended to run:
`make -j <cpu>`
where `<cpu>` is the number of cores actually allocated to your build environment.
This prevents unexpected behavior by avoiding the total core count detected by the system, which might be higher than what you can really use.
doc.eresearch.unige.ch[doc.eresearch.unige.ch](https://doc.eresearch.unige.ch/hpc/hpc_clusters#how_do_our_clusters_work:)
### hpc:hpc_clusters [eResearch Doc][hpc:hpc_clusters [eResearch Doc]](https://doc.eresearch.unige.ch/hpc/hpc_clusters#how_do_our_clusters_work:)
> - a login node (aka headnode) allowing users to connect and submit jobs to the cluster. Aach user is limited to 2 CPU cores and 8 GB of RAM on the login node.
