# Problem with matlab compilator mcc on bamboo

- Source: https://hpc-community.unige.ch/t/3763

- Created: 2024-12-09T17:22:59.537Z

- Tags: bamboo

- Posts: 4

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Nicolas.Cuny (2024-12-09T17:22:59.577Z)

## Primary informations
Username: cuny
Cluster: bamboo
I try to compile a matlab code in an executable with mcc on bamboo but get the following error:
```
mclCopyDir failed: Failed copying '/home/users/c/cuny/.matlab/R2022a/matlab.prf' to '/tmp/mcc1340-7eef-fa78-08e6.tmp/user/matlab.prf':  reason 'fl:filesystem:SystemError'
Failed to copy user-level MATLAB settings from "/home/users/c/cuny/.matlab/R2022a" into temporary directory "/tmp/mcc1340-7eef-fa78-08e6.tmp/user". Cannot continue without risk of corrupting MATLAB settings.
```
Here is the command I use:
```
mcc -m -v -R '-nojvm, -nodisplay' -I functions_dynamics/ -o epithelia_closed Epithelia_Dynamics.m
```


## Post 2 by @Yann.Sagon (2024-12-10T08:18:41.934Z)

Dear @Nicolas.Cuny[@Nicolas.Cuny](https://hpc-community.unige.ch/u/nicolas.cuny)
on which server do you issue this command? login1 or a compute node?


## Post 3 by @Nicolas.Cuny (2024-12-10T17:29:13.973Z)

I did try it both on login1 and on compute node and did get the same result.


## Post 4 by @Yann.Sagon (2024-12-11T15:08:58.756Z)

I was able to reproduce the error. I tried like that and I went a little bit futher:
```
(bamboo)-[cuny@login1 TestYS]$ mkdir my_pref
(bamboo)-[cuny@login1 TestYS]$ MATLAB_PREFDIR=$PWD/my_pref mcc -m -v -R '-nojvm, -nodisplay' -I functions_dynamics/ -o epithelia_closed Epithelia_Dynamics.m
Opening log file:  /home/users/c/cuny/java.log.2348
Compiler version: 8.4 (R2022a)
Analyzing file dependencies.
Parsing file "/home/users/c/cuny/EpithelialShell/TestYS/Epithelia_Dynamics.m"
        (referenced from command line).
Failed to build CTF file: 'epithelia_closed.ctf'. Details: Error in closing CTF file 'epithelia_closed.ctf'. Details: 'Could not commit changes:  copy failed: fl:filesystem:SystemError: /tmp/.f1e3-1657-7cd0-11d3.tmp Invalid cross-device link'
```
Then, I tried again overriding `$TMP`. You have to empty the `my_pref` directory first.
```
(bamboo)-[cuny@login1 TestYS]$ rm my_pref/*
(bamboo)-[cuny@login1 TestYS]$ TMP=$PWD/my_tmp MATLAB_PREFDIR=$PWD/my_pref mcc -m -v -R '-nojvm, -nodisplay' -I functions_dynamics/ -o epithelia_closed Epithelia_Dynamics.m
Opening log file:  /home/users/c/cuny/java.log.23732
Compiler version: 8.4 (R2022a)
Analyzing file dependencies.
Parsing file "/home/users/c/cuny/EpithelialShell/TestYS/Epithelia_Dynamics.m"
        (referenced from command line).
Generating file "/home/users/c/cuny/EpithelialShell/TestYS/readme.txt".
Generating file "run_epithelia_closed.sh".
```
Seems to work?!
