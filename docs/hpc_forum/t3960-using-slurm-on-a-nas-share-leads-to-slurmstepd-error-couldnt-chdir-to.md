# Using slurm on a NAS share leads to slurmstepd: error: couldn't chdir to?

- Source: https://hpc-community.unige.ch/t/3960

- Created: 2025-05-23T13:34:24.199Z

- Posts: 13

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Matthieu.Stigler (2025-05-23T13:34:24.245Z)

I am trying to use slurm to run an R script hosted on a NASS share point. On the compute node, I am able to mount the NASS drive, `cd` to it and run the script with `R CMD BATCH` (similar to `Rscript`), but running the same command with `srun` in front will fail, and hence will fail within a bash script too?
Here’s what I do:
```
module purge               
dbus-launch bash
gio mount smb://isis.unige.ch/nasac/gsem/stigler-smb/STIGLER < .credentials
cd /run/user/451690/gvfs/smb-share:server=isis.unige.ch,share=nasac/gsem/stigler-smb/shared_projects

module load GCC/12.3.0 OpenMPI/4.1.5 R/4.3.2  R-bundle-CRAN/2023.12
R CMD BATCH /run/user/451690/gvfs/smb-share:server=isis.unige.ch,share=nasac/gsem/stigler-smb/shared_projects/test.R
```
This works. But now the same with srun:
```
srun R CMD BATCH /run/user/451690/gvfs/smb-share:server=isis.unige.ch,share=nasac/gsem/stigler-smb/shared_projects/test.R
```
leads to:
> slurmstepd: error: couldn’t chdir to `/run/user/451690/gvfs/smb-share:server=isis.unige.ch,share=nasac/gsem/stigler-smb/shared_projects': No such file or directory: going to /tmp instead slurmstepd: error: couldn't chdir to `/run/user/451690/gvfs/smb-share:server=isis.unige.ch,share=nasac/gsem/stigler-smb/shared_projects’: Permission denied: going to /tmp instead


## Post 2 by @Yann.Sagon (2025-05-27T15:16:40.582Z)

Dear @Matthieu.Stigler[@Matthieu.Stigler](https://hpc-community.unige.ch/u/matthieu.stigler)
You have to “mount” the NAS share on every host you use. In your first code block you probably executed the command on the login node to have access to the storage. This is fine. But as soon as you submit a job using sbatch or srun, your job is executed in a compute node where the filesystem isn’t mounted, you need to mount it first.
My advice is to create an sbatch[sbatch](https://doc.eresearch.unige.ch/hpc/slurm#batch_mode_sbatch) script instead or using directly srun, in this script you mount the filesystem and execute then your rscript.


## Post 3 by @Matthieu.Stigler (2025-05-27T15:53:48.213Z)

Thanks @Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon)
This is what I am doing already, mounting the NAS share on the compute node, not the login node.
The same problem happens with an sbatch. I just ran 18336658 and got the same error.
Here’s the .sh script:
```
#!/bin/sh
#SBATCH --job-name R_test           # this is a parameter to help you sort your job when listing it
#SBATCH --error slurm_logs/job_%j_%x_error.txt     # optional. By default a file slurm-{jobid}.out will be created
#SBATCH --output slurm_logs/job_%j_%x_out.txt      # optional. By default the error and output files are merged
#SBATCH --ntasks 1                    # number of tasks in your job. One by default
#SBATCH --cpus-per-task 1             # number of cpus for each task. One by default
#SBATCH --partition debug-cpu         # the partition to use. By default debug-cpu
#SBATCH --time 15:00                  # maximum run time.
 
module purge               
dbus-launch bash
gio mount smb://isis.unige.ch/nasac/gsem/XXX < .credentials
cd /run/user/451690/gvfs/smb-share:server=isis.unige.ch,share=nasac/gsem/XXX

module load GCC/12.3.0 OpenMPI/4.1.5 R/4.3.2  R-bundle-CRAN/2023.12
srun R CMD BATCH /run/user/451690/gvfs/smb-share:server=isis.unige.ch,share=nasac/gsem/XXX/test.R
```


## Post 4 by @Yann.Sagon (2025-05-28T08:26:20.926Z)

Where is located the sbatch script please? Do you launch it with `sbatch mybatch.sh` or similar?


## Post 5 by @Matthieu.Stigler (2025-05-28T08:56:52.497Z)

I run the following command, from the login node:
```
sbatch --ntasks=1 --partition=shared-cpu --time=0-00:02:00 bash_test.sh
```
the file `bash_test.sh`, shown above, is in my user home. Its goal is to run `R CMD BATCH` on a script that is on the NAS share.
In case that can help, I just ran the batch job 18349423, which gave in the log the error message:
> slurmstepd: error: couldn’t chdir to `/run/user/451690/gvfs/smb-share:server=isis.unige.ch,share=nasac/gsem/XXX/shared_projects’: Permission denied: going to /tmp instead
Thanks!


## Post 6 by @Adrien.Albert (2025-05-28T10:19:01.275Z)

Hi @Matthieu.Stigler[@Matthieu.Stigler](https://hpc-community.unige.ch/u/matthieu.stigler)
Here are my tests. The script in the SMB share contains only:
```
echo toto
```
### Reproducing the current situation (Not working):
```
(baobab)-[alberta@login1 ~]$ cat !$
cat sbatch_test.sh
#!/bin/sh
#SBATCH --job-name test_gio
#SBATCH --cpus-per-task 1
#SBATCH --time 00:05:00
#SBATCH --partition debug-cpu

dbus-launch bash
gio  mount smb://isis.unige.ch/nasac/hpc_exchange/backup < .credentials
bash /var/run/user/401775/gvfs/smb-share\:server\=isis.unige.ch\,share\=nasac/hpc_exchange/backup/tutu

(baobab)-[alberta@login1 ~]$ sbatch !$
sbatch sbatch_test.sh
Submitted batch job 18350710
    
     
(baobab)-[alberta@login1 ~]$ cat slurm-18350710.out 
Error creating proxy: Could not connect: No such file or directory (g-io-error-quark, 1)
gio: smb://isis.unige.ch/nasac/hpc_exchange/backup: volume doesn’t implement mount
bash: /var/run/user/401775/gvfs/smb-share:server=isis.unige.ch,share=nasac/hpc_exchange/backup/tutu: No such file or directory
```
## Resolution
After analyzing the logs on the compute node, I noticed that `dbus-launch` had not fully completed its initialization when `gio` was executed.
To address this, I added a `sleep 5` after the `dbus-launch` command, and it appears to resolve the issue.
The reason is that `dbus-launch` returns immediately and sets up the D-Bus session environment in the background. This means that `gio` can start running before the D-Bus session is fully ready, which results in an error. By adding a short delay (e.g., `sleep 5`), we give D-Bus enough time to complete its initialization, allowing `gio mount` to work as expected.
```
#!/bin/sh
#SBATCH --job-name test_gio
#SBATCH --cpus-per-task 1
#SBATCH --time 00:05:00
#SBATCH --partition debug-cpu

dbus-launch bash
sleep 5
gio  mount smb://isis.unige.ch/nasac/hpc_exchange/backup < .credentials
bash /var/run/user/401775/gvfs/smb-share\:server\=isis.unige.ch\,share\=nasac/hpc_exchange/backup/tutu
```
```
(baobab)-[alberta@login1 ~]$ !sbat
sbatch sbatch_test.sh
Submitted batch job 18350763
(baobab)-[alberta@login1 ~]$ cat slurm-18350763.out 
Authentication Required
Enter user and password for share “nasac” on “isis.unige.ch”:
User [alberta]: Domain [SAMBA]: Password: 

toto   <------ Result
```


## Post 7 by @Matthieu.Stigler (2025-05-31T10:03:00.293Z)

Thanks, this helps a lot!
I am still facing a complicated issue:
- I can run: `srun R CMD BATCH /run/user/451690/gvfs/smb-share:server=isis.unige.ch,share=nasac/gsem/XXX/test.R`
- I cannot cd to that directory: `cd /run/user/451690/gvfs/smb-share:server=isis.unige.ch,share=nasac/gsem/XXX`
If I try to cd to the NAS share:
- There is an error message in the logs: `slurmstepd: error: couldn't chdir to `/run/user/451690/gvfs/smb-share:server=isis.unige.ch,share=nasac/gsem/XXX: Permission denied: going to /tmp instead`
- It seems though the `cd` call is actually working, as: `output=$(ls); echo "$output"` will print the correct files
- But the call `srun R CMD BATCH test.R` won’t work anymore…
Any idea of what is happening? Why slurm is complaining about this `cd`? The problem is that I need to run `R CMD BATCH` in the NAS folder, and also using `srun --chdir /run/user/451690/gvfs/smb-share:server=isis.unige.ch,share=nasac/gsem/XXX/ R CMD BATCH test.R` won’t work.
Thanks!


## Post 8 by @Yann.Sagon (2025-06-02T14:48:31.089Z)

Matthieu.Stigler:
> There is an error message in the logs: `slurmstepd: error: couldn't chdir to `/run/user/451690/gvfs/smb-share:server=isis.unige.ch,share=nasac/gsem/XXX: Permission denied: going to /tmp instead`
Your are probably running the sbatch from inside the nas share. As this path doesn’t exist yet on the compute node, this won’t work until the share is mounted.
Matthieu.Stigler:
> output=$(ls); echo “$output”
Or you can simply write `ls`. The cd step is probably working.
Matthieu.Stigler:
> But the call `srun R CMD BATCH test.R` won’t work anymore…
This I don’t understand, it should work. What is the error message of this step? Maybe try without the `srun` in front of the command?
Best regards


## Post 9 by @Matthieu.Stigler (2025-06-03T10:04:43.955Z)

Thanks @Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon) for your message. I don’t understand what you mean by  “Your are probably running the sbatch from inside the nas share.”. I am running the sbatch file from the login node, and as you see below (script printed) my sbatch file does the `gio mount` NAS first, then `cd` to NAS, then `srun R CMD BATCH` on NAS.
I wanted to try your suggestion to run wihtout `srun`, but baobab is under maintenance. I tried then on bamboo, and get there a different error message running the same script:
> Error creating proxy: Could not connect: Connection refused (g-io-error-quark, 39)
> gio: smb://isis.unige.ch/nasac/gsem/stigler-smb/STIGLER: volume doesn’t implement mount
Importantly, I am running exactly the same steps, and the command `gio mount` will work when conducted outside of sbatch?
```
#!/bin/sh
#SBATCH --job-name R_test           # this is a parameter to help you sort your job when listing it
#SBATCH --error slurm_logs/job_%j_%x_error.txt     # optional. By default a file slurm-{jobid}.out will be created
#SBATCH --output slurm_logs/job_%j_%x_out.txt      # optional. By default the error and output files are merged
#SBATCH --ntasks 1                    # number of tasks in your job. One by default
#SBATCH --cpus-per-task 1             # number of cpus for each task. One by default
#SBATCH --partition debug-cpu         # the partition to use. By default debug-cpu
#SBATCH --time 15:00                  # maximum run time.
 
module purge             
dbus-launch bash
sleep 5
gio mount smb://isis.unige.ch/nasac/gsem/XXX < .credentials

#sleep 1
cd /run/user/451690/gvfs/smb-share:server=isis.unige.ch,share=nasac/gsem/XXX
ls
#output=$(ls)
#echo "$output"

module load GCC/12.3.0 OpenMPI/4.1.5 R/4.3.2  R-bundle-CRAN/2023.12
#srun R CMD BATCH /run/user/451690/gvfs/smb-share:server=isis.unige.ch,share=nasac/gsem/XXX/test.R
##srun --chdir /run/user/451690/gvfs/smb-share:server=isis.unige.ch,share=nasac/gsem/XXX R CMD BATCH test.R
R CMD BATCH test.R

echo "Done!!"
```


## Post 10 by @Yann.Sagon (2025-06-11T12:40:05.387Z)

Matthieu.Stigler:
> I don’t understand what you mean by “Your are probably running the sbatch from inside the nas share.”
I meant that maybe you did something as `cd /run/user/451690/gvfs/smb-share:server=isis.unige.ch,share=nasac/gsem/XXX` and then `sbatch myscript.sh`. As this path isn’t available on the compute node before it is mounted, this will produce an error such as `slurmstepd: error: couldn't chdir to /run/user/451690/gvfs/smb-share:server=isis.unige.ch,share=nasac/gsem/XXX: Permission denied: going to /tmp instead`
Can you confirm this is the case? If yes, try to do a `cd somewherelse` and do the sbatch from there.


## Post 11 by @Matthieu.Stigler (2025-06-16T10:05:29.580Z)

Hi Yann
No, I run the bash script from the login node, at the root level, and only the bash script does the cd part, as you can see in my previous message containing the script.
It seems at that point the easiest would be if you could do some testing, as Adrien started? You would just need to do:
- Create `test.R` on the NASS share: `echo 'pdf("test.pdf"); plot(1); dev.off()' > test.R`
- Run the bash script above, adjusting the path?
Also, can you comment on the error on bamboo versus baobab? The `gio mount` step is working on baobab, but not bamboo, with error:
> Error creating proxy: Could not connect: Connection refused (g-io-error-quark, 39)
> gio: smb://isis.unige.ch/nasac/gsem/stigler-smb/STIGLER: volume doesn’t implement mount
Thanks a lot!


## Post 12 by @Yann.Sagon (2025-06-19T15:08:14.064Z)

Matthieu.Stigler:
> Also, can you comment on the error on bamboo versus baobab?
Is this still the case after the Bamboo maintenance this week?


## Post 13 by @Matthieu.Stigler (2025-06-23T14:23:57.309Z)

@Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon) this is indeed not the case anymore on bamboo, I don’t get the error message. But I still get the same error, though it seems they are spurious warning messages now: the child script is well executed.
