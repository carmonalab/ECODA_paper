# Issue with matlab on Baobab

- Source: https://hpc-community.unige.ch/t/3963

- Created: 2025-06-02T06:43:23.158Z

- Tags: baobab

- Posts: 8

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @nicolas.clairis (2025-06-02T06:43:23.221Z)

I have a matlab script (using the SPM toolbox) that normally works well but suddenly when I launched it in Baobab on the shared-cpu partition, I get the following error messages as if the cluster didn’t manage to open matlab:
```
std::exception: Error loading /unige/matlab2021b/bin/glnxa64/builtins/lxe/lxe/mwm_lxe_builtinimpl.so. /unige/matlab2021b/bin/glnxa64/builtins/lxe/lxe/mwm_lxe_builtinimpl.so: cannot open shared object file: Remote I/O error: Remote I/O error: Remote I/O error

srun: error: cpu329: task 0: Exited with exit code 1
```
or for another subject:
```
std::exception: libmwkeybrd_impl.so: cannot open shared object file: No such file or directory: No such file or directory

srun: error: cpu331: task 0: Exited with exit code 1
```
or for another subject:
```
Caught "std::exception" Exception message is:
Error loading /unige/matlab2021b/bin/glnxa64/matlab_startup_plugins/addons/registry/mwaddons_registry_startup_plugin.so. /unige/matlab2021b/bin/glnxa64/matlab_startup_plugins/addons/registry/mwaddons_registry_startup_plugin.so: cannot open shared object file: Remote I/O error: Remote I/O error: Remote I/O error
Caught "std::exception" Exception message is:
Error loading /unige/matlab2021b/bin/glnxa64/matlab_startup_plugins/connector/worker/libmwconnectorworkerstartupbundle.so. /unige/matlab2021b/bin/glnxa64/matlab_startup_plugins/connector/worker/libmwconnectorworkerstartupbundle.so: cannot open shared object file: Remote I/O error: Connection timed out: Connection timed out
Caught "std::exception" Exception message is:
Error loading /unige/matlab2021b/bin/glnxa64/matlab_startup_plugins/graphics/graphics_init_plugin/mwgraphics_init_plugin.so. /unige/matlab2021b/bin/glnxa64/matlab_startup_plugins/graphics/graphics_init_plugin/mwgraphics_init_plugin.so: cannot open shared object file: Remote I/O error: Connection timed out: Connection timed out
Caught "std::exception" Exception message is:
Error loading /unige/matlab2021b/bin/glnxa64/matlab_startup_plugins/startup_plugins/startup_folder/mwStartupFolderPlugin.so. /unige/matlab2021b/bin/glnxa64/matlab_startup_plugins/startup_plugins/startup_folder/mwStartupFolderPlugin.so: cannot open shared object file: Remote I/O error: Connection timed out: Connection timed out
Caught "std::exception" Exception message is:
Message Catalog MATLAB:legacy_two_part was not loaded from the file. Please check file location, format or contents
Caught "std::exception" Exception message is:
Message Catalog MATLAB:legacy_two_part was not loaded from the file. Please check file location, format or contents
Caught "std::exception" Exception message is:
Error loading /unige/matlab2021b/bin/glnxa64/builtins/matlab_toolbox_lang_builtins/mwdiagnostic_builtinimpl.so. /unige/matlab2021b/bin/glnxa64/builtins/matlab_toolbox_lang_builtins/mwdiagnostic_builtinimpl.so: cannot open shared object file: Remote I/O error: Invalid argument: Invalid argument
```
Can you explain me what is going on?


## Post 2 by @nicolas.clairis (2025-06-02T07:19:39.548Z)

The error occurred when I launched 3 independent scripts in parallel with each launching a parallel over 74 subjects (so I ask for 74*3=222 nodes and matlab licences in total). Weirdly, when I launch the first one alone, at least for now, there is no error message. Could the error come from the fact that I ask for too many matlab licences in parallel? I am afraid to launch the 2nd script now because I don’t want to mess up with the fact that the first one runs well… Any idea about what may have happened there so that I know if I can launch my 2 other scripts or if I need to wait for the first to finish before doing so?


## Post 3 by @nicolas.clairis (2025-06-02T08:06:00.419Z)

FYI: I still tried to launch the 2nd script. The first task is still running without any issue, but my 2nd task immediately crashed with the following error code:
```
std::exception: libmwkeybrd_impl.so: cannot open shared object file: No such file or directory: No such file or directory

srun: error: cpu325: task 0: Exited with exit code 1
```
Can you explain me why? What is libmwkeybrd_impl and why is it called and not accessible? This is not an error related to my script but to the cluster from what I understand


## Post 4 by @nicolas.clairis (2025-06-02T10:35:50.873Z)

As I said before, when launching the scripts sequentially (launching my script launch_preprocessing_T1.sbatch only once launch_preprocessing_T0.sbatch has finished running rather than at the same time) the issue seems to be solved and the problem does not emerge anymore, although they are completely independent in principle. I’m guessing that the problem was coming from the fact that too many matlab licences were asked at the same time(?) In any case, if that’s really the issue it’s a bit annoying as it kills the whole interest to use the cluster and parallelize but at least my issue is solved and the scripts are running right now :slight_smile:


## Post 5 by @Yann.Sagon (2025-06-02T14:42:19.438Z)

Dear @nicolas.clairis[@nicolas.clairis](https://hpc-community.unige.ch/u/nicolas.clairis)
first things not obvious: we have two Matlab installation: one hosted on the cluster and one hosted on an external machine. The one you are using is the legacy version on an external machine. We’ll remove it as there is now not reason to keep it.
```
(baobab)-[sagon@login1 ~]$ ml spider MATLAB

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  MATLAB:
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
     Versions:
        MATLAB/2021a
        MATLAB/2021b << this version is hosted outside of Baobab >>
        MATLAB/2022a
```
Can you please try using the version `MATLAB/2022a`?
About your scripts: please share them here so we can give advises.


## Post 6 by @nicolas.clairis (2025-06-02T15:05:59.715Z)

Thanks for the reply, but I never specified to use that old version so not sure why the issue happened. You can find my script below, as far as I can understand there is no direct call to any specific version of Matlab so I guess that the script just tried to use any matlab licence that was available and when the more recent was not available it used the older one which created the bug. Could that explain my former issue? In any case, the scripts worked now by launching them sequentially
```
#!/bin/bash
#SBATCH --job-name=wololo_preprocessing

# Send an email when the job is completed
#SBATCH --mail-user=nicolas.clairis@unige.ch
#SBATCH --mail-type=END

# Select partition
#SBATCH --partition=shared-cpu
#SBATCH --mem-per-cpu=64G # define RAM to use for each node

# Request CPU resource for a serial job
#SBATCH --array=1-74 # This will tell the batch that there are 74 subjects (i.e. tasks) to process
#SBATCH --ntasks-per-node=1 # number of tasks (i.e. subjects) per node

# Set your minimum acceptable walltime, format: day-hours:minutes:seconds
#SBATCH --time=12:00:00

# Set the working directory - put your MATLAB script there
#SBATCH --chdir=/home/users/c/clairis/scripts/
#SBATCH --output=logs/preprocessing/preprocessing_matlab_output_T0_%j_%a.txt 
#SBATCH --error=logs/preprocessing/preprocessing_matlab_error_T0_%j_%a.txt 

# load relevant softwares on Baobab (Matlab)
module load MATLAB

# load list of subjects to preprocess (the launch_preprocessing batch makes the loop inside Matlab, but launch_preprocessing2 will make it inside the batch directly)
SUBJECTS=("E002" "E003" "E004" "E006" "E007" "E010" "E011" "E012" "E014" "E015" "E017" "E019" "E020" "E021" "E022" "E023" "E028" "E035" "E036" "E037" "E043" "E044" "E048" "E050" "E052" "E055" "E056" "E059" "E061" "E063" "E064" "E068" "E070" "E076" "E078" "E084" "E090" "E092" "E093" "E100" "E101" "E107" "E108" "E109" "E111" "E112" "E114" "E115" "E117" "E123" "E124" "E125" "E126" "E128" "E130" "E131" "E133" "E135" "E136" "E137" "E138" "E139" "E140" "E141" "E142" "E145" "E146" "E148" "E149" "E150" "E151" "E153" "E158" "E161");

# perform the loop
SUBJECT=${SUBJECTS[$SLURM_ARRAY_TASK_ID-1]} # $SLURM_ARRAY_TASK_ID retrieves the job array index. We subtract 1 ($SLURM_ARRAY_TASK_ID-1) because Bash arrays are zero-indexed while Slurm start at 1 so that the first subject would be skipped if not doing so)

# Afficher le participant en cours (debugging)
echo "Processing subject $SUBJECT on node $SLURMD_NODENAME"

# lancer le script
srun matlab -nodesktop -softwareopengl -nosplash -nodisplay -r "addpath('/home/users/c/clairis/scripts'); preprocessing_loop_T0('$SUBJECT'); exit;" # adds script folder path + launches the script
```


## Post 7 by @Yann.Sagon (2025-06-02T15:20:03.360Z)

nicolas.clairis:
> `module load MATLAB`
It is a good practice to always specify the version of the soft you use, this is more reproducible. And by default this is version 2021b that is loaded.
```
module load MATLAB/2022a
```
You can see which one is the default version:
```
(baobab)-[sagon@login1 ~]$ ml avail MATLAB

---------------------------------------------------------------------------------------------------------- /opt/modulefiles/Linux ----------------------------------------------------------------------------------------------------------
   MATLAB/2021b (D)

--------------------------------------------------------------------------------------------------------- /opt/ebmodules/all/Core ----------------------------------------------------------------------------------------------------------
   MATLAB/2021a    MATLAB/2022a

  Where:
   D:  Default Module
```
Anyway, as I said in my previous post, we should get ride of this version.
In your sbatch script, you are requesting 64GB per job. Each job has only one CPU. This means that on a 64GB compute node with 16 cores, you can only run one job as the memory will be exhauted. Are you sure you need that much memory per job? If yes, can you add more CPUs per task to have a faster processing?
Each time you run a matlab instance, it consumes a license. You need to constraint the number of parallel instances by requesting a license matlab as explained here: hpc:applications_and_libraries [eResearch Doc][hpc:applications_and_libraries [eResearch Doc]](https://doc.eresearch.unige.ch/hpc/applications_and_libraries#matlab)
A good practice is to compile your matlab code. Doing so, you don’t need to use a license when running your job. hpc:applications_and_libraries [eResearch Doc][hpc:applications_and_libraries [eResearch Doc]](https://doc.eresearch.unige.ch/hpc/applications_and_libraries#compile_your_matlab_code)


## Post 8 by @nicolas.clairis (2025-06-02T15:58:48.699Z)

ok thanks for the tips. I will force the scripts to use Matlab 2022a from now on.
