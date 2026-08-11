# Cannot compile (Palabos) with apptainer in Baobab or Yggdrasil

- Source: https://hpc-community.unige.ch/t/3283

- Created: 2024-02-01T15:07:55.333Z

- Posts: 6

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Francesco.Marson (2024-02-01T15:07:55.391Z)

Hello HPC team,
I tried to build several NVIDIA HPC SDK Docker images with Apptainer and then tried to compile in a computing node with (Ampere) GPUs. Even if I can compile using the same container both in the laptop and the working station, the same procedure leads to compilation errors both in Baobab and Yggdrasil independently on the version of the NVIDIA HPC SDK container.
what did you try:
```
salloc --partition=shared-gpu --gres=gpu:ampere --time=30:00
apptainer build nvhpc-23.1-devel.sif docker://nvcr.io/nvidia/nvhpc:23.11-devel-cuda12.3-ubuntu22.04
apptainer shell --nv ~/projects/nvhpc-23.1-devel.sif
```
→ build normally as I do in my PC using cmake
what didn’t work:
compilation errors without any information, the compilation errors are different for different containers and always unclear, for example:
```
Building CXX object CMakeFiles/palabos.dir/home/marson/projects/palabos-hybrid-deployed/src/io/multiBlockWriter2D.cpp.o
c++: fatal error: Killed signal terminated program cc1plus
compilation terminated.
make[2]: *** [CMakeFiles/palabos.dir/build.make:132: CMakeFiles/palabos.dir/home/marson/projects/palabos-hybrid-deployed/src/atomicBlock/atomicBlock3D.cpp.o] Error 1
make[2]: *** Waiting for unfinished jobs....
c++: fatal error: Killed signal terminated program cc1plus
compilation terminated.
make[2]: *** [CMakeFiles/palabos.dir/build.make:160: CMakeFiles/palabos.dir/home/marson/projects/palabos-hybrid-deployed/src/atomicBlock/atomicBlockOperations3D.cpp.o] Error 1
c++: fatal error: Killed signal terminated program cc1plus
compilation terminated.
make[2]: *** [CMakeFiles/palabos.dir/build.make:146: CMakeFiles/palabos.dir/home/marson/projects/palabos-hybrid-deployed/src/atomicBlock/atomicBlockOperations2D.cpp.o] Error 1
```
What am I doing wrong?
Thanks for the help,
Francesco
admin edit: code block


## Post 2 by @Yann.Sagon (2024-02-01T16:37:32.544Z)

Hi,
I checked quickly the status of your latest jobs on Baobab: one terminated due to timeout and the other due to out of memory.
You should maybe try to request more memory and or cores? By default you get 3GB of RAM per core.
If you still have the issue, please let us know the job id of the faulty job.
Best


## Post 3 by @Francesco.Marson (2024-02-01T17:39:31.322Z)

Hi Yann,
Indeed, I was running also without enough memory, I increased to --mem=30GB and now I am getting the same errors I was getting the other times I tried
```
ptxas fatal   : Could not open output file '/scratch/tmpxft_000369b3_0000000a'
NVC++-F-0155-Compiler failed to translate accelerator region (see -Minfo messages): Device compiler exited with error status code (/home/marson/projects/palabos-hybrid-deployed/src/algorithm/empiricalData.cpp: 416)
NVC++/x86-64 Linux 23.11-0: compilation aborted
make[2]: *** [CMakeFiles/palabos.dir/build.make:90: CMakeFiles/palabos.dir/home/marson/projects/palabos-hybrid-deployed/src/algorithm/empiricalData.cpp.o] Error 2
make[2]: *** Waiting for unfinished jobs....
```
any clue about the reason for that?
Thank you!!
Francesco


## Post 4 by @Yann.Sagon (2024-02-02T07:41:18.016Z)

Francesco.Marson:
> ptxas fatal : Could not open output file ‘/scratch/tmpxft_000369b3_0000000a’
Now we can see something interesting! As you are compiling from inside a container and you try to write outside your home directory, you must mount bind explicitly `/scratch` when launch the container.
ps: please take some time to edit your previous post to format it properly using code block, it is easier to read.
Best
Yann


## Post 5 by @Francesco.Marson (2024-02-02T10:25:58.341Z)

Hi Yann,
it worked like a charm, thanks a lot. I am just a bit confused about why it is using the “/scratch” folder for the compilation, could you give me an insight on how does this work, so next time I will be more careful about that.
Thank you for your help :slight_smile: ,
Francesco


## Post 6 by @Yann.Sagon (2024-02-02T13:04:49.570Z)

Hi,
Francesco.Marson:
> could you give me an insight on how does this work
This is because we set the variable TMPDIR:
```
(baobab)-[sagon@login2 ~]$ srun env | grep -i tmp
srun: job 7151716 queued and waiting for resources
srun: job 7151716 has been allocated resources
TMPDIR=/scratch
```
