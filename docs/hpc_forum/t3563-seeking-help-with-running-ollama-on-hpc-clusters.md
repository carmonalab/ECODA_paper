# Seeking Help with Running Ollama on HPC Clusters

- Source: https://hpc-community.unige.ch/t/3563

- Created: 2024-07-24T09:15:13.540Z

- Tags: baobab

- Posts: 8

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Jamil.Zaghir (2024-07-24T09:15:13.592Z)

## Primary informations
Username: zaghir
Cluster: baobab
## Description
Hello everyone,
I’m currently working with HPC clusters and I’m facing some challenges with running Ollama (https://ollama.com/[https://ollama.com/](https://ollama.com/)), which involves instantiating an AI model (e.g. Llama3) and performing inferences using GPUs.
The primary issue is that the installation of Ollama requires both sudo privileges and service-related permissions with systemctl. Obtaining such permissions is typically problematic (Adrien from the HPC team tried to find a way to install it through Easybuild, but in vain).
We attempted to use Apptainer as a workaround, but we encountered issues related to the TMPDIR configuration. Our TMPDIR is set to /scratch instead of /tmp, and Ollama does not seem to function properly with this setup.
```
Apptainer> ollama serve &
[1] 2497999
Apptainer> Couldn't find '/home/users/z/zaghir/.ollama/id_ed25519'. Generating new private key.
Your new public key is:

ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIGQnzA7pbOXXkts3USdarkYZgJ862sfht0v/NxRBiIfH

2024/07/25 10:27:09 routes.go:1100: INFO server config env="map[CUDA_VISIBLE_DEVICES: GPU_DEVICE_ORDINAL: HIP_VISIBLE_DEVICES: HSA_OVERRIDE_GFX_VERSION: OLLAMA_DEBUG:false OLLAMA_FLASH_ATTENTION:false OLLAMA_HOST:http://0.0.0.0:11434 OLLAMA_INTEL_GPU:false OLLAMA_KEEP_ALIVE:5m0s OLLAMA_LLM_LIBRARY: OLLAMA_MAX_LOADED_MODELS:0 OLLAMA_MAX_QUEUE:512 OLLAMA_MODELS:/home/users/z/zaghir/.ollama/models OLLAMA_NOHISTORY:false OLLAMA_NOPRUNE:false OLLAMA_NUM_PARALLEL:0 OLLAMA_ORIGINS:[http://localhost https://localhost http://localhost:* https://localhost:* http://127.0.0.1 https://127.0.0.1 http://127.0.0.1:* https://127.0.0.1:* http://0.0.0.0 https://0.0.0.0 http://0.0.0.0:* https://0.0.0.0:* app://* file://* tauri://*] OLLAMA_RUNNERS_DIR: OLLAMA_SCHED_SPREAD:false OLLAMA_TMPDIR: ROCR_VISIBLE_DEVICES:]"
time=2024-07-25T10:27:09.612+02:00 level=INFO source=images.go:784 msg="total blobs: 0"
time=2024-07-25T10:27:09.613+02:00 level=INFO source=images.go:791 msg="total unused blobs removed: 0"
time=2024-07-25T10:27:09.617+02:00 level=INFO source=routes.go:1147 msg="Listening on 0.0.0.0:11434 (version 0.2.8)"
Error: unable to initialize llm library failed to generate tmp dir: stat /scratch: no such file or directory
```
Has anyone successfully run Ollama on HPC clusters, or does anyone have suggestions for other efficient implementations that might be more compatible with HPC environments?
For those unfamiliar with Ollama, it works the following way:
Step 1- instantiate the serv (if Ollama is installed):
`ollama serve`
or daemon version: `ollama serve &`
Step 2- Instantiate and make Llama3 expect prompt through CLI:
`ollama run llama3`
(I usually do Step 2 with python integration using GitHub - ollama/ollama-python: Ollama Python library[GitHub - ollama/ollama-python: Ollama Python library](https://github.com/ollama/ollama-python), which is already installed. But the Step 1 still requires Ollama).
Any guidance or advice would be greatly appreciated.
Thank you!
EDIT: Added outputted error from the Apptainer solution.


## Post 2 by @Jamil.Zaghir (2024-07-25T15:26:40.058Z)

I found a way to deal with the tmp dir. I needed to customize my .sif file by adding a symbolic link between scratch and tmp.
- Make a Docker def file, let’s name it ollama.def:
```
Bootstrap: docker
From: ollama/ollama:0.2.8

%post
    # Create the symbolic link
    ln -s /tmp /scratch

%runscript
    # Default command to run when the container is executed
    exec /bin/bash
```
- Build the sif file based on the def file: `apptainer build ollama.sif ollama.def`
The tmp issue is gone. Now I am trying to figure out how to use the Apptainer with the GPU nodes.
EDIT: `ollama run llama3` works !


## Post 3 by @Berk.Gercek (2024-07-26T11:56:53.406Z)

Good info! For running apptainer with GPU access I’ve found (albeit a few months ago) that using the `--nv` flag on `apptainer run`, `apptainer exec`, etc. works fine. Just make sure to use the `module` command to load CUDA into the job env first, as the apptainer integrations rely on the native CUDA for communicating with GPUs.


## Post 4 by @Jamil.Zaghir (2024-07-26T12:28:23.532Z)

That is right ! Ollama does not have access to GPUs without this flag. Thank you for completing this thread with this crucial info.


## Post 5 by @Adrien.Albert (2024-07-26T12:31:13.042Z)

Hello
I am glade to see it’s working, thank you for your research and, @Berk.Gercek[@Berk.Gercek](https://hpc-community.unige.ch/u/berk.gercek), your help.
Jamil.Zaghir:
```
%runscript
    # Default command to run when the container is executed
    exec /bin/bash
```
In the runscript section, did you try to exec ollama dierectly ?


## Post 6 by @Jamil.Zaghir (2024-07-26T12:35:54.578Z)

I did not. I worked using apptainer shell to run ollama (and debug).


## Post 7 by @ahmad.alhineidi (2024-08-06T14:28:34.538Z)

I was able to run it by downloading the binary of ollama manually from here: ollama/docs/linux.md at main · ollama/ollama · GitHub[ollama/docs/linux.md at main · ollama/ollama · GitHub](https://github.com/ollama/ollama/blob/main/docs/linux.md)
so no need for sudo privilege on the HPC to install it. Then in an interactive slurm job or a normal sbatch one, you can ollama serve &, then run any ollama models for whatever. I hope this helps


## Post 8 by @Yann.Sagon (2024-08-15T13:53:45.641Z)

my two cents:
Instead of modifying the sif image, you could as well use the option mount and bind option of apptainer.
https://apptainer.org/docs/user/main/bind_paths_and_mounts.html#user-defined-bind-paths[https://apptainer.org/docs/user/main/bind_paths_and_mounts.html#user-defined-bind-paths](https://apptainer.org/docs/user/main/bind_paths_and_mounts.html#user-defined-bind-paths)
```
apptainer exec --bind /tmp:/scratch my_container.sif
```
or unset `TMPDIR`
```
TMPDIR='' apptainer .. 
```
