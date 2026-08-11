# [Tutorial] Create a Custom Kernel for JupyterLab Using cotainr

- Source: https://hpc-community.unige.ch/t/4251

- Created: 2026-03-16T11:05:36.711Z

- Posts: 9

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yann.Sagon (2026-03-16T11:05:36.867Z)

This tutorial explains how to build a GPU‑ready container image using cotainr[cotainr](https://doc.eresearch.unige.ch/hpc/applications_and_libraries#conda), install PyTorch and other Python packages inside it, and expose it as a custom Jupyter kernel in an HPC environment using Open OnDemand[Open OnDemand](https://doc.eresearch.unige.ch/hpc/how_to_use_openondemand).
---
## 1. Create the Conda Environment File (`env.yml`)
Create a file named `env.yml` with the following content:
```
name: torch-gpu
channels:
  - pytorch
  - nvidia
  - conda-forge
  - defaults

dependencies:
  - python=3.10
  - pip
  - pytorch
  - torchvision
  - torchaudio
  - pytorch-cuda=12.1
  - ipykernel
 ### optional packages! ###
  - pip:
      - timm
      - torchinfo
      - thop
      - "flwr[simulation]"
      - flwr-datasets
```
Save it in your home directory, e.g. `~/tutorial/env.yml`.
---
## 2. Build the Apptainer `.sif` Image Using cotainr
Use a CUDA runtime image so cotainr can install Conda cleanly:
```
cotainr build pytorchkernel.sif   --base-image=docker://nvidia/cuda:12.1.1-runtime-ubuntu22.04   --accept-licenses   --conda-env=~/tutorial/env.yml
```
cotainr will install Miniforge, create the Conda environment `torch-gpu`, and install all required packages.
---
## 3. Create the Kernel Specification
Create the folder:
```
mkdir -p ~/.local/share/jupyter/kernels/torch-gpu-apptainer
```
Create `kernel.json` inside it:
```
{
  "argv": [
    "/home/YOURUSER/.local/share/jupyter/kernels/torch-gpu-apptainer/init_kernel.sh",
    "-f",
    "{connection_file}"
  ],
  "display_name": "PyTorch GPU (Apptainer)",
  "language": "python"
}
```
Replace `YOURUSER` with your username.
---
## 4. Create the `init_kernel.sh` Script
Create the file:
```
~/.local/share/jupyter/kernels/torch-gpu-apptainer/init_kernel.sh
```
Add:
```
#!/bin/bash

SIF="/home/YOURUSER/tutorial/pytorchkernel.sif"

exec apptainer exec --cleanenv --nv "$SIF"     conda run -n torch-gpu python -m ipykernel "$@"
```
Make it executable:
```
chmod +x ~/.local/share/jupyter/kernels/torch-gpu-apptainer/init_kernel.sh
```
---
## 5. Launch JupyterLab via Open OnDemand
- Open your Open OnDemand portal.
- Go to Interactive Apps → JupyterLab.
- Start a session with GPU resources.
---
## 6. Create a New Notebook and Select the Custom Kernel
In JupyterLab:
- Open the Launcher.
- Choose “PyTorch GPU (Apptainer)” (see picture below)
- Run a test:
```
import torch
print(torch.cuda.is_available())
print(torch.cuda.get_device_name(0))
```
image
image603×297 11.9 KB
[image603×297 11.9 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/b5a815a157f1205a39d56c848e7b1b0778f45c5a.png)
If everything works, your containerized kernel is ready.
---
## :tada: Done!
You now have a reproducible, GPU-enabled Apptainer image and a custom Jupyter kernel fully integrated with Open OnDemand.
ref:
- Enable CUDA in jupyterlab environment? - #2 by Yann.Sagon[Enable CUDA in jupyterlab environment? - #2 by Yann.Sagon](https://hpc-community.unige.ch/t/enable-cuda-in-jupyterlab-environment/3733/2)
- [tutorial] load modules from JuypterLab launched on OpenOnDemand[[tutorial] load modules from JuypterLab launched on OpenOnDemand](https://hpc-community.unige.ch/t/tutorial-load-modules-from-juypterlab-launched-on-openondemand/3241)
- Baobab OnDemand: jupyterLab environment-specific kernel never connects[Baobab OnDemand: jupyterLab environment-specific kernel never connects](https://hpc-community.unige.ch/t/baobab-ondemand-jupyterlab-environment-specific-kernel-never-connects/3348)
- Installing specific versions of pytorch and transformers via pip - #5 by Adrien.Albert[Installing specific versions of pytorch and transformers via pip - #5 by Adrien.Albert](https://hpc-community.unige.ch/t/installing-specific-versions-of-pytorch-and-transformers-via-pip/4065/5)


## Post 2 by @Sajad.Hamzenejadi (2026-04-16T07:52:08.274Z)

Hello @Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon)
Many thanks for your instructions. I followed the exact instructions (baobab), and did the step 6 to check if my containerized kernel is ready
> import torch
> print(torch.cuda.is_available())
> print(torch.cuda.get_device_name(0))
But I got the following error. I do not know which step is broken:
```
---------------------------------------------------------------------------
ModuleNotFoundError                       Traceback (most recent call last)
Cell In[1], line 1
----> 1 import torch
      2 print(torch.cuda.is_available())
      3 print(torch.cuda.get_device_name(0))

ModuleNotFoundError: No module named 'torch'
```


## Post 3 by @Sajad.Hamzenejadi (2026-04-20T08:30:07.641Z)

Hello!
If anyone from the HPC team has suggestions on how to resolve the mentioned issue, I would greatly appreciate your input!
Thank you in advance for your help and support!


## Post 4 by @Yann.Sagon (2026-04-21T08:16:44.116Z)

Sajad.Hamzenejadi:
> I followed the exact instructions (baobab), and did the step 6 to check if my containerized kernel is ready
Did you do the previous steps too? I’ve checked on Baobab in your home directory, the scripts and directory from the previous steps aren’t there.


## Post 6 by @Sajad.Hamzenejadi (2026-04-21T09:09:44.586Z)

Hello Yann,
Many thanks for your response. Yes, I did all steps (1-6).
What do you suggest at this point? Should I re-do steps 1-6?


## Post 7 by @Yann.Sagon (2026-04-21T11:18:04.135Z)

Sajad.Hamzenejadi:
> Yes, I did all steps (1-6).
Then where is this file?
```
~/.local/share/jupyter/kernels/torch-gpu-apptainer/init_kernel.sh
```


## Post 8 by @Sajad.Hamzenejadi (2026-04-21T11:34:46.406Z)

You’re right, I cannot find it as well :frowning:


## Post 9 by @Sajad.Hamzenejadi (2026-04-23T10:30:30.496Z)

Hello @Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon),
I have repeated all of the above steps, now I have the kernel “PyTorch GPU (Apptainer)“ available, but the notebook cannot connect to kernel and the kernel status changes from to “connecting“ to “disconnected”.
I would greatly appreciate your time and support to help me solving this issue.


## Post 10 by @Yann.Sagon (2026-04-23T12:35:33.519Z)

Hello,
can you please modify your file `init_kernel.sh` with the following content:
```
#!/bin/bash

# Emplacement de ton conteneur cotainr
SIF="/home/users/h/hamzenej/pytorchkernel.sif"

# Lancer le kernel Python dans le conteneur
exec apptainer exec --cleanenv --nv "$SIF" \
    python -m ipykernel "$@"
```
