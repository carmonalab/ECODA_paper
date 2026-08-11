# Multiple GPU nodes doesn't detect CUDA

- Source: https://hpc-community.unige.ch/t/3713

- Created: 2024-11-07T16:43:27.813Z

- Tags: baobab

- Posts: 5

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Nizar.Michaud (2024-11-07T16:43:27.864Z)

## Primary Information
- Username: michaun3
- Cluster: Baobab
## Description
I’m unable to access GPUs on the shared-gpu partition. Across several nodes, `nvidia-smi` does not detect any GPUs, and CUDA is also not accessible within PyTorch.
## Steps to Reproduce
- Request a GPU node on the shared-gpu partition:
```
salloc --partition=shared-gpu --gres=gpu --time=00:30:00 --nodelist=gpu046
```
- Once allocated, check GPU visibility:
```
nvidia-smi
```
- Run a container with NVIDIA support:
```
apptainer run --nv sgwh_hpc.simg
```
- Inside the container, check `nvidia-smi` again and test CUDA access in PyTorch:
```
nvidia-smi
python -c "import torch; print(torch.cuda.is_available())"
```
## Expected Result
- `nvidia-smi` should list available GPUs.
- `torch.cuda.is_available()` should return `True`.
## Actual Result
- `nvidia-smi` outputs “No devices were found” across multiple nodes, including gpu046, gpu022, gpu032, gpu029, gpu033, gpu045, gpu017, gpu021, gpu026, gpu035.
- Inside the container, `torch.cuda.is_available()` returns `False`, indicating no CUDA access.
Could you please verify the NVIDIA driver installation and configuration on the shared-gpu partition? It seems the GPUs may not be correctly configured or exposed across nodes or in the container.
Thank you.


## Post 2 by @Adrien.Albert (2024-11-07T18:58:07.380Z)

Hi @Nizar.Michaud[@Nizar.Michaud](https://hpc-community.unige.ch/u/nizar.michaud)
On baobab, cpu046 i can list gpu with nvidia-smi:
```
(baobab)-[alberta@login1 ~]$ salloc --partition=shared-gpu --gres=gpu --time=00:30:00 --nodelist=gpu046
salloc: Pending job allocation 13543949
salloc: job 13543949 queued and waiting for resources
salloc: job 13543949 has been allocated resources
salloc: Granted job allocation 13543949
salloc: Waiting for resource configuration
salloc: Nodes gpu046 are ready for job
(baobab)-[alberta@gpu046 ~]$ nvidia-smi
Thu Nov  7 19:51:25 2024       
+-----------------------------------------------------------------------------------------+
| NVIDIA-SMI 565.57.01              Driver Version: 565.57.01      CUDA Version: 12.7     |
|-----------------------------------------+------------------------+----------------------+
| GPU  Name                 Persistence-M | Bus-Id          Disp.A | Volatile Uncorr. ECC |
| Fan  Temp   Perf          Pwr:Usage/Cap |           Memory-Usage | GPU-Util  Compute M. |
|                                         |                        |               MIG M. |
|=========================================+========================+======================|
|   0  NVIDIA RTX A5500               On  |   00000000:81:00.0 Off |                  Off |
| 30%   31C    P8              8W /  230W |       7MiB /  24564MiB |      0%      Default |
|                                         |                        |                  N/A |
+-----------------------------------------+------------------------+----------------------+
                                                                                         
+-----------------------------------------------------------------------------------------+
| Processes:                                                                              |
|  GPU   GI   CI        PID   Type   Process name                              GPU Memory |
|        ID   ID                                                               Usage      |
|=========================================================================================|
|  No running processes found                                                             |
+-----------------------------------------------------------------------------------------+
```
Can you confirm that you don’t any gpu with this procedure ?
Next, please  could you give me a procedure to the reproduce with the container you are using


## Post 3 by @Nizar.Michaud (2024-11-08T08:42:34.106Z)

Hello @Adrien.Albert[@Adrien.Albert](https://hpc-community.unige.ch/u/Adrien.Albert)
Indeed, here’s my procedure for gpu046:
```
(baobab)-[michaun3@login1 ~]$ salloc --partition=shared-gpu --gres=VramPerGpu:10G --nodelist=gpu046 --time=00:30:00
salloc: Pending job allocation 13548579
salloc: job 13548579 queued and waiting for resources
salloc: job 13548579 has been allocated resources
salloc: Granted job allocation 13548579
salloc: Waiting for resource configuration
salloc: Nodes gpu046 are ready for job
(baobab)-[michaun3@gpu046 ~]$ nvidia-smi
No devices were found
(baobab)-[michaun3@gpu046 ~]$ apptainer run --nv sgwh_hpc.simg

==========
== CUDA ==
==========

CUDA Version 12.6.2

Container image Copyright (c) 2016-2023, NVIDIA CORPORATION & AFFILIATES. All rights reserved.

This container image and its contents are governed by the NVIDIA Deep Learning Container License.
By pulling and using the container, you accept the terms and conditions of this license:
https://developer.nvidia.com/ngc/nvidia-deep-learning-container-license

A copy of this license is made available in this container at /NGC-DL-CONTAINER-LICENSE for your convenience.

Apptainer> nvidia-smi
No devices were found
Apptainer> python
Python 3.10.12 (main, Sep 11 2024, 15:47:36) [GCC 11.4.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import torch
>>> torch.cuda.is_available()
False
```
And here is how to compile the singularity image:
```
salloc --partition=shared-cpu --time=01:00:00
export APPTAINER_TMPDIR=/scratch/michaun3/tmp_apptainer
mkdir -p "$APPTAINER_TMPDIR"
apptainer build --tmpdir $APPTAINER_TMPDIR sgwh_hpc.simg docker://endomorphin/sgwh_hpc:v1
apptainer run --nv sgwh_hpc.simg
```
Hope that helps


## Post 4 by @Adrien.Albert (2024-11-08T09:20:00.973Z)

Your salloc is not correct, you forget to ask 1 gpu
```
-(baobab)-[alberta@login1 ~]$ salloc --partition=shared-gpu --gres=VramPerGpu:10G --nodelist=gpu046 --time=00:30:00
+(baobab)-[alberta@login1 ~]$ salloc --partition=shared-gpu --gres=VramPerGpu:10GB --gpus=1 --nodelist=gpu046 --time=00:30:00
```
```
(baobab)-[alberta@login1 ~]$ salloc --partition=shared-gpu --gres=VramPerGpu:10GB --gpus=1 --nodelist=gpu046 --time=00:30:00
salloc: Pending job allocation 13548691
salloc: job 13548691 queued and waiting for resources
salloc: job 13548691 has been allocated resources
salloc: Granted job allocation 13548691
salloc: Waiting for resource configuration
salloc: Nodes gpu046 are ready for job
(baobab)-[alberta@gpu046 ~]$ nvidia-smi
Fri Nov  8 10:17:29 2024       
+-----------------------------------------------------------------------------------------+
| NVIDIA-SMI 565.57.01              Driver Version: 565.57.01      CUDA Version: 12.7     |
|-----------------------------------------+------------------------+----------------------+
| GPU  Name                 Persistence-M | Bus-Id          Disp.A | Volatile Uncorr. ECC |
| Fan  Temp   Perf          Pwr:Usage/Cap |           Memory-Usage | GPU-Util  Compute M. |
|                                         |                        |               MIG M. |
|=========================================+========================+======================|
|   0  NVIDIA RTX A5500               On  |   00000000:41:00.0 Off |                  Off |
| 30%   31C    P8              8W /  230W |       7MiB /  24564MiB |      0%      Default |
|                                         |                        |                  N/A |
+-----------------------------------------+------------------------+----------------------+
                                                                                         
+-----------------------------------------------------------------------------------------+
| Processes:                                                                              |
|  GPU   GI   CI        PID   Type   Process name                              GPU Memory |
|        ID   ID                                                               Usage      |
|=========================================================================================|
|  No running processes found                                                             |
+-----------------------------------------------------------------------------------------+
```


## Post 5 by @Nizar.Michaud (2024-11-08T09:35:58.952Z)

:man_facepalming:
Thank you so much!!
