# Bamboo GPU008: RuntimeError: CUDA driver initialization failed, you might not have a CUDA gpu

- Source: https://hpc-community.unige.ch/t/4121

- Created: 2025-10-20T09:34:56.640Z

- Tags: bamboo

- Posts: 9

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Raphael.Rubino (2025-10-20T09:34:56.712Z)

## Primary informations
Username: rubinor
Cluster: bamboo
## Description
At least one GPU is not working on bamboo:gpu008.
## Steps to Reproduce
When trying to run anything on GPU, for instance in Python:
```
import torch
print(torch.cuda.current_device())
```
```
RuntimeError: CUDA driver initialization failed, you might not have a CUDA gpu.
```
## Expected Result
No error.
## Actual Result
```
RuntimeError: CUDA driver initialization failed, you might not have a CUDA gpu.
```


## Post 2 by @Adrien.Albert (2025-10-20T16:30:32.387Z)

Hello @Raphael.Rubino[@Raphael.Rubino](https://hpc-community.unige.ch/u/raphael.rubino)
## Tests
For me everything is working, with nvidia-smi you can seel all gpu:
```
(bamboo)-[alberta@gpu008 ~]$ nvidia-smi 
Mon Oct 20 17:20:25 2025       
+-----------------------------------------------------------------------------------------+
| NVIDIA-SMI 580.95.05              Driver Version: 580.95.05      CUDA Version: 13.0     |
+-----------------------------------------+------------------------+----------------------+
| GPU  Name                 Persistence-M | Bus-Id          Disp.A | Volatile Uncorr. ECC |
| Fan  Temp   Perf          Pwr:Usage/Cap |           Memory-Usage | GPU-Util  Compute M. |
|                                         |                        |               MIG M. |
|=========================================+========================+======================|
|   0  NVIDIA RTX PRO 6000 Blac...    On  |   00000000:04:00.0 Off |                    0 |
| N/A   33C    P8             34W /  600W |       0MiB /  97887MiB |      0%      Default |
|                                         |                        |             Disabled |
+-----------------------------------------+------------------------+----------------------+

+-----------------------------------------------------------------------------------------+
| Processes:                                                                              |
|  GPU   GI   CI              PID   Type   Process name                        GPU Memory |
|        ID   ID                                                               Usage      |
|=========================================================================================|
|  No running processes found                                                             |
+-----------------------------------------------------------------------------------------+
```
However it’s the new GPU model  `NVIDIA RTX PRO 6000 Blacwell` and I suspect the CUDA used with pytorch is not compatible:
PyTorch/2.1.2-CUDA-12.1.1
Only `CUDA/12.8.0` or newer seems working with these GPU model
I tested with a first  containr of `pytorch-cuda-12.8` but unfortunetly there is no support cuda capabilities for `sm_120`:
```
(bamboo)-[alberta@gpu008 ~]$ singularity exec --nv docker://pytorch/pytorch:2.7.0-cuda12.8-cudnn9-runtime bash
INFO:    Converting OCI blobs to SIF format
INFO:    Starting build...
INFO:    Fetching OCI image...
4.0GiB / 4.0GiB [================================================================================================================================================================================] 100 % 0.0 b/s 0s
6.9MiB / 6.9MiB [================================================================================================================================================================================] 100 % 0.0 b/s 0s
29.0MiB / 29.0MiB [==============================================================================================================================================================================] 100 % 0.0 b/s 0s
INFO:    Extracting OCI image...
INFO:    Inserting Apptainer configuration...
INFO:    Creating SIF file...
[========================================================================================================================================================================================================] 100 % 0s
Apptainer> python --version
Python 3.11.12
Apptainer> python
Python 3.11.12 | packaged by conda-forge | (main, Apr 10 2025, 22:23:25) [GCC 13.3.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import torch
>>> print(torch.cuda.current_device())
/home/users/a/alberta/.local/lib/python3.11/site-packages/torch/cuda/__init__.py:173: UserWarning: 
NVIDIA RTX PRO 6000 Blackwell Server Edition with CUDA capability sm_120 is not compatible with the current PyTorch installation.
The current PyTorch install supports CUDA capabilities sm_37 sm_50 sm_60 sm_70 sm_75 sm_80 sm_86.
If you want to use the NVIDIA RTX PRO 6000 Blackwell Server Edition GPU with PyTorch, please check the instructions at https://pytorch.org/get-started/locally/

  warnings.warn(incompatible_device_warn.format(device_name, capability, " ".join(arch_list), device_name))
0
```
So I tried with `nvidia-ml-py3` and it seems working, (that means my GPU is here):
```
(bamboo)-[alberta@gpu008 ~]$ ml

Currently Loaded Modules:
  1) GCCcore/14.3.0   3) binutils/2.44   5) CUDA/12.8.0   7) ncurses/6.5       9) libtommath/1.3.0  11) SQLite/3.50.1  13) libffi/3.5.1  15) Python/3.13.5
  2) zlib/1.3.1       4) GCC/14.3.0      6) bzip2/1.0.8   8) libreadline/8.2  10) Tcl/9.0.1         12) XZ/5.8.1       14) OpenSSL/3

(bamboo)-[alberta@gpu008 ~]$ pip install nvidia-ml-py3
Defaulting to user installation because normal site-packages is not writeable
Collecting nvidia-ml-py3
  Downloading nvidia-ml-py3-7.352.0.tar.gz (19 kB)
  Preparing metadata (setup.py) ... done
Building wheels for collected packages: nvidia-ml-py3
  DEPRECATION: Building 'nvidia-ml-py3' using the legacy setup.py bdist_wheel mechanism, which will be removed in a future version. pip 25.3 will enforce this behaviour change. A possible replacement is to use the standardized build interface by setting the `--use-pep517` option, (possibly combined with `--no-build-isolation`), or adding a `pyproject.toml` file to the source tree of 'nvidia-ml-py3'. Discussion can be found at https://github.com/pypa/pip/issues/6334
  Building wheel for nvidia-ml-py3 (setup.py) ... done
  Created wheel for nvidia-ml-py3: filename=nvidia_ml_py3-7.352.0-py3-none-any.whl size=19208 sha256=d59c1625dc66ce844f0db0784f6524df62914b0a6922352d992c5f7443bfb527
  Stored in directory: /home/users/a/alberta/.cache/pip/wheels/ea/47/38/29179ca914d95f79296647a42943b8e576dc9d318f94bad57a
Successfully built nvidia-ml-py3
Installing collected packages: nvidia-ml-py3
Successfully installed nvidia-ml-py3-7.352.0

[notice] A new release of pip is available: 25.1.1 -> 25.2
[notice] To update, run: pip install --upgrade pip
(bamboo)-[alberta@gpu008 ~]$ python
Python 3.13.5 (main, Oct 14 2025, 11:34:05) [GCC 14.3.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import pynvml
... pynvml.nvmlInit()
... count = pynvml.nvmlDeviceGetCount()
... print(f"GPU detected : {count}")
... for i in range(count):
...     handle = pynvml.nvmlDeviceGetHandleByIndex(i)
...     name = pynvml.nvmlDeviceGetName(handle)
...     print(f"GPU {i}: {name.decode()}")
...     
GPU detected : 1
GPU 0: NVIDIA RTX PRO 6000 Blackwell Server Edition
```
For now it’s not possible to run pytorch-cuda on these gpu nodes until we re-compile pytorch available with `module`.
## Work Arround | solution
But I found a container provided by `Nvidia` and it seems working:
```
(bamboo)-[alberta@gpu008 ~]$ singularity exec --nv docker://nvcr.io/nvidia/pytorch:25.09-py3 bash
INFO:    Using cached SIF image
Apptainer> python
Python 3.12.3 (main, Aug 14 2025, 17:47:21) [GCC 13.3.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import torch
/usr/local/lib/python3.12/dist-packages/torch/cuda/__init__.py:63: FutureWarning: The pynvml package is deprecated. Please install nvidia-ml-py instead. If you did not install pynvml directly, please report this to the maintainers of the package that installed pynvml for you.
  import pynvml  # type: ignore[import]

>>> print(torch.cuda.current_device())
0
>>> torch.cuda.device_count()
1
>>> torch.cuda.get_device_name(0)
'NVIDIA RTX PRO 6000 Blackwell Server Edition'
```
Let me know how it’s working for you.


## Post 3 by @Raphael.Rubino (2025-10-20T16:37:25.445Z)

Thank you for checking gpu008.
Sorry, I should have said that I could run jobs on 2 GPUs of gpu008, but my jobs failed on at least one GPU of gpu008. In other words, not all GPUs on gpu008 are faulty but there is at least one faulty GPU on gpu008.


## Post 4 by @Adrien.Albert (2025-10-22T09:31:20.570Z)

Hi @Raphael.Rubino[@Raphael.Rubino](https://hpc-community.unige.ch/u/raphael.rubino)
The server is draining; I will run a test.
Thank you for reporting


## Post 5 by @Ramon.CalvoGonzalez (2025-11-05T10:51:00.661Z)

Hi,
I am seeing the same problem in gpu008. The GPUs are detected (reported by `nvidia-smi` and also `torch.cuda.device_count()`). But when I try to actually create a CUDA context, it fails:
```
Python 3.12.11 (main, Jul 11 2025, 22:43:48) [Clang 20.1.4 ] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import torch
>>> torch.cuda.is_available()
/home/users/c/calvogon/ssl-jax/.venv/lib/python3.12/site-packages/torch/cuda/__init__.py:182: UserWarning: CUDA initialization: CUDA driver initialization failed, you might not have a CUDA gpu. (Triggered internally at /pytorch/c10/cuda/CUDAFunctions.cpp:119.)
  return torch._C._cuda_getDeviceCount() > 0
False
```
The same job I am trying to run works on gpu007, which afaict is an equivalent node.


## Post 6 by @Ramon.CalvoGonzalez (2025-11-05T11:06:19.982Z)

I have run `torch.cuda.is_available()'` with every possible combination of `CUDA_VISIBLE_DEVICES`, and I get this CUDA error only when there are 3 or more GPUs in `CUDA_VISIBLE_DEVICES`.
The topology (`nvidia-smi topo -m`) is different to the one of gpu007.
gpu007:
```
        GPU0    GPU1    GPU2    GPU3    NIC0    CPU Affinity    NUMA Affinity   GPU NUMA ID
GPU0     X      SYS     SYS     SYS     PHB     96-107  3               N/A
GPU1    SYS      X      SYS     SYS     SYS     64-75   2               N/A
GPU2    SYS     SYS      X      SYS     SYS     0-11    0               N/A
GPU3    SYS     SYS     SYS      X      SYS     32-43   1               N/A
NIC0    PHB     SYS     SYS     SYS      X
```
gpu008:
```
        GPU0    GPU1    GPU2    GPU3    NIC0    CPU Affinity    NUMA Affinity   GPU NUMA ID
GPU0     X      PIX     SYS     SYS     SYS     48-63   3               N/A
GPU1    PIX      X      SYS     SYS     SYS     48-63   3               N/A
GPU2    SYS     SYS      X      PIX     SYS     64-79   4               N/A
GPU3    SYS     SYS     PIX      X      SYS     64-79   4               N/A
NIC0    SYS     SYS     SYS     SYS      X
```


## Post 7 by @Yann.Sagon (2025-11-06T10:05:29.837Z)

Dear @Ramon.CalvoGonzalez[@Ramon.CalvoGonzalez](https://hpc-community.unige.ch/u/ramon.calvogonzalez) and @Raphael.Rubino[@Raphael.Rubino](https://hpc-community.unige.ch/u/raphael.rubino)
I was able to reproduce the issue even without slurm and pytorch. It isn’t possible to use more than two gpus at the same time. This is not the case on gpu007, probably because this node has only one physical cpu.
I’ve implemented a workaround on gpu008 which consists to disable hmm as this is a known issue: CUDA initialization failiure on latest release drivers · Issue #797 · NVIDIA/open-gpu-kernel-modules · GitHub[CUDA initialization failiure on latest release drivers · Issue #797 · NVIDIA/open-gpu-kernel-modules · GitHub](https://github.com/NVIDIA/open-gpu-kernel-modules/issues/797)


## Post 8 by @Ramon.CalvoGonzalez (2025-11-06T21:21:43.440Z)

Thanks!
Just FYI, both nodes are currently inaccessible due to these reasons:
- gpu007 (reserved): reserved for `sagon_1288`
- gpu008 (down): Node unexpectedly rebooted [slurm@2025-11-06T09:10:34]
Best,
Ramon.


## Post 9 by @Yann.Sagon (2025-11-07T07:54:49.946Z)

Ramon.CalvoGonzalez:
> both nodes are currently inaccessible
@Ramon.CalvoGonzalez[@Ramon.CalvoGonzalez](https://hpc-community.unige.ch/u/ramon.calvogonzalez) oups, fixed.
