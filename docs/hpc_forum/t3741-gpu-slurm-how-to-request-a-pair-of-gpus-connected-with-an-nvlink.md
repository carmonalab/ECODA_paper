# [GPU][SLURM] How to request a pair of GPUs connected with an NVLINK?

- Source: https://hpc-community.unige.ch/t/3741

- Created: 2024-11-21T13:27:24.678Z

- Posts: 6

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @maciej.falkiewicz (2024-11-21T13:27:24.742Z)

Dear @support[@support](https://hpc-community.unige.ch/groups/support) ,
I want to use the gpu048 node of Baobab announced here[here](https://hpc-community.unige.ch//hpc-community.unige.ch/t/new-compute-installed-gpu048-baobab/3640), and would like to benefit from the NVLINK when running my job on two GPUs.
How should I schedule my job in order to get GPUs that are connected?
Thank you in advance for your support,
Kind regards,
Maciej Falkiewicz


## Post 2 by @Yann.Sagon (2024-11-28T13:43:35.892Z)

Dear @maciej.falkiewicz[@maciej.falkiewicz](https://hpc-community.unige.ch/u/maciej.falkiewicz)
You can check the NVlink status like that:
```
 nvidia-smi nvlink -s

(baobab)-[root@gpu048 ~]$ nvidia-smi nvlink -s
GPU 0: NVIDIA RTX A6000 (UUID: GPU-1b574915-9e55-3096-1fdb-d32d0754017c)
         Link 0: 14.062 GB/s
         Link 1: 14.062 GB/s
         Link 2: 14.062 GB/s
         Link 3: 14.062 GB/s
GPU 1: NVIDIA RTX A6000 (UUID: GPU-f0d29aba-53cd-2500-951e-21767ef0871d)
         Link 0: 14.062 GB/s
         Link 1: 14.062 GB/s
         Link 2: 14.062 GB/s
         Link 3: 14.062 GB/s
GPU 2: NVIDIA RTX A6000 (UUID: GPU-73e74607-9f67-c6a1-a713-81072f9aa686)
         Link 0: 14.062 GB/s
         Link 1: 14.062 GB/s
         Link 2: 14.062 GB/s
         Link 3: 14.062 GB/s
GPU 3: NVIDIA RTX A6000 (UUID: GPU-d10ff217-cca3-6530-0d49-4257fd2f221b)
         Link 0: 14.062 GB/s
         Link 1: 14.062 GB/s
         Link 2: 14.062 GB/s
         Link 3: 14.062 GB/s
GPU 4: NVIDIA RTX A6000 (UUID: GPU-0606b155-03ff-5431-a941-c5edff4461a2)
         Link 0: 14.062 GB/s
         Link 1: 14.062 GB/s
         Link 2: 14.062 GB/s
         Link 3: 14.062 GB/s
GPU 5: NVIDIA RTX A6000 (UUID: GPU-cf13b5b7-9865-4e81-e2c5-818a37ce63e8)
         Link 0: 14.062 GB/s
         Link 1: 14.062 GB/s
         Link 2: 14.062 GB/s
         Link 3: 14.062 GB/s
GPU 6: NVIDIA RTX A6000 (UUID: GPU-ea7d0f44-8883-589a-ca19-f92cdebccabd)
         Link 0: 14.062 GB/s
         Link 1: 14.062 GB/s
         Link 2: 14.062 GB/s
         Link 3: 14.062 GB/s
GPU 7: NVIDIA RTX A6000 (UUID: GPU-971d3b23-dc8f-0f96-5271-3eb4060cff19)
         Link 0: 14.062 GB/s
         Link 1: 14.062 GB/s
         Link 2: 14.062 GB/s
         Link 3: 14.062 GB/s
```
We’ve now installed `CUDA-Samples/12.1-CUDA-12.1.1` that provides a tool named `deviceQuery`.
This is the output from it:
```
deviceQuery Starting...

 CUDA Device Query (Runtime API) version (CUDART static linking)

Detected 8 CUDA Capable device(s)

Device 0: "NVIDIA RTX A6000"
  CUDA Driver Version / Runtime Version          12.7 / 12.1
  CUDA Capability Major/Minor version number:    8.6
  Total amount of global memory:                 48790 MBytes (51160023040 bytes)
  (084) Multiprocessors, (128) CUDA Cores/MP:    10752 CUDA Cores
  GPU Max Clock rate:                            1800 MHz (1.80 GHz)
  Memory Clock rate:                             8001 Mhz
  Memory Bus Width:                              384-bit
  L2 Cache Size:                                 6291456 bytes
  Maximum Texture Dimension Size (x,y,z)         1D=(131072), 2D=(131072, 65536), 3D=(16384, 16384, 16384)
  Maximum Layered 1D Texture Size, (num) layers  1D=(32768), 2048 layers
  Maximum Layered 2D Texture Size, (num) layers  2D=(32768, 32768), 2048 layers
  Total amount of constant memory:               65536 bytes
  Total amount of shared memory per block:       49152 bytes
  Total shared memory per multiprocessor:        102400 bytes
  Total number of registers available per block: 65536
  Warp size:                                     32
  Maximum number of threads per multiprocessor:  1536
  Maximum number of threads per block:           1024
  Max dimension size of a thread block (x,y,z): (1024, 1024, 64)
  Max dimension size of a grid size    (x,y,z): (2147483647, 65535, 65535)
  Maximum memory pitch:                          2147483647 bytes
  Texture alignment:                             512 bytes
  Concurrent copy and kernel execution:          Yes with 3 copy engine(s)
  Run time limit on kernels:                     No
  Integrated GPU sharing Host Memory:            No
  Support host page-locked memory mapping:       Yes
  Alignment requirement for Surfaces:            Yes
  Device has ECC support:                        Disabled
  Device supports Unified Addressing (UVA):      Yes
  Device supports Managed Memory:                Yes
  Device supports Compute Preemption:            Yes
  Supports Cooperative Kernel Launch:            Yes
  Supports MultiDevice Co-op Kernel Launch:      Yes
  Device PCI Domain ID / Bus ID / location ID:   0 / 1 / 0
  Compute Mode:
     < Default (multiple host threads can use ::cudaSetDevice() with device simultaneously) >

[...]

> Peer access from NVIDIA RTX A6000 (GPU0) -> NVIDIA RTX A6000 (GPU1) : Yes
> Peer access from NVIDIA RTX A6000 (GPU0) -> NVIDIA RTX A6000 (GPU2) : Yes
> Peer access from NVIDIA RTX A6000 (GPU0) -> NVIDIA RTX A6000 (GPU3) : Yes
> Peer access from NVIDIA RTX A6000 (GPU0) -> NVIDIA RTX A6000 (GPU4) : Yes
> Peer access from NVIDIA RTX A6000 (GPU0) -> NVIDIA RTX A6000 (GPU5) : Yes
> Peer access from NVIDIA RTX A6000 (GPU0) -> NVIDIA RTX A6000 (GPU6) : Yes
> Peer access from NVIDIA RTX A6000 (GPU0) -> NVIDIA RTX A6000 (GPU7) : Yes
> Peer access from NVIDIA RTX A6000 (GPU1) -> NVIDIA RTX A6000 (GPU0) : Yes
> Peer access from NVIDIA RTX A6000 (GPU1) -> NVIDIA RTX A6000 (GPU2) : Yes
> Peer access from NVIDIA RTX A6000 (GPU1) -> NVIDIA RTX A6000 (GPU3) : Yes
> Peer access from NVIDIA RTX A6000 (GPU1) -> NVIDIA RTX A6000 (GPU4) : Yes
> Peer access from NVIDIA RTX A6000 (GPU1) -> NVIDIA RTX A6000 (GPU5) : Yes
> Peer access from NVIDIA RTX A6000 (GPU1) -> NVIDIA RTX A6000 (GPU6) : Yes
> Peer access from NVIDIA RTX A6000 (GPU1) -> NVIDIA RTX A6000 (GPU7) : Yes
> Peer access from NVIDIA RTX A6000 (GPU2) -> NVIDIA RTX A6000 (GPU0) : Yes
> Peer access from NVIDIA RTX A6000 (GPU2) -> NVIDIA RTX A6000 (GPU1) : Yes
> Peer access from NVIDIA RTX A6000 (GPU2) -> NVIDIA RTX A6000 (GPU3) : Yes
> Peer access from NVIDIA RTX A6000 (GPU2) -> NVIDIA RTX A6000 (GPU4) : Yes
> Peer access from NVIDIA RTX A6000 (GPU2) -> NVIDIA RTX A6000 (GPU5) : Yes
> Peer access from NVIDIA RTX A6000 (GPU2) -> NVIDIA RTX A6000 (GPU6) : Yes
> Peer access from NVIDIA RTX A6000 (GPU2) -> NVIDIA RTX A6000 (GPU7) : Yes
> Peer access from NVIDIA RTX A6000 (GPU3) -> NVIDIA RTX A6000 (GPU0) : Yes
> Peer access from NVIDIA RTX A6000 (GPU3) -> NVIDIA RTX A6000 (GPU1) : Yes
> Peer access from NVIDIA RTX A6000 (GPU3) -> NVIDIA RTX A6000 (GPU2) : Yes
> Peer access from NVIDIA RTX A6000 (GPU3) -> NVIDIA RTX A6000 (GPU4) : Yes
> Peer access from NVIDIA RTX A6000 (GPU3) -> NVIDIA RTX A6000 (GPU5) : Yes
> Peer access from NVIDIA RTX A6000 (GPU3) -> NVIDIA RTX A6000 (GPU6) : Yes
> Peer access from NVIDIA RTX A6000 (GPU3) -> NVIDIA RTX A6000 (GPU7) : Yes
> Peer access from NVIDIA RTX A6000 (GPU4) -> NVIDIA RTX A6000 (GPU0) : Yes
> Peer access from NVIDIA RTX A6000 (GPU4) -> NVIDIA RTX A6000 (GPU1) : Yes
> Peer access from NVIDIA RTX A6000 (GPU4) -> NVIDIA RTX A6000 (GPU2) : Yes
> Peer access from NVIDIA RTX A6000 (GPU4) -> NVIDIA RTX A6000 (GPU3) : Yes
> Peer access from NVIDIA RTX A6000 (GPU4) -> NVIDIA RTX A6000 (GPU5) : Yes
> Peer access from NVIDIA RTX A6000 (GPU4) -> NVIDIA RTX A6000 (GPU6) : Yes
> Peer access from NVIDIA RTX A6000 (GPU4) -> NVIDIA RTX A6000 (GPU7) : Yes
> Peer access from NVIDIA RTX A6000 (GPU5) -> NVIDIA RTX A6000 (GPU0) : Yes
> Peer access from NVIDIA RTX A6000 (GPU5) -> NVIDIA RTX A6000 (GPU1) : Yes
> Peer access from NVIDIA RTX A6000 (GPU5) -> NVIDIA RTX A6000 (GPU2) : Yes
> Peer access from NVIDIA RTX A6000 (GPU5) -> NVIDIA RTX A6000 (GPU3) : Yes
> Peer access from NVIDIA RTX A6000 (GPU5) -> NVIDIA RTX A6000 (GPU4) : Yes
> Peer access from NVIDIA RTX A6000 (GPU5) -> NVIDIA RTX A6000 (GPU6) : Yes
> Peer access from NVIDIA RTX A6000 (GPU5) -> NVIDIA RTX A6000 (GPU7) : Yes
> Peer access from NVIDIA RTX A6000 (GPU6) -> NVIDIA RTX A6000 (GPU0) : Yes
> Peer access from NVIDIA RTX A6000 (GPU6) -> NVIDIA RTX A6000 (GPU1) : Yes
> Peer access from NVIDIA RTX A6000 (GPU6) -> NVIDIA RTX A6000 (GPU2) : Yes
> Peer access from NVIDIA RTX A6000 (GPU6) -> NVIDIA RTX A6000 (GPU3) : Yes
> Peer access from NVIDIA RTX A6000 (GPU6) -> NVIDIA RTX A6000 (GPU4) : Yes
> Peer access from NVIDIA RTX A6000 (GPU6) -> NVIDIA RTX A6000 (GPU5) : Yes
> Peer access from NVIDIA RTX A6000 (GPU6) -> NVIDIA RTX A6000 (GPU7) : Yes
> Peer access from NVIDIA RTX A6000 (GPU7) -> NVIDIA RTX A6000 (GPU0) : Yes
> Peer access from NVIDIA RTX A6000 (GPU7) -> NVIDIA RTX A6000 (GPU1) : Yes
> Peer access from NVIDIA RTX A6000 (GPU7) -> NVIDIA RTX A6000 (GPU2) : Yes
> Peer access from NVIDIA RTX A6000 (GPU7) -> NVIDIA RTX A6000 (GPU3) : Yes
> Peer access from NVIDIA RTX A6000 (GPU7) -> NVIDIA RTX A6000 (GPU4) : Yes
> Peer access from NVIDIA RTX A6000 (GPU7) -> NVIDIA RTX A6000 (GPU5) : Yes
> Peer access from NVIDIA RTX A6000 (GPU7) -> NVIDIA RTX A6000 (GPU6) : Yes

deviceQuery, CUDA Driver = CUDART, CUDA Driver Version = 12.7, CUDA Runtime Version = 12.1, NumDevs = 8
Result = PASS
```
As far as I understand, each GPU is connected to each GPU.
Best


## Post 3 by @maciej.falkiewicz (2024-11-29T12:04:27.924Z)

Thank you @Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon) for the reply.
First of all, YOU can check it by being the root user, I cannot :slight_smile:
Then, I am afraid that what you say is not true. NVLink connects only two GPUs, please check the documentation here: NVIDIA RTX A6000 Datasheet[NVIDIA RTX A6000 Datasheet](https://resources.nvidia.com/en-us-briefcase-for-datasheets/proviz-print-nvidia-1?ncid=no-ncid)
Could you please show the output of `nvidia-smi topo -m`? `nvidia-smi` documentation just in case: https://docs.nvidia.com/deploy/nvidia-smi/index.html[https://docs.nvidia.com/deploy/nvidia-smi/index.html](https://docs.nvidia.com/deploy/nvidia-smi/index.html)
Thank you in advance,
Kind regards,
Maciej Falkiewicz


## Post 4 by @Yann.Sagon (2024-11-29T14:38:04.315Z)

Hello @maciej.falkiewicz[@maciej.falkiewicz](https://hpc-community.unige.ch/u/maciej.falkiewicz)
maciej.falkiewicz:
> YOU can check it by being the root user, I cannot
You are right, this is why I pasted the output in the post:)
maciej.falkiewicz:
> Then, I am afraid that what you say is not true
I was trying to understand the output, but yes, this makes sense, only two GPUs can be connected through NVLink, but each GPU is connected through PCI or other, probably with less BW.
maciej.falkiewicz:
> Could you please show the output of `nvidia-smi topo -m`
Yes sure.
```
(baobab)-[root@gpu048 ~]$ nvidia-smi topo -m
        GPU0    GPU1    GPU2    GPU3    GPU4    GPU5    GPU6    GPU7    NIC0    NIC1    CPU Affinity    NUMA Affinity   GPU NUMA ID
GPU0     X      NV4     NODE    NODE    SYS     SYS     SYS     SYS     SYS     SYS     0-63    0               N/A
GPU1    NV4      X      NODE    NODE    SYS     SYS     SYS     SYS     SYS     SYS     0-63    0               N/A
GPU2    NODE    NODE     X      NV4     SYS     SYS     SYS     SYS     SYS     SYS     0-63    0               N/A
GPU3    NODE    NODE    NV4      X      SYS     SYS     SYS     SYS     SYS     SYS     0-63    0               N/A
GPU4    SYS     SYS     SYS     SYS      X      NV4     NODE    NODE    NODE    NODE    64-127  1               N/A
GPU5    SYS     SYS     SYS     SYS     NV4      X      NODE    NODE    NODE    NODE    64-127  1               N/A
GPU6    SYS     SYS     SYS     SYS     NODE    NODE     X      NV4     PHB     PHB     64-127  1               N/A
GPU7    SYS     SYS     SYS     SYS     NODE    NODE    NV4      X      NODE    NODE    64-127  1               N/A
NIC0    SYS     SYS     SYS     SYS     NODE    NODE    PHB     NODE     X      PIX
NIC1    SYS     SYS     SYS     SYS     NODE    NODE    PHB     NODE    PIX      X

Legend:

  X    = Self
  SYS  = Connection traversing PCIe as well as the SMP interconnect between NUMA nodes (e.g., QPI/UPI)
  NODE = Connection traversing PCIe as well as the interconnect between PCIe Host Bridges within a NUMA node
  PHB  = Connection traversing PCIe as well as a PCIe Host Bridge (typically the CPU)
  PXB  = Connection traversing multiple PCIe bridges (without traversing the PCIe Host Bridge)
  PIX  = Connection traversing at most a single PCIe bridge
  NV#  = Connection traversing a bonded set of # NVLinks

NIC Legend:

  NIC0: mlx5_0
  NIC1: mlx5_1
```
I’ve enabled a feature on gpu048: we are now discovering the GPUs instead of hard code their spec. The advantage is that this adds two extra flags: Cores and Links. See ref[ref](https://slurm.schedmd.com/gres.conf.html).
This is how slurm sees the GPUs now:
```
(baobab)-[root@gpu048 ~]$ slurmd -G
slurmd: gpu/nvml: _get_system_gpu_list_nvml: 8 GPU system device(s) detected
slurmd: Gres Name=gpu Type=nvidia_rtx_a6000 Count=1 Index=0 ID=7696487 File=/dev/nvidia0 Cores=0-63 CoreCnt=128 Links=-1,4,0,0,0,0,0,0 Flags=HAS_FILE,HAS_TYPE,ENV_NVML
slurmd: Gres Name=gpu Type=nvidia_rtx_a6000 Count=1 Index=1 ID=7696487 File=/dev/nvidia1 Cores=0-63 CoreCnt=128 Links=4,-1,0,0,0,0,0,0 Flags=HAS_FILE,HAS_TYPE,ENV_NVML
slurmd: Gres Name=gpu Type=nvidia_rtx_a6000 Count=1 Index=2 ID=7696487 File=/dev/nvidia2 Cores=0-63 CoreCnt=128 Links=0,0,-1,4,0,0,0,0 Flags=HAS_FILE,HAS_TYPE,ENV_NVML
slurmd: Gres Name=gpu Type=nvidia_rtx_a6000 Count=1 Index=3 ID=7696487 File=/dev/nvidia3 Cores=0-63 CoreCnt=128 Links=0,0,4,-1,0,0,0,0 Flags=HAS_FILE,HAS_TYPE,ENV_NVML
slurmd: Gres Name=gpu Type=nvidia_rtx_a6000 Count=1 Index=4 ID=7696487 File=/dev/nvidia4 Cores=64-127 CoreCnt=128 Links=0,0,0,0,-1,4,0,0 Flags=HAS_FILE,HAS_TYPE,ENV_NVML
slurmd: Gres Name=gpu Type=nvidia_rtx_a6000 Count=1 Index=5 ID=7696487 File=/dev/nvidia5 Cores=64-127 CoreCnt=128 Links=0,0,0,0,4,-1,0,0 Flags=HAS_FILE,HAS_TYPE,ENV_NVML
slurmd: Gres Name=gpu Type=nvidia_rtx_a6000 Count=1 Index=6 ID=7696487 File=/dev/nvidia6 Cores=64-127 CoreCnt=128 Links=0,0,0,0,0,0,-1,4 Flags=HAS_FILE,HAS_TYPE,ENV_NVML
slurmd: Gres Name=gpu Type=nvidia_rtx_a6000 Count=1 Index=7 ID=7696487 File=/dev/nvidia7 Cores=64-127 CoreCnt=128 Links=0,0,0,0,0,0,4,-1 Flags=HAS_FILE,HAS_TYPE,ENV_NVML
slurmd: Gres Name=VramPerGpu Type=(null) Count=51539607552 ID=3033812246 Links=(null) Flags=CountOnly
```
What I understand: when you request two GPU, Slurm will try to allocate two GPUs connected together.
HTH
Best
Yann


## Post 5 by @maciej.falkiewicz (2024-11-29T14:39:39.908Z)

Thank you very much for the reply and the solution!


## Post 6 by @Yann.Sagon (2024-11-29T14:58:45.939Z)

Feedbacks welcome as we have no experience with those parameters.
