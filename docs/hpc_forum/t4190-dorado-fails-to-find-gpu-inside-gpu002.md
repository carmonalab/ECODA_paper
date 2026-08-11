# Dorado fails to find GPU inside GPU002

- Source: https://hpc-community.unige.ch/t/4190

- Created: 2026-01-13T09:05:58.334Z

- Tags: yggdrasil

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Juan.Diazmiyar (2026-01-13T09:05:58.389Z)

## Primary informations
Username: diazmiya
Cluster: Yggdrasil
## Description
When I run `dorado basecall` on GPU002, it does not find the GPU and exits on error. The same code works on any other GPU in the public-gpu partition.
## Steps to Reproduce
You will need to download a dorado model first with `dorado download --model dna_r10.4.1_e8.2_400bps_sup@v5.0.0 --models-directory $MODELDIR`
Then, run `dorado basecall` on any pod5 file:
```
MODELDIR=results/PAO_SHIFT_SMK/dorado_model/dna_r10.4.1_e8.2_400bps_sup@v5.0.0

POD5=${SCRATCH}/pao-shift-20250909/20250909_1455_MN24748_FBD81765_9dfd643d/pod5/FBD81765_9dfd643d_ec484651_10.pod5

srun --nodelist=gpu002 --gres=gpu:1 --partition=shared-gpu dorado basecaller \
     --device cuda:0 --kit-name SQK-NBD114-24 \
     --trim none $MODELDIR $POD5 > test_gpu002cuda0.bam
```
## Expected Result
Basecalling starts (log from a run on gpu004):
```
[2026-01-13 09:50:40.939] [info] Running: "basecaller" "--device" "cuda:0" "--kit-name" "SQK-NBD114-24" "--trim" "none" "results/PAO_SHIFT_SMK/dorado_model/dna_r10.4.1_e8.2_400bps_sup@v5.0.0" "/srv/beegfs/scratch/users/d/diazmiya/pao-shift-20250909/20250909_1455_MN24748_FBD81765_9dfd643d/pod5/FBD81765_9dfd643d_ec484651_10.pod5"
[2026-01-13 09:50:41.587] [info] > Creating basecall pipeline
[2026-01-13 09:50:41.785] [info] Using CUDA devices:
[2026-01-13 09:50:41.785] [info] cuda:0 - NVIDIA TITAN RTX
[2026-01-13 09:50:43.544] [info] Calculating optimized batch size for GPU "NVIDIA TITAN RTX" and model dna_r10.4.1_e8.2_400bps_sup@v5.0.0. Full benchmarking will run for this device, which may take some time.
[2026-01-13 09:50:53.104] [info] cuda:0 using chunk size 12288, batch size 288
[2026-01-13 09:50:54.267] [info] cuda:0 using chunk size 6144, batch size 512
[...]
```
## Actual Result
Error looks like this:
```
[2026-01-13 09:47:52.579] [info] Running: "basecaller" "--device" "cuda:0" "--kit-name" "SQK-NBD114-24" "--trim" "none" "results/PAO_SHIFT_SMK/dorado_model/dna_r10.4.1_e8.2_400bps_sup@v5.0.0" "/srv/beegfs/scratch/users/d/diazmiya/pao-shift-20250909/20250909_1455_MN24748_FBD81765_9dfd643d/pod5/FBD81765_9dfd643d_ec484651_10.pod5"
[2026-01-13 09:47:52.790] [info] > Creating basecall pipeline
[2026-01-13 09:47:52.877] [info] Using CUDA devices:
[2026-01-13 09:47:52.878] [info] cuda:0 -  fX�
[2026-01-13 09:47:52.883] [error] CUDA error: invalid device ordinal
GPU device may be out of range, do you have enough GPUs?
CUDA kernel errors might be asynchronously reported at some other API call, so the stacktrace below might be incorrect.
For debugging consider passing CUDA_LAUNCH_BLOCKING=1
Device-side assertions were explicitly omitted for this error check; the error probably arose while initializing the DSA handlers.
Exception raised from c10_cuda_check_implementation at /builds/machine-learning/torch-builds/pytorch/c10/cuda/CUDAException.cpp:44 (most recent call first):
frame #0: c10::Error::Error(c10::SourceLocation, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >) + 0x88 (0x7f865cf5a7d8 in /opt/ebsofts/dorado/1.3.0/lib/libc10.so)
frame #1: <unknown function> + 0x5e8aa (0x7f865d0298aa in /opt/ebsofts/dorado/1.3.0/lib/libc10_cuda.so)
frame #2: c10::cuda::c10_cuda_check_implementation(int, char const*, char const*, int, bool) + 0x1c9 (0x7f865d029579 in /opt/ebsofts/dorado/1.3.0/lib/libc10_cuda.so)
frame #3: c10::cuda::CUDAKernelLaunchRegistry::CUDAKernelLaunchRegistry() + 0xba (0x7f865d02850a in /opt/ebsofts/dorado/1.3.0/lib/libc10_cuda.so)
frame #4: c10::cuda::CUDAKernelLaunchRegistry::get_singleton_ref() + 0x4a (0x7f865d02879a in /opt/ebsofts/dorado/1.3.0/lib/libc10_cuda.so)
frame #5: c10::cuda::c10_cuda_check_implementation(int, char const*, char const*, int, bool) + 0x45 (0x7f865d0293f5 in /opt/ebsofts/dorado/1.3.0/lib/libc10_cuda.so)
frame #6: c10::cuda::ExchangeDevice(signed char) + 0x9a (0x7f865d029c5a in /opt/ebsofts/dorado/1.3.0/lib/libc10_cuda.so)
frame #7: /opt/ebsofts/dorado/1.3.0/bin/dorado() [0xea92aa]
frame #8: /opt/ebsofts/dorado/1.3.0/bin/dorado() [0xfa0211]
frame #9: /opt/ebsofts/dorado/1.3.0/bin/dorado() [0xeba4db]
frame #10: /opt/ebsofts/dorado/1.3.0/bin/dorado() [0xeb7bb7]
frame #11: /opt/ebsofts/dorado/1.3.0/bin/dorado() [0x57985b]
frame #12: <unknown function> + 0x8f3a8 (0x7f865568f3a8 in /lib64/libc.so.6)
frame #13: /opt/ebsofts/dorado/1.3.0/bin/dorado() [0xeb812a]
frame #14: /opt/ebsofts/dorado/1.3.0/bin/dorado() [0x57af03]
frame #15: <unknown function> + 0xc2b23 (0x7f8655aa6b23 in /opt/ebsofts/dorado/1.3.0/lib/libstdc++.so.6)
frame #16: <unknown function> + 0x8a19a (0x7f865568a19a in /lib64/libc.so.6)
frame #17: <unknown function> + 0x10f100 (0x7f865570f100 in /lib64/libc.so.6)
```
If you remove the `--device` argument or use `cuda:all`(which auto-detects the available gpus), it produces the same error. I ran `srun --nodelist=gpu002 --gres=gpu:1 --partition=shared-gpu nvidia-smi`and it says that the GPU is alive and well and mapped to the “0” position.


## Post 2 by @Yann.Sagon (2026-01-22T14:08:06.977Z)

@Juan.Diazmiyar[@Juan.Diazmiyar](https://hpc-community.unige.ch/u/juan.diazmiyar)
One of the GPUs was faulty and we are now preventing jobs to use it. It should be working now.
Best regards
Yann
