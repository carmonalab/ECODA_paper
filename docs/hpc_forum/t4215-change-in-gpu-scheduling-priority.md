# Change in GPU scheduling priority

- Source: https://hpc-community.unige.ch/t/4215

- Created: 2026-02-06T07:59:33.462Z

- Posts: 9

- Category: 6

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yann.Sagon (2026-02-06T07:59:33.546Z)

Dear uses,
We have made a change to the way GPUs are allocated.
Previously, if you requested a GPU without specifying any constraints, you might have been allocated either a high- or low-end GPU.
To prevent this from happening, we have assigned a weight[weight](https://slurm.schedmd.com/slurm.conf.html#OPT_Weight) to each GPU.  Low-end GPUs have a low weight and high-end GPUs have a high weight. When you request a GPU, those with a lower weight are prioritised. It was already the case for some old GPUs, but now it is up to date and based on the GPU billing value.
This has two benefits: you’ll get lower-cost GPUs if you don’t need a more powerful one, and we can reserve the high-end GPUs for those who need them.
Best regards


## Post 2 by @Ludovic.Dumoulin (2026-02-06T10:17:50.142Z)

Hello,
Thanks for the update. My job array 6695592_[71-128] is currently pending with (BadConstraints).
Is this a temporary state that will resolve on its own, or is a re-submission required following the system changes?
Best,
Ludovic Dumoulin


## Post 3 by @Ludovic.Dumoulin (2026-02-06T10:32:51.861Z)

I saw the update on the doc for the new slurm name, would it be possible to use the same format for all the node ?
A100
Ampere
40GB
8.0
nvidia_a100_40gb_pcie
ampere
3
gpu[027]
A100
Ampere
40GB
8.0
nvidia_a100-pcie-40gb
ampere
6
gpu[022]
A100
Ampere
40GB
8.0
nvidia_a100-pcie-40gb
ampere
1
gpu[028]
A100
Ampere
40GB
8.0
nvidia_a100-pcie-40gb
ampere
4
gpu[020,030-031]
A100
Ampere
80GB
8.0
nvidia_a100-pcie-80gb
ampere
4
gpu[029]
A100
Ampere
80GB
8.0
nvidia_a100-pcie-80gb
ampere
3
gpu[032-033]
A100
Ampere
80GB
8.0
nvidia_a100-pcie-80gb
ampere
the node GPU022 and GPU027 have the same gpus but one needs to write nvidia_a100_40gb_pcie and nvidia_a100_pcie_40gb ?
There is a similar issue on Bamboo
A100
Ampere
80GB
8.0
nvidia_a100_80gb_pcie
4
gpu[003]
YES
H100
Hopper
94GB
9.0
nvidia_h100_nvl
1
gpu[004]
NO
H200
Hopper
144GB
9.0
nvidia_h200_nvl
4
gpu[005]
NO
H200
Hopper
144GB
9.0
nvidia_h200_nvl
4
gpu[006]
YES
the GPU003 is also 80gb_pcie and not pcie_80gb.
In the doc there is still the old constraints:
A100
40GB
nvidia_a100-pcie-40gb
ampere
COMPUTE_TYPE_AMPERE
COMPUTE_CAPABILITY_8_0
11.0
DOUBLE_PRECISION_GPU
COMPUTE_MODEL_A100_40G
60
A100
80GB
nvidia_a100_80gb_pcie
ampere
COMPUTE_TYPE_AMPERE
COMPUTE_CAPABILITY_8_0
11.0
DOUBLE_PRECISION_GPU
COMPUTE_MODEL_A100_80G
70
then I don’t get why my jobs are not running:
`#!/bin/env bash`
`#SBATCH --array=1-128%40`
`#SBATCH --partition=private-kruse-gpu,shared-gpu`
`#SBATCH --time=0-12:00:00`
`#SBATCH --output=%J.out`
`#SBATCH --mem=8000`
`#SBATCH --gpus=1`
`#SBATCH --constraint=DOUBLE_PRECISION_GPU,COMPUTE_TYPE_AMPERE`
Thank you for your help.


## Post 4 by @Yann.Sagon (2026-02-06T11:02:19.564Z)

Well, this is a mess:)
This happens because we manually updated by copy-paste the doc, not very efficient!
We are using slurm autodetect for gpus, thus the GPU name are automatically set and not very coherent. Anyway, the good news is that for a single GPU model, the name is always the same!
You can list all the GPU models with their correct name like this:
```
sinfo -N -o "%N %G"
```
We’ll soon autogenerate the documentation with a script to be more accurate.
Then in your case, you probably wanted to request an A100 with 40 or 80GB? This is another change we did, not announced, because it is for our internal purpose, but, unfortunately I just figured out that the “DOUBLE_PRECISION_GPU” Feature (to be used with `--constraint`) is missing. We’ll correct that this afternoon I hope.
In the meantime what you can do is to identify every compute node which has an A100 and specify this nodelist.
Sorry for the inconvenience, I’ll keep you posted.


## Post 5 by @Yann.Sagon (2026-02-06T12:18:20.214Z)

Me again!
We now expose the GPU model in each node’s `Feature` field, allowing users to request jobs matching several specific GPU models. Previously, filtering was coarse (e.g., “Ampere with ≥40GB”), which could result in receiving an unintended GPU type such as an RTX A6000. Since users can now combine multiple GPU models in a single constraint expression, we removed the outdated `DOUBLE_PRECISION` constraint. Researchers should now explicitly select the GPU models they want.
Example:
```
srun --partition=shared-gpu --gpus 1 --constraint "nvidia_a100-pcie-40gb|nvidia_a100_80gb_pcie" nvidia-smi
```


## Post 6 by @Yann.Sagon (2026-02-06T12:35:52.290Z)

I’ve updated your pending jobs to use the A100 40GB or 80GB, they are running now!


## Post 7 by @Ludovic.Dumoulin (2026-02-06T14:24:54.160Z)

Thank you,
Yann.Sagon:
> Researchers should now explicitly select the GPU models they want.
I agree, I also think it is better like this.
Before I update our lab’s internal job submission interface to match these changes, I have a few suggestions regarding the naming conventions. Currently, there are some inconsistencies that make automated scripting difficult for users:
- Vendor Prefix Consistency: Most GPUs start with `nvidia_`, but the P100s are still named `tesla_p100-pcie-12gb`. Standardizing all of them to start with the vendor name (e.g., `nvidia_p100...`) would make it easier to parse and select GPUs programmatically.
- Naming Format: Most models follow the `vendor_model_interface_vram` pattern, but the A100s use a different order (`nvidia_A100_vram_pcie`). A unified string format across all models would be ideal.
- Relevance of Interface (PCIe/NVL) vs. P2P: Since there isn’t currently a choice between PCIe and NVL for the same model, including these tags might be redundant. However, knowing if Peer-to-Peer (P2P) communication is supported is crucial. For example, node `gpu[006]` on Bamboo supports P2P, while `pu[004-005]`do not, despite similar names. Perhaps a `_p2p` suffix would be more informative than the bus type?
For our lab’s workflow, we typically filter GPUs based on:
- Vendor: (NVIDIA/AMD/Intel) depending on the framework (CUDA, ROCm, OneAPI).
- Model/Compute Capability: (e.g., CC 8.0 for double-precision FMA using Tensor Cores).
- VRAM: Based on the system size.
- P2P Support: Essential for multi-GPU simulations.
I have already submitted the survey, but I wanted to share these specific thoughts on naming schemes. While these reflect our lab’s needs, I believe a more coherent naming convention would benefit the entire user base.
Thanks


## Post 8 by @Yann.Sagon (2026-02-06T14:48:48.376Z)

Ludovic.Dumoulin:
> Vendor Prefix Consistency: Most GPUs start with `nvidia_`, but the P100s are still named `tesla_p100-pcie-12gb`. Standardizing all of them to start with the vendor name (e.g., `nvidia_p100...`) would make it easier to parse and select GPUs programmatically.
Yann.Sagon:
> We are using slurm autodetect for gpus, thus the GPU name are automatically set and not very coherent.
As I said, we rely on autodetect naming from nvml/slurm Slurm Workload Manager - Generic Resource (GRES) Scheduling[Slurm Workload Manager - Generic Resource (GRES) Scheduling](https://slurm.schedmd.com/gres.html#AutoDetect) and thus have no control over the naming scheme. But yes I fully agree that is is incoherent! There is already a case[case](https://support.schedmd.com/show_bug.cgi?id=21741) opened at schedmd. But if they do something, you’ll have to update your configuration again:). I’ll post in this issue.
Ludovic.Dumoulin:
> However, knowing if Peer-to-Peer (P2P) communication is supported is crucial. For example, node `gpu[006]` on Bamboo supports P2P, while `pu[004-005]`do not, despite similar names. Perhaps a `_p2p` suffix would be more informative than the bus type?
As the word ‘crucial’ is present in the text, it is clear that it was generated by ChatGPT. :smile: no offense.
One option would be to add a flag, such as “NVLINK”. However, the issue then arises that there are no guarantees that, if you request two GPUs on the compute node, they will be connected together. Why is that? If the GPUs are connected in pairs, perhaps the other GPUs are already in use and there are no free pairs. I have no idea how to ensure that two GPUs are linked together, and as far as I know, SLURM doesn’t have a specific flag to enforce this request.
edit: I confirm slurm hasn’t yet a way to force to have a link between gpus. 15995 – Allocate pairs of GPUs with NVLINK[15995 – Allocate pairs of GPUs with NVLINK](https://support.schedmd.com/show_bug.cgi?id=15995)


## Post 9 by @Ludovic.Dumoulin (2026-02-09T09:24:07.766Z)

Yann.Sagon:
> As I said, we rely on autodetect naming from nvml/slurm Slurm Workload Manager - Generic Resource (GRES) Scheduling[Slurm Workload Manager - Generic Resource (GRES) Scheduling](https://slurm.schedmd.com/gres.html#AutoDetect) and thus have no control over the naming scheme.
Oh okay, I didn’t quite understand. I thought you were relying on the default Slurm naming for now, but that you might update the names later on.
Yann.Sagon:
> But if they do something, you’ll have to update your configuration again:). I’ll post in this issue.
Thank you. Then, I will probably go with parsing**`sinfo`** to directly select the desired nodes each time a job is submitted via our interface.
Yann.Sagon:
> As the word ‘crucial’ is present in the text, it is clear that it was generated by ChatGPT. :smile: no offense.
Yes almost all my messages are rewritten by Gemini/ChatGPT :slight_smile: . It helps turn my rough English into something much more readable and fluid.
Yann.Sagon:
> SLURM doesn’t have a specific flag to enforce this request
Oh, I hope it will be updated !
Thank you for your reply. I have a better understanding of the constraints now.
