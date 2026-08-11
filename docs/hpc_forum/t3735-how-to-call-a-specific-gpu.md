# How to call a specific GPU

- Source: https://hpc-community.unige.ch/t/3735

- Created: 2024-11-20T10:06:26.010Z

- Posts: 5

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Denis.Mongin (2024-11-20T10:06:26.045Z)

Dear team
I am wondering how to call a specific GPU in the bash call.
Here is my situation: Is it possible to actually call the RTX 3080 GPU on baobab ?
In the doc, you can specify the VRAM per GPU, the precision, so I most of the time do:
```
#SBATCH --time=00:05:00
#SBATCH --gpus-per-task=2
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1
#SBATCH --partition=shared-gpu
#SBATCH --gres=VramPerGpu:10G
#SBATCH --mem=20G
#SBATCH --constraint="COMPUTE_CAPABILITY_8_6"
```
To have 2 GPU RTX 3080 of 10G.
But I noticed that sometime I have a a single RTX 3090 instead.
Is there any option to call a precise model of GPU ? I did not find in the doc.
I tried
```
#SBATCH --constraint=COMPUTE_MODEL_RTX_3080_10G
```
and
```
#SBATCH --constraint=COMPUTE_MODEL_RTX_3080
```
But it did not work:
```
sbatch: error: Batch job submission failed: Invalid feature specification
```
Second question: is there a place where I can have an idea of the disponibility of the GPU (to adapt my call) ?
Thank you for your help
Denis


## Post 2 by @Yann.Sagon (2024-11-20T15:19:39.572Z)

Denis.Mongin:
> To have 2 GPU RTX 3080 of 10G.
> But I noticed that sometime I have a a single RTX 3090 instead.
Dear @Denis.Mongin[@Denis.Mongin](https://hpc-community.unige.ch/u/denis.mongin) this shouldn’t happen. You should get two GPUs with compute capability of 8.6 and at least 10G of GPU ram.
Do you have a job id where this has happened or a log?


## Post 3 by @Denis.Mongin (2024-11-21T10:25:06.148Z)

Dear @Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon) , you are right, I had two RTX 3090.
Yet is there a way to call a specific GPU ?
And were can I get a glance of GPU utilisation, to adapt my sbatch call ?
The reason:
I actually want to launch job where my interest is the total Vram, that I can share between N GPU. I am trying to get the most disponibility on baobab, to be able to launch multiple jobs, for example ~100 jobs with lets say 30G of Vram.
3x10G would work, but if the RTX 3080 are not available, then I jobs use 3x25, or 3x48 (on node 48), which is a pity and way to much for my job (I would rather use 1x48 for the jobs, and thus get more jobs done for same dispo of Gpus).
Thus my desire to call a specific GPU: In my example, per job, I would call 3x10 RTX 3080 if those are available, or 2x25G RTX A5500 or RTX 3090 if these are the GPU available, or 1x48 RTX A6000 if those are available.
But
```
#SBATCH --constraint=COMPUTE_MODEL_RTX_3080_10G
```
Does not work, although you seemed to propose this here:
How to ask for a specific GPU?[How to ask for a specific GPU?](https://hpc-community.unige.ch/t/how-to-ask-for-a-specific-gpu/2116)
> Hello, 
Now that there are RTX 3080 and 3090, when I ask for an ampere card I can get these GPUs. It is really annoying, is it possible to ask only for A100 gpus ? 
There is probably a solution but I don’t know SLURM. 
I would try to use --nodelist=gpu[020,022,027,028] but I don’t want to affect negatively the cluster then I prefer to ask first. 
Thank you, 
Ludovic


## Post 4 by @Yann.Sagon (2024-11-21T15:35:30.034Z)

I would like to better understand your workflow to give you a more accurate answer.
Denis.Mongin:
> I actually want to launch job where my interest is the total Vram, that I can share between N GPU. I am trying to get the most disponibility on baobab, to be able to launch multiple jobs, for example ~100 jobs with lets say 30G of Vram.
If I understand correctly, each of your job is able to use let say three GPUs with 10GB RAM at the same time and it “sees” 30GB of GPU RAM? Is that correct?
When you are talking about 100 jobs: are you referring to 100 sbatch instances or job array or one sbatch with 100 job inside using a for loop for example?
Denis.Mongin:
> But
```
#SBATCH --constraint=COMPUTE_MODEL_RTX_3080_10G
```
> Does not work, although you seemed to propose this here:
This was never implemented as at the end the usage of `VramPerGpu` did the trick for the user.
As this may be needed for your use case, I have implemented that this afternoon on Baobab. You can check on the documentation[documentation](https://doc.eresearch.unige.ch/hpc/hpc_clusters#gpus_models_on_the_clusters) the constraint name to use to target a specific GPU.
Denis.Mongin:
> And were can I get a glance of GPU utilisation, to adapt my sbatch call ?
You can use `sinfo` to have a realtime overview:
```
(baobab)-[root@login1 ~]$ sinfo -p shared-gpu -N --Format="NodeList:10,CPUsState,Memory:10,StateLong:10,Gres:45,GresUsed:40,Features:80"
NODELIST  CPUS(A/I/O/T)       MEMORY    STATE     GRES                                         GRES_USED                               AVAIL_FEATURES
gpu002    9/3/0/12            257000    mixed     gpu:titan:6,VramPerGpu:no_consume:12G        gpu:titan:6(IDX:0-5),VramPerGpu:0       E5-2643V3,V5,COMPUTE_CAPABILITY_6_1,COMPUTE_TYPE_TITAN,SIMPLE_PRECISION_GPU
gpu004    12/8/0/20           128000    mixed     gpu:pascal:6,VramPerGpu:no_consume:12G       gpu:pascal:3(IDX:0-2),VramPerGpu:0      E5-2630V4,V6,COMPUTE_CAPABILITY_6_0,COMPUTE_TYPE_PASCAL,DOUBLE_PRECISION_GPU
gpu005    16/4/0/20           128000    mixed     gpu:pascal:5,VramPerGpu:no_consume:12G       gpu:pascal:2(IDX:0-1),VramPerGpu:0      E5-2630V4,V6,COMPUTE_CAPABILITY_6_0,COMPUTE_TYPE_PASCAL,DOUBLE_PRECISION_GPU
[...]
```
In the output you can see have a look a the columns GRES and GRES_USED that will display the number of GPUs vs the number of allocated GPUs.
One more thing: if you specify that you want a GPU with 10GB of RAM and there’s none available, you’ll get a model with more RAM, but the model will be chosen based on its “weight”. The heavier it is, the less likely you’ll get it. Check our documentation to see the weight associated with each GPU.


## Post 5 by @Denis.Mongin (2024-11-22T07:17:26.026Z)

Dear @Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon)
Thanks for the detailed answer !
> If I understand correctly, each of your job is able to use let say three GPUs with 10GB RAM at the same time and it “sees” 30GB of GPU RAM? Is that correct?
Yes! I am using llms locally, and I load the model on various GPU.
> When you are talking about 100 jobs: are you referring to 100 sbatch instances or job array or one sbatch with 100 job inside using a for loop for example?
Speaking about job array here.
> One more thing: if you specify that you want a GPU with 10GB of RAM and there’s none available, you’ll get a model with more RAM, but the model will be chosen based on its “weight”. The heavier it is, the less likely you’ll get it. Check our documentation to see the weight associated with each GPU.
Yes, and that is why I would like to avoid that: I don’t want to have a use multiple “heavy” GPU when light ones would be enough: it penalises me for future jobs (longer queues).
Thank you so much for the help!
Denis
