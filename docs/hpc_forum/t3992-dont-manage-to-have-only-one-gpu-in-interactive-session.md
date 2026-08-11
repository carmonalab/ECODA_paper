# Don't manage to have only one GPU in interactive session

- Source: https://hpc-community.unige.ch/t/3992

- Created: 2025-06-25T08:37:16.722Z

- Posts: 10

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Denis.Mongin (2025-06-25T08:37:16.765Z)

## Primary informations
Username: mongin
Cluster: baobab
## Description
I would like to perform tests with a specific GPU on baobab.
I manage, thanks to @Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon) (How to call a specific GPU - #5 by Denis.Mongin[How to call a specific GPU - #5 by Denis.Mongin](https://hpc-community.unige.ch/t/how-to-call-a-specific-gpu/3735/5)) to ask for a specific GPU, but my problem is I don’t get the number of GPU I want, I get too many GPU.
Here is the command I use:
```
salloc --gpus=1 --partition=shared-gpu  --time=01:00:00 --constraint=COMPUTE_MODEL_RTX_A6000_48G  --mem=20G --ntasks=1 --cpus-per-task=1
```
I ask here for 1 GPU.
But when checking the memory:
```
>>> get_all_mem_info()
2025-06-25 10h-35m-47s : GPU NVIDIA RTX A6000
2025-06-25 10h-35m-47s : memory info: used 7.69G; total 51.16G; free 43.47G
2025-06-25 10h-35m-47s : GPU NVIDIA RTX A6000
2025-06-25 10h-35m-47s : memory info: used 9.06G; total 51.16G; free 42.1G
2025-06-25 10h-35m-47s : GPU NVIDIA RTX A6000
2025-06-25 10h-35m-47s : memory info: used 9.06G; total 51.16G; free 42.1G
2025-06-25 10h-35m-47s : GPU NVIDIA RTX A6000
2025-06-25 10h-35m-47s : memory info: used 9.06G; total 51.16G; free 42.1G
2025-06-25 10h-35m-47s : GPU NVIDIA RTX A6000
2025-06-25 10h-35m-47s : memory info: used 9.06G; total 51.16G; free 42.1G
2025-06-25 10h-35m-47s : GPU NVIDIA RTX A6000
2025-06-25 10h-35m-47s : memory info: used 9.06G; total 51.16G; free 42.1G
2025-06-25 10h-35m-47s : GPU NVIDIA RTX A6000
2025-06-25 10h-35m-47s : memory info: used 9.06G; total 51.16G; free 42.1G
2025-06-25 10h-35m-47s : GPU NVIDIA RTX A6000
2025-06-25 10h-35m-47s : memory info: used 5.74G; total 51.16G; free 45.42G
```
I get the 8 GPU of the cluster.
How do I restrict so that I only get one, to mimic a situation where I only have access to one GPU ?


## Post 2 by @Yann.Sagon (2025-06-25T08:52:58.678Z)

Hi,
can you please share the code you do in python to have this output so I can reproduce the issue?


## Post 3 by @Denis.Mongin (2025-06-25T09:51:23.865Z)

```
def log_mess(message):
   plouf = datetime.datetime.now() 
   print(plouf.strftime("%Y-%m-%d %Hh-%Mm-%Ss") + " : " + message)

def get_mem_info(i):
  meminfo = torch.cuda.mem_get_info(i)
  free = np.round(meminfo[0]/1e9,2)
  total = np.round(meminfo[1]/1e9,2)
  used = np.round(total - free,2)
  log_mess("memory info: used " + str(used) + "G; total " + str(total)+ "G; free " + str(free) + "G")

# print devices

def get_all_mem_info():
   for i in range(torch.cuda.device_count()):
      log_mess("GPU " + torch.cuda.get_device_name(i))
      get_mem_info(i)
```


## Post 4 by @Yann.Sagon (2025-06-25T09:52:21.948Z)

Hi, thanks for the code, in the meantime we figured out what the issue is, we are fixing it right now.


## Post 5 by @Denis.Mongin (2025-06-25T09:52:33.966Z)

I just tested using sbatch, and I don’t have the problem in this situation:
something like:
```
#SBATCH --job-name simple_load
#SBATCH --time=05:00:00
#SBATCH --gpus-per-task=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --partition=shared-gpu
#SBATCH --constraint=COMPUTE_MODEL_RTX_A6000_48G
#SBATCH --mem=20G
#SBATCH --array=1
```
gives only one single RTX A6000, using the same python function.


## Post 6 by @Yann.Sagon (2025-06-25T09:53:23.479Z)

The issue is when you request the resources and connect to the compute node using ssh.


## Post 7 by @Denis.Mongin (2025-06-25T09:53:35.843Z)

Thanks for the reactivity, you are the best


## Post 8 by @Denis.Mongin (2025-06-25T10:44:45.359Z)

yes, exactly, when I am using interactive session with ssh, it looks like I have access to all the GPU of the node, not only to the ressources I asked


## Post 9 by @Denis.Mongin (2025-06-25T10:51:35.654Z)

The same happens using ` nvidia-smi` command:
```
(baobab)-[mongin@gpu048 ~]$ nvidia-smi
Wed Jun 25 12:48:57 2025
+-----------------------------------------------------------------------------------------+
| NVIDIA-SMI 575.51.03              Driver Version: 575.51.03      CUDA Version: 12.9     |
|-----------------------------------------+------------------------+----------------------+
| GPU  Name                 Persistence-M | Bus-Id          Disp.A | Volatile Uncorr. ECC |
| Fan  Temp   Perf          Pwr:Usage/Cap |           Memory-Usage | GPU-Util  Compute M. |
|                                         |                        |               MIG M. |
|=========================================+========================+======================|
|   0  NVIDIA RTX A6000               On  |   00000000:01:00.0 Off |                  Off |
| 30%   29C    P8             22W /  300W |       1MiB /  49140MiB |      0%      Default |
|                                         |                        |                  N/A |
+-----------------------------------------+------------------------+----------------------+
|   1  NVIDIA RTX A6000               On  |   00000000:22:00.0 Off |                  Off |
| 30%   27C    P8             28W /  300W |       1MiB /  49140MiB |      0%      Default |
|                                         |                        |                  N/A |
+-----------------------------------------+------------------------+----------------------+
|   2  NVIDIA RTX A6000               On  |   00000000:41:00.0 Off |                  Off |
| 31%   60C    P2            191W /  300W |   40738MiB /  49140MiB |     51%      Default |
|                                         |                        |                  N/A |
+-----------------------------------------+------------------------+----------------------+
|   3  NVIDIA RTX A6000               On  |   00000000:61:00.0 Off |                  Off |
| 30%   27C    P8             18W /  300W |       1MiB /  49140MiB |      0%      Default |
|                                         |                        |                  N/A |
+-----------------------------------------+------------------------+----------------------+
|   4  NVIDIA RTX A6000               On  |   00000000:81:00.0 Off |                  Off |
| 30%   28C    P8             19W /  300W |       1MiB /  49140MiB |      0%      Default |
|                                         |                        |                  N/A |
+-----------------------------------------+------------------------+----------------------+
|   5  NVIDIA RTX A6000               On  |   00000000:A1:00.0 Off |                  Off |
| 30%   26C    P8             24W /  300W |       1MiB /  49140MiB |      0%      Default |
|                                         |                        |                  N/A |
+-----------------------------------------+------------------------+----------------------+
|   6  NVIDIA RTX A6000               On  |   00000000:C1:00.0 Off |                  Off |
| 30%   28C    P8             24W /  300W |       1MiB /  49140MiB |      0%      Default |
|                                         |                        |                  N/A |
+-----------------------------------------+------------------------+----------------------+
|   7  NVIDIA RTX A6000               On  |   00000000:E1:00.0 Off |                  Off |
| 30%   27C    P8             23W /  300W |    1108MiB /  49140MiB |      0%      Default |
|                                         |                        |                  N/A |
+-----------------------------------------+------------------------+----------------------+

+-----------------------------------------------------------------------------------------+
| Processes:                                                                              |
|  GPU   GI   CI              PID   Type   Process name                        GPU Memory |
|        ID   ID                                                               Usage      |
|=========================================================================================|
|    2   N/A  N/A         1699106      C   ...ab_python_env_LLM3/bin/python      40730MiB |
|    7   N/A  N/A         1672598      C   python3                                1098MiB |
+-----------------------------------------------------------------------------------------+
```


## Post 10 by @Yann.Sagon (2025-06-25T11:48:53.402Z)

This has been fixed in the three clusters, BUT if a user requested a GPU before the fix and connected to the compute node using SSH, and is using a different GPU to the one allocated (a lot of “and”), this could cause an issue until their job finishes. Unfortunately, it is not easy to detect this situation.
