# TensorFlow not recognising GPU

- Source: https://hpc-community.unige.ch/t/3743

- Created: 2024-11-26T09:41:08.345Z

- Tags: yggdrasil

- Posts: 3

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yuzheng.Kang (2024-11-26T09:41:08.393Z)

Hi,
I am trying to run TensorFlow with a GPU node, turns out it seems my TensorFlow cannot recognising the GPU. Could anyone help me with it?
So here is the .sh script where I request a GPU node from Yggdrasil
```
#!/bin/sh
#SBATCH --job-name=cosmopower         # Job name
#SBATCH --partition=shared-gpu             # Partition (queue) name
#SBATCH --nodes=1                   # Number of nodes
#SBATCH --ntasks=1                  # Number of tasks (processes)
#SBATCH --gres=gpu:1                # Request 1 GPU
#SBATCH --time=00:15:00             # Time limit (hh:mm:ss)
#SBATCH --cpus-per-task=1

conda activate cp_env
module load GCC/10.3.0  OpenMPI/4.1.1 TensorFlow/2.6.0
module load cuDNN/8.6 CUDA/11.8.0
nvidia-smi
python gpu_test.py
```
where the gpu_test.py is just print out number of GPUs available like `gpus = tf.config.list_physical_devices('GPU') print(f"Num GPUs Available: {len(gpus)}")`
From the Nvidia-smi, I can see there is a GPU requested, but the output from the python code says “Num GPUs Available: 0”. Does anyone know where is the issue?
```
+-----------------------------------------------------------------------------------------+
| NVIDIA-SMI 560.35.03              Driver Version: 560.35.03      CUDA Version: 12.6     |
|-----------------------------------------+------------------------+----------------------+
| GPU  Name                 Persistence-M | Bus-Id          Disp.A | Volatile Uncorr. ECC |
| Fan  Temp   Perf          Pwr:Usage/Cap |           Memory-Usage | GPU-Util  Compute M. |
|                                         |                        |               MIG M. |
|=========================================+========================+======================|
|   0  NVIDIA TITAN RTX               On  |   00000000:3D:00.0 Off |                  N/A |
| 41%   31C    P8             13W /  280W |       1MiB /  24576MiB |      0%      Default |
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


## Post 2 by @Adrien.Albert (2024-11-26T13:27:35.444Z)

Hi @Yuzheng.Kang[@Yuzheng.Kang](https://hpc-community.unige.ch/u/yuzheng.kang)
You are using a version of tensorflow not built with cuda:
```
(baobab)-[alberta@login1 tensorflow]$ ml spider tensor
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  TensorFlow:
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Description:
      An open-source software library for Machine Intelligence

     Versions:
        TensorFlow/1.7.0-Python-3.6.4
        TensorFlow/1.10.1-Python-2.7.15
        TensorFlow/1.15.0-Python-3.7.4
        TensorFlow/2.0.0-Python-3.7.2
        TensorFlow/2.0.0-Python-3.7.4
        TensorFlow/2.1.0-Python-3.7.4
        TensorFlow/2.4.1
        TensorFlow/2.5.0
        TensorFlow/2.5.3
        TensorFlow/2.6.0-CUDA-11.3.1
        TensorFlow/2.6.0
        TensorFlow/2.7.1-CUDA-11.4.1
        TensorFlow/2.11.0-CUDA-11.7.0
        TensorFlow/2.11.0
```
I tested and it’s working on my side:
### Replacing module compatible with GPU:
```
(baobab)-[alberta@login1 tensorflow]$ cat sbatch.sh 
#!/bin/sh
#SBATCH --job-name=cosmopower         # Job name
#SBATCH --partition=shared-gpu             # Partition (queue) name
#SBATCH --nodes=1                   # Number of nodes
#SBATCH --ntasks=1                  # Number of tasks (processes)
#SBATCH --gres=gpu:1                # Request 1 GPU
#SBATCH --time=00:01:00             # Time limit (hh:mm:ss)
#SBATCH --cpus-per-task=1

conda activate cp_env
module load GCC/11.3.0  OpenMPI/4.1.4 TensorFlow/2.11.0-CUDA-11.7.0
module load cuDNN/8.4.1.50-CUDA-11.7.0
nvidia-smi
python gpu_test.py
```
### Test script:
```
(baobab)-[alberta@login1 tensorflow]$ cat gpu_test.py 
import tensorflow as tf
print("Num GPUs Available: ", len(tf.config.list_physical_devices('GPU')))
```
### Submission and result:
```
(baobab)-[alberta@login1 tensorflow]$ sbatch !$
sbatch sbatch.sh
Submitted batch job 13780842
```
```
(baobab)-[alberta@login1 tensorflow]$ cat slurm-13780842.out
Tue Nov 26 14:21:48 2024       
+-----------------------------------------------------------------------------------------+
| NVIDIA-SMI 565.57.01              Driver Version: 565.57.01      CUDA Version: 12.7     |
|-----------------------------------------+------------------------+----------------------+
| GPU  Name                 Persistence-M | Bus-Id          Disp.A | Volatile Uncorr. ECC |
| Fan  Temp   Perf          Pwr:Usage/Cap |           Memory-Usage | GPU-Util  Compute M. |
|                                         |                        |               MIG M. |
|=========================================+========================+======================|
|   0  NVIDIA TITAN X (Pascal)        On  |   00000000:04:00.0 Off |                  N/A |
| 25%   39C    P8             10W /  250W |       8MiB /  12288MiB |      0%      Default |
|                                         |                        |                  N/A |
+-----------------------------------------+------------------------+----------------------+
                                                                                         
+-----------------------------------------------------------------------------------------+
| Processes:                                                                              |
|  GPU   GI   CI        PID   Type   Process name                              GPU Memory |
|        ID   ID                                                               Usage      |
|=========================================================================================|
|  No running processes found                                                             |
+-----------------------------------------------------------------------------------------+
2024-11-26 14:21:48.921272: I tensorflow/core/platform/cpu_feature_guard.cc:193] This TensorFlow binary is optimized with oneAPI Deep Neural Network Library (oneDNN) to use the following CPU instructions in performance-critical operations:  SSE4.1 SSE4.2 AVX AVX2 FMA
To enable them in other operations, rebuild TensorFlow with the appropriate compiler flags.

Num GPUs Available:  1 <=== Here what we want !
```


## Post 3 by @Yuzheng.Kang (2024-11-26T13:38:13.124Z)

Oh, I did not notice this. Thanks for the help!
