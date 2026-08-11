# SLURM auto-requeueing enabled. Setting signal handlers

- Source: https://hpc-community.unige.ch/t/3319

- Created: 2024-02-18T00:36:28.741Z

- Posts: 4

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Simon.Gabay (2024-02-18T00:36:28.823Z)

Hello,
I am using Kraken[Kraken](https://github.com/mittagessen/kraken), `version 4.3.13`. Prior training, I have compiled my data in binary format[binary format](https://github.com/mittagessen/kraken/blob/main/docs/ketos.rst#binary-datasets). I use a pretty big dataset (a bit more than 2000 pages to train an OCR model), which is why I cannot share it here:
```
(kraken-env) (yggdrasil)-[gabays@login1 ~]$ stat dataset.arrow 
  File: dataset.arrow
  Size: 18305066626	Blocks: 35752084   IO Block: 524288 regular file
```
Before using `sbatch`, I played a bit to find the correct configuration to understand how it works. I did `salloc`:
```
alloc --partition=shared-gpu --time=12:00:00 --gpus=2 --mem=24GB --cpus-per-task=12 --gres=gpu:2,VramPerGpu:24GB
```
Things go fine for the first epoch, then it suddenly stops:
```
(kraken-env) (yggdrasil)-[gabays@gpu002 ~]$ ketos train -f binary -d cuda:0 -B 16 --workers 1 -r 0.0004 -u NFC dataset.arrow --freq 1
scikit-learn version 1.2.2 is not supported. Minimum required version: 0.17. Maximum required version: 1.1.2. Disabling scikit-learn conversion API.
Torch version 2.0.1+cu117 has not been tested with coremltools. You may run into unexpected errors. Torch 2.0.0 is the most recent version that has been tested.
GPU available: True (cuda), used: True
TPU available: False, using: 0 TPU cores
IPU available: False, using: 0 IPUs
HPU available: False, using: 0 HPUs
`Trainer(val_check_interval=1.0)` was configured so validation will run at the end of the training epoch..
LOCAL_RANK: 0 - CUDA_VISIBLE_DEVICES: [0,1]
┏━━━━┳━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃    ┃ Name      ┃ Type                     ┃ Params ┃                 In sizes ┃                Out sizes ┃
┡━━━━╇━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━┩
│ 0  │ val_cer   │ _CharErrorRate           │      0 │                        ? │                        ? │
│ 1  │ val_wer   │ _WordErrorRate           │      0 │                        ? │                        ? │
│ 2  │ net       │ MultiParamSequential     │  4.1 M │  [[1, 1, 120, 400], '?'] │   [[1, 258, 1, 50], '?'] │
│ 3  │ net.C_0   │ ActConv2D                │  1.3 K │  [[1, 1, 120, 400], '?'] │ [[1, 32, 120, 400], '?'] │
│ 4  │ net.Do_1  │ Dropout                  │      0 │ [[1, 32, 120, 400], '?'] │ [[1, 32, 120, 400], '?'] │
│ 5  │ net.Mp_2  │ MaxPool                  │      0 │ [[1, 32, 120, 400], '?'] │  [[1, 32, 60, 200], '?'] │
│ 6  │ net.C_3   │ ActConv2D                │ 40.0 K │  [[1, 32, 60, 200], '?'] │  [[1, 32, 60, 200], '?'] │
│ 7  │ net.Do_4  │ Dropout                  │      0 │  [[1, 32, 60, 200], '?'] │  [[1, 32, 60, 200], '?'] │
│ 8  │ net.Mp_5  │ MaxPool                  │      0 │  [[1, 32, 60, 200], '?'] │  [[1, 32, 30, 100], '?'] │
│ 9  │ net.C_6   │ ActConv2D                │ 55.4 K │  [[1, 32, 30, 100], '?'] │  [[1, 64, 30, 100], '?'] │
│ 10 │ net.Do_7  │ Dropout                  │      0 │  [[1, 64, 30, 100], '?'] │  [[1, 64, 30, 100], '?'] │
│ 11 │ net.Mp_8  │ MaxPool                  │      0 │  [[1, 64, 30, 100], '?'] │   [[1, 64, 15, 50], '?'] │
│ 12 │ net.C_9   │ ActConv2D                │  110 K │   [[1, 64, 15, 50], '?'] │   [[1, 64, 15, 50], '?'] │
│ 13 │ net.Do_10 │ Dropout                  │      0 │   [[1, 64, 15, 50], '?'] │   [[1, 64, 15, 50], '?'] │
│ 14 │ net.S_11  │ Reshape                  │      0 │   [[1, 64, 15, 50], '?'] │   [[1, 960, 1, 50], '?'] │
│ 15 │ net.L_12  │ TransposedSummarizingRNN │  1.9 M │   [[1, 960, 1, 50], '?'] │   [[1, 400, 1, 50], '?'] │
│ 16 │ net.Do_13 │ Dropout                  │      0 │   [[1, 400, 1, 50], '?'] │   [[1, 400, 1, 50], '?'] │
│ 17 │ net.L_14  │ TransposedSummarizingRNN │  963 K │   [[1, 400, 1, 50], '?'] │   [[1, 400, 1, 50], '?'] │
│ 18 │ net.Do_15 │ Dropout                  │      0 │   [[1, 400, 1, 50], '?'] │   [[1, 400, 1, 50], '?'] │
│ 19 │ net.L_16  │ TransposedSummarizingRNN │  963 K │   [[1, 400, 1, 50], '?'] │   [[1, 400, 1, 50], '?'] │
│ 20 │ net.Do_17 │ Dropout                  │      0 │   [[1, 400, 1, 50], '?'] │   [[1, 400, 1, 50], '?'] │
│ 21 │ net.O_18  │ LinSoftmax               │  103 K │   [[1, 400, 1, 50], '?'] │   [[1, 258, 1, 50], '?'] │
└────┴───────────┴──────────────────────────┴────────┴──────────────────────────┴──────────────────────────┘
Trainable params: 4.1 M                                                                                                        
Non-trainable params: 0                                                                                                        
Total params: 4.1 M                                                                                                            
Total estimated model params size (MB): 16                                                                                     
stage 0/∞ ━━━━━━━━━━━━━━━━━━━━━━━━━ 10992/10992 0:31:20 • 0:00:00 6.07it/s val_accuracy: 0.765       early_stopping: 0/10      
                                                                           val_word_accuracy: 0.328  0.76487                   
stage 1/∞ ━━━━━━━╸━━━━━━━━━━━━━━━━━ 3515/10992 0:10:01 • 0:21:18 5.85it/s val_accuracy: 0.765       early_stopping: 0/10       
stage 1/∞ ━━━━━━━━━━━━━━━━━━━━━━━╸━ 10519/10992 0:30:16 • 0:01:22 5.81it/s val_accuracy: 0.765       early_stopping: 0/10      
                                                                           val_word_accuracy: 0.328  0.76487                   Killed
```
I try to restart the training from the first model trained (`--load model_0.mlmodel`), diminishing the size of `-B` from 16 to 8. Same issue, but faster:
```
(kraken-env) (yggdrasil)-[gabays@gpu002 ~]$ ketos train -f binary -d cuda:0 -B 8 --workers 1 -r 0.0003 -u NFC dataset.arrow --freq 1  --load model_0.mlmodel
scikit-learn version 1.2.2 is not supported. Minimum required version: 0.17. Maximum required version: 1.1.2. Disabling scikit-learn conversion API.
Torch version 2.0.1+cu117 has not been tested with coremltools. You may run into unexpected errors. Torch 2.0.0 is the most recent version that has been tested.
GPU available: True (cuda), used: True
TPU available: False, using: 0 TPU cores
IPU available: False, using: 0 IPUs
HPU available: False, using: 0 HPUs
`Trainer(val_check_interval=1.0)` was configured so validation will run at the end of the training epoch..
LOCAL_RANK: 0 - CUDA_VISIBLE_DEVICES: [0,1]
┏━━━━┳━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃    ┃ Name      ┃ Type                     ┃ Params ┃                 In sizes ┃                Out sizes ┃
┡━━━━╇━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━┩
│ 0  │ val_cer   │ _CharErrorRate           │      0 │                        ? │                        ? │
│ 1  │ val_wer   │ _WordErrorRate           │      0 │                        ? │                        ? │
│ 2  │ net       │ MultiParamSequential     │  4.1 M │  [[1, 1, 120, 400], '?'] │   [[1, 258, 1, 50], '?'] │
│ 3  │ net.C_0   │ ActConv2D                │  1.3 K │  [[1, 1, 120, 400], '?'] │ [[1, 32, 120, 400], '?'] │
│ 4  │ net.Do_1  │ Dropout                  │      0 │ [[1, 32, 120, 400], '?'] │ [[1, 32, 120, 400], '?'] │
│ 5  │ net.Mp_2  │ MaxPool                  │      0 │ [[1, 32, 120, 400], '?'] │  [[1, 32, 60, 200], '?'] │
│ 6  │ net.C_3   │ ActConv2D                │ 40.0 K │  [[1, 32, 60, 200], '?'] │  [[1, 32, 60, 200], '?'] │
│ 7  │ net.Do_4  │ Dropout                  │      0 │  [[1, 32, 60, 200], '?'] │  [[1, 32, 60, 200], '?'] │
│ 8  │ net.Mp_5  │ MaxPool                  │      0 │  [[1, 32, 60, 200], '?'] │  [[1, 32, 30, 100], '?'] │
│ 9  │ net.C_6   │ ActConv2D                │ 55.4 K │  [[1, 32, 30, 100], '?'] │  [[1, 64, 30, 100], '?'] │
│ 10 │ net.Do_7  │ Dropout                  │      0 │  [[1, 64, 30, 100], '?'] │  [[1, 64, 30, 100], '?'] │
│ 11 │ net.Mp_8  │ MaxPool                  │      0 │  [[1, 64, 30, 100], '?'] │   [[1, 64, 15, 50], '?'] │
│ 12 │ net.C_9   │ ActConv2D                │  110 K │   [[1, 64, 15, 50], '?'] │   [[1, 64, 15, 50], '?'] │
│ 13 │ net.Do_10 │ Dropout                  │      0 │   [[1, 64, 15, 50], '?'] │   [[1, 64, 15, 50], '?'] │
│ 14 │ net.S_11  │ Reshape                  │      0 │   [[1, 64, 15, 50], '?'] │   [[1, 960, 1, 50], '?'] │
│ 15 │ net.L_12  │ TransposedSummarizingRNN │  1.9 M │   [[1, 960, 1, 50], '?'] │   [[1, 400, 1, 50], '?'] │
│ 16 │ net.Do_13 │ Dropout                  │      0 │   [[1, 400, 1, 50], '?'] │   [[1, 400, 1, 50], '?'] │
│ 17 │ net.L_14  │ TransposedSummarizingRNN │  963 K │   [[1, 400, 1, 50], '?'] │   [[1, 400, 1, 50], '?'] │
│ 18 │ net.Do_15 │ Dropout                  │      0 │   [[1, 400, 1, 50], '?'] │   [[1, 400, 1, 50], '?'] │
│ 19 │ net.L_16  │ TransposedSummarizingRNN │  963 K │   [[1, 400, 1, 50], '?'] │   [[1, 400, 1, 50], '?'] │
│ 20 │ net.Do_17 │ Dropout                  │      0 │   [[1, 400, 1, 50], '?'] │   [[1, 400, 1, 50], '?'] │
│ 21 │ net.O_18  │ LinSoftmax               │  103 K │   [[1, 400, 1, 50], '?'] │   [[1, 258, 1, 50], '?'] │
└────┴───────────┴──────────────────────────┴────────┴──────────────────────────┴──────────────────────────┘
Trainable params: 4.1 M                                                                                                        
Non-trainable params: 0                                                                                                        
Total params: 4.1 M                                                                                                            
Total estimated model params size (MB): 16                                                                                     
stage 0/∞ ━━━━━━━━━━━━━━━━━━━━━━━━━ 21984/21984 0:38:37 • 0:00:00 9.69it/s val_accuracy: 0.855       early_stopping: 0/10      
                                                                           val_word_accuracy: 0.508  0.85462                   
stage 1/∞ ━━━━━━━━━━━━━━━━━━╺━━━━━━ 16021/21984 0:30:43 • 3:05:24 0.54it/s val_accuracy: 0.855       early_stopping: 0/10      
                                                                           val_word_accuracy: 0.508  0.85462                   Killed
```
Then I start training a model with `sbatch` (I keep diminishing the batch size because it is the main reason for `killed`):
```
#!/bin/env bash
#SBATCH --partition=shared-gpu
#SBATCH --time=12:00:00
#SBATCH --gpus=1
#SBATCH --output=kraken-%j.out
#SBATCH --mem=24GB
#SBATCH --cpus-per-task=12
#SBATCH --gres=gpu:1,VramPerGpu:24GB

module load fosscuda/2020b Python/3.8.6
source ~/kraken-env/bin/activate

echo "KETOS training"
srun ketos train -f binary -d cuda:0 -B 2  -r 0.0001 -u NFC dataset.arrow
```
and I get this error:
```
(yggdrasil)-[gabays@login1 ~]$ cat kraken-31170677.out
KETOS training
scikit-learn version 1.2.2 is not supported. Minimum required version: 0.17. Maximum required version: 1.1.2. Disabling scikit-learn conversion API.
Torch version 2.0.1+cu117 has not been tested with coremltools. You may run into unexpected errors. Torch 2.0.0 is the most recent version that has been tested.
GPU available: True (cuda), used: True
TPU available: False, using: 0 TPU cores
IPU available: False, using: 0 IPUs
HPU available: False, using: 0 HPUs
`Trainer(val_check_interval=1.0)` was configured so validation will run at the end of the training epoch..
LOCAL_RANK: 0 - CUDA_VISIBLE_DEVICES: [0,1]
┏━━━━┳━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━┳━━━━━━━━┳━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━┓
┃    ┃ Name      ┃ Type            ┃ Params ┃       In sizes ┃       Out sizes ┃
┡━━━━╇━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━╇━━━━━━━━╇━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━┩
│ 0  │ val_cer   │ _CharErrorRate  │      0 │              ? │               ? │
│ 1  │ val_wer   │ _WordErrorRate  │      0 │              ? │               ? │
│ 2  │ net       │ MultiParamSequ… │  4.1 M │   [[1, 1, 120, │    [[1, 258, 1, │
│    │           │                 │        │     400], '?'] │       50], '?'] │
│ 3  │ net.C_0   │ ActConv2D       │  1.3 K │   [[1, 1, 120, │   [[1, 32, 120, │
│    │           │                 │        │     400], '?'] │      400], '?'] │
│ 4  │ net.Do_1  │ Dropout         │      0 │  [[1, 32, 120, │   [[1, 32, 120, │
│    │           │                 │        │     400], '?'] │      400], '?'] │
│ 5  │ net.Mp_2  │ MaxPool         │      0 │  [[1, 32, 120, │    [[1, 32, 60, │
│    │           │                 │        │     400], '?'] │      200], '?'] │
│ 6  │ net.C_3   │ ActConv2D       │ 40.0 K │   [[1, 32, 60, │    [[1, 32, 60, │
│    │           │                 │        │     200], '?'] │      200], '?'] │
│ 7  │ net.Do_4  │ Dropout         │      0 │   [[1, 32, 60, │    [[1, 32, 60, │
│    │           │                 │        │     200], '?'] │      200], '?'] │
│ 8  │ net.Mp_5  │ MaxPool         │      0 │   [[1, 32, 60, │    [[1, 32, 30, │
│    │           │                 │        │     200], '?'] │      100], '?'] │
│ 9  │ net.C_6   │ ActConv2D       │ 55.4 K │   [[1, 32, 30, │    [[1, 64, 30, │
│    │           │                 │        │     100], '?'] │      100], '?'] │
│ 10 │ net.Do_7  │ Dropout         │      0 │   [[1, 64, 30, │    [[1, 64, 30, │
│    │           │                 │        │     100], '?'] │      100], '?'] │
│ 11 │ net.Mp_8  │ MaxPool         │      0 │   [[1, 64, 30, │    [[1, 64, 15, │
│    │           │                 │        │     100], '?'] │       50], '?'] │
│ 12 │ net.C_9   │ ActConv2D       │  110 K │   [[1, 64, 15, │    [[1, 64, 15, │
│    │           │                 │        │      50], '?'] │       50], '?'] │
│ 13 │ net.Do_10 │ Dropout         │      0 │   [[1, 64, 15, │    [[1, 64, 15, │
│    │           │                 │        │      50], '?'] │       50], '?'] │
│ 14 │ net.S_11  │ Reshape         │      0 │   [[1, 64, 15, │    [[1, 960, 1, │
│    │           │                 │        │      50], '?'] │       50], '?'] │
│ 15 │ net.L_12  │ TransposedSumm… │  1.9 M │   [[1, 960, 1, │    [[1, 400, 1, │
│    │           │                 │        │      50], '?'] │       50], '?'] │
│ 16 │ net.Do_13 │ Dropout         │      0 │   [[1, 400, 1, │    [[1, 400, 1, │
│    │           │                 │        │      50], '?'] │       50], '?'] │
│ 17 │ net.L_14  │ TransposedSumm… │  963 K │   [[1, 400, 1, │    [[1, 400, 1, │
│    │           │                 │        │      50], '?'] │       50], '?'] │
│ 18 │ net.Do_15 │ Dropout         │      0 │   [[1, 400, 1, │    [[1, 400, 1, │
│    │           │                 │        │      50], '?'] │       50], '?'] │
│ 19 │ net.L_16  │ TransposedSumm… │  963 K │   [[1, 400, 1, │    [[1, 400, 1, │
│    │           │                 │        │      50], '?'] │       50], '?'] │
│ 20 │ net.Do_17 │ Dropout         │      0 │   [[1, 400, 1, │    [[1, 400, 1, │
│    │           │                 │        │      50], '?'] │       50], '?'] │
│ 21 │ net.O_18  │ LinSoftmax      │  103 K │   [[1, 400, 1, │    [[1, 258, 1, │
│    │           │                 │        │      50], '?'] │       50], '?'] │
└────┴───────────┴─────────────────┴────────┴────────────────┴─────────────────┘
Trainable params: 4.1 M                                                         
Non-trainable params: 0                                                         
Total params: 4.1 M                                                             
Total estimated model params size (MB): 16                                      
SLURM auto-requeueing enabled. Setting signal handlers.
```
And nothing happens.
Do you have any idea what is happening?


## Post 2 by @Simon.Gabay (2024-02-18T00:37:14.126Z)

The forum does not want me to put more than two links so I add here:
I see Kraken uses Lightning[Lightning](https://kraken.re/4.0/api.html#training) and there are some similar issues on GitHub[GitHub](https://github.com/Lightning-AI/pytorch-lightning/issues/17066). Maybe it helps.
It is possible that it is not a problem with the main memory, but the shared memory, but not clue on how to solve that…


## Post 3 by @Simon.Gabay (2024-02-21T23:49:16.691Z)

The developer of Kraken told me:
It is possible that `cuda:0` is not the GPU that you’ve got assigned by the cluster (which would get your batch killed if something else is running on there already). Usually they set an environment variable as one node might have more than one GPU and masking them out for each individual batch requires some special setup.
But anyway with a batch size of 2 you should be using less than 3Gb of GPU memory. And basically no CPU memory as you don’t have any separate loading workers. So something is awry.
The no data loading workers also means you demanding more CPU cores from the cluster should not have any impact in the crashes at all.
The developer read the doc of Yggdrasil. Second answer:
You can’t just use `cuda:0` on that cluster as there are half a dozen or more GPUs per node. There’s an environment variable called `$CUDA_VISIBLE_DEVICES` that tells you which one is allocated to you.
They’re wrapping the OOM kill message with something (kernel messages look different) so it isn’t directly visible were it comes from. It is possible that your batches just get killed because you’re using somebody else’s GPU and the 50Gb memory allocation thing is just a coincidence.


## Post 4 by @Yann.Sagon (2024-02-27T07:15:29.908Z)

Hello @Simon.Gabay[@Simon.Gabay](https://hpc-community.unige.ch/u/simon.gabay)
You can check yourself if the job was killed for an out of memory reason.
Example:
```
(yggdrasil)-[root@admin1 ~]$ sacct -u gabays -S  2024-02-15 -E 2024-02-19 -o 'Stat,Start,jobid%20' | grep OUT_OF
OUT_OF_ME+ 2024-02-16T04:38:07 31130319.interactive
OUT_OF_ME+ 2024-02-16T20:41:06 31142559.interactive
OUT_OF_ME+ 2024-02-17T15:51:26 31163091.interactive
OUT_OF_ME+ 2024-02-18T10:38:54           31170682.0
OUT_OF_ME+ 2024-02-18T13:00:26 31171825.interactive
```
You can see the detail of your individual jobs like this:
```
(yggdrasil)-[gabays@login1 ~]$ sacct --format=Start,AveCPU,State,MaxRSS,JobID,NodeList,ReqMem --units=G -j 31142559
              Start     AveCPU      State     MaxRSS JobID               NodeList     ReqMem
------------------- ---------- ---------- ---------- ------------ --------------- ----------
2024-02-16T20:41:06               TIMEOUT            31142559              gpu007        24G
2024-02-16T20:41:06   06:44:01 OUT_OF_ME+     95.03G 31142559.in+          gpu007
2024-02-16T20:41:06   00:00:00  COMPLETED      0.00G 31142559.ex+          gpu007
2024-02-16T22:17:24   00:00:00     FAILED      0.00G 31142559.0            gpu007
2024-02-16T22:22:31   00:00:04     FAILED      0.32G 31142559.1            gpu007
2024-02-16T22:24:29   00:00:04 CANCELLED+      0.33G 31142559.2            gpu007
2024-02-16T22:25:28   00:00:04     FAILED      0.33G 31142559.3            gpu007
2024-02-16T22:28:42   00:00:00 CANCELLED+      0.00G 31142559.4            gpu007
2024-02-16T22:28:53   03:27:07 CANCELLED+     31.38G 31142559.5            gpu007
```
So indeed you had a couple of jobs killed due to out of memory, I don’t know if this is the one you are talking about in your post?
Unless there is a bug, you can’t use the GPU from someone else. Once in your interactive session, just type the command `nvidia-smi`, you should ses only the GPUs allocated to you.
```
(yggdrasil)-[gabays@gpu001 ~]$ nvidia-smi
Tue Feb 27 08:08:11 2024
+---------------------------------------------------------------------------------------+
| NVIDIA-SMI 545.23.08              Driver Version: 545.23.08    CUDA Version: 12.3     |
|-----------------------------------------+----------------------+----------------------+
| GPU  Name                 Persistence-M | Bus-Id        Disp.A | Volatile Uncorr. ECC |
| Fan  Temp   Perf          Pwr:Usage/Cap |         Memory-Usage | GPU-Util  Compute M. |
|                                         |                      |               MIG M. |
|=========================================+======================+======================|
|   0  NVIDIA TITAN RTX               On  | 00000000:1A:00.0 Off |                  N/A |
| 41%   30C    P8               2W / 280W |      1MiB / 24576MiB |      0%      Default |
|                                         |                      |                  N/A |
+-----------------------------------------+----------------------+----------------------+
|   1  NVIDIA TITAN RTX               On  | 00000000:1B:00.0 Off |                  N/A |
| 41%   33C    P8              15W / 280W |      1MiB / 24576MiB |      0%      Default |
|                                         |                      |                  N/A |
+-----------------------------------------+----------------------+----------------------+

+---------------------------------------------------------------------------------------+
| Processes:                                                                            |
|  GPU   GI   CI        PID   Type   Process name                            GPU Memory |
|        ID   ID                                                             Usage      |
|=======================================================================================|
|  No running processes found                                                           |
+---------------------------------------------------------------------------------------+
(yggdrasil)-[gabays@gpu001 ~]$ echo $CUDA_VISIBLE_DEVICES
0,1
```
By the way in your example, you are requesting twice the GPUs (probably harmless, but I’m not sure which pragma is in fact used):
```
#SBATCH --gpus=1
#SBATCH --gres=gpu:1,VramPerGpu:24GB
```
If you want to specify the VramPerGpu, remove the first line. Please try again and check the memory usage. You can connect a second time to the compute node and launch `htop` while your job is running.
