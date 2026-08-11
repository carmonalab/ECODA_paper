# Libvips library yggdrasil

- Source: https://hpc-community.unige.ch/t/3589

- Created: 2024-08-09T08:40:01.122Z

- Tags: yggdrasil

- Posts: 6

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Maxime.Humeau (2024-08-09T08:40:01.165Z)

Hello @support[@support](https://hpc-community.unige.ch/groups/support) !
## Primary informations
Username: humeau
Cluster: Yggdrasil
## Description
I want to use a python script requiring the pyvips library. To use this library I need libvips which is the original library, pyvips is only a binding.
I tried to find the module with spider, but it was not successful. Is it possible to find this library via another module?
## Actual Result
OSError: cannot load library ‘libvips.so.42’: libvips.so.42: cannot open shared object file: No such file or directory.  Additionally, ctypes.util.find_library() did not manage to locate a library called ‘libvips.so.42’
Any help would be greatly appreciated, thank you !
Maxime


## Post 2 by @Adrien.Albert (2024-08-09T09:38:17.242Z)

Hi @Maxime.Humeau[@Maxime.Humeau](https://hpc-community.unige.ch/u/maxime.humeau)
To find out more:
Did you install pyvips in your home, with the Python version you use ?
Could you give me a small code to reproduce the issue (code + sbatch) ?


## Post 3 by @Maxime.Humeau (2024-08-09T10:17:43.399Z)

Thanks you for your support !
for sbatch
```
#!/bin/env bash

#SBATCH --partition=debug-gpu
#SBATCH --time=00:15:00
#SBATCH --output=journal-%j.out
#SBATCH --mem=8000
#SBATCH --gres=gpu:1,VramPerGpu:20GB
#SBATCH --cpus-per-task=10
module load cuDNN/8.6.0.163-CUDA-11.8.0 GCCcore/10.3.0 Python/3.9.5 OpenSSL/1.1.1q

echo "I: CUDA_VISIBLE_DEVICES: ${CUDA_VISIBLE_DEVICES}"

echo "====="

srun nvidia-smi

echo "====="

#srun python3  -m venv yaltaienv
source yaltaienv/bin/activate

#wget https://github.com/FoNDUE-HTR/Documentation/releases/download/v.0.9/fondue_emmental.mlmodel -O htr_model.mlmodel
wget https://zenodo.org/records/10972956/files/CapricciosaX.pt?download=1 -O seg_model.pt

pip install -r requirements.txt
pip install yaltai~=1.0.0 --extra-index-url https://download.pytorch.org/whl/cu118

srun python3 run.py
```
I’m installing pyvips with `requirements.txt`
For the code, i’m reusing this (GitHub - PonteIneptique/rtk at a9902aa1d21b4c57cd28fba71f81e1023ffeb14a[GitHub - PonteIneptique/rtk at a9902aa1d21b4c57cd28fba71f81e1023ffeb14a](https://github.com/PonteIneptique/rtk/tree/a9902aa1d21b4c57cd28fba71f81e1023ffeb14a)) in same environment. I don’t separate yaltai and kraken


## Post 4 by @Adrien.Albert (2024-08-09T12:09:41.659Z)

Maxime.Humeau:
> I’m installing pyvips with `requirements.txt`
You don’t pull the git repository to make this instruction works, it’s this one ? :
Maxime.Humeau:
> (GitHub - PonteIneptique/rtk at a9902aa1d21b4c57cd28fba71f81e1023ffeb14a [GitHub - PonteIneptique/rtk at a9902aa1d21b4c57cd28fba71f81e1023ffeb14a](https://github.com/PonteIneptique/rtk/tree/a9902aa1d21b4c57cd28fba71f81e1023ffeb14a))


## Post 5 by @Maxime.Humeau (2024-08-09T12:17:55.546Z)

Adrien.Albert:
> You don’t pull the git repository to make this instruction works, it’s this one ? :
I’m using an older version at the moment. But in any case, the new version still uses libvips via pipvips binding. I download in zip and i transfert with ssh.


## Post 6 by @Yann.Sagon (2024-08-15T06:05:53.078Z)

Hi,
libvips isn’t provided by easybuild (module). Did you tried to install libvips in you home by following the doc: libvips[libvips](https://www.libvips.org/install.html) => Building libvips from source ?
