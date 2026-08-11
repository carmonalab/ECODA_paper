# Compiling and running PhantomHFE on Baobab

- Source: https://hpc-community.unige.ch/t/3966

- Created: 2025-06-02T14:29:39.685Z

- Posts: 6

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yann.Sagon (2025-06-02T14:29:39.771Z)

A user asked us to install PhantomHFE. As this software isn’t available in EasyBuild below is a step-by-step guide to build the library and run it using SLURM as user.
This setup is based on the official PhantomHFE documentation[official PhantomHFE documentation](https://encryptorion-lab.gitbook.io/phantom-fhe/usages/python-applications).
:package: Step 1: Clone the Repository
```
git clone --recursive https://github.com/Encryptorion-Lab/phantom-fhe.git
cd phantom-fhe
```
:gear: Step 2: Load Required Modules
```
module load foss/2023a
module load CUDA/12.3.0
module load CMake/3.26.3
module load Python/3.11.3
```
:brick: Step 3: Configure and Build
```
mkdir build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_PYTHON_BINDINGS=ON
make -j 4
```
The compiled Python module pyPhantom will be located in build/lib.
:test_tube: Step 4: Use the Python Module
To run a Python script with the module, export the path:
```
export PYTHONPATH=$PWD/lib:$PYTHONPATH
```
Test the module:
```
python3 -c "import pyPhantom; print('Phantom module loaded successfully')"
```
:receipt: Step 5: Run with SLURM (sbatch)
Create a SLURM job script, e.g. phantom_job.sbatch:
```
#!/bin/bash
#SBATCH --job-name=phantom_ckks
#SBATCH --output=phantom_ckks.out
#SBATCH --error=phantom_ckks.err
#SBATCH --time=00:10:00
#SBATCH --partition=shared-gpu
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4

module load foss/2023a
module load CUDA/12.3.0
module load CMake/3.26.3
module load Python/3.11.3

cd /path/to/phantom-fhe/
export PYTHONPATH=$PYTHONPATH:build/lib

srun python3 python/examples/ckks.py
```
Replace /path/to/phantom-fhe with the actual path to your PhantomHFE repository.
:white_check_mark: Example Used: python/examples/ckks.py


## Post 2 by @walter.jauch (2025-06-04T07:13:30.824Z)

Hi !
Thank you for creating this topic.
I have tried to follow these exact steps but I keep getting the following error when trying to import the module in Python :
```
ModuleNotFoundError: No module named 'pyPhantom'
```
I tried to import the libPhantom.so in another way but it only changes the problem. I tried the following :
```
>>> import sys
>>> sys.path.insert(0, "lib")
>>> import libPhantom
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
ImportError: dynamic module does not define module export function (PyInit_libPhantom)
```
This was done in the PhantomFHE build folder.
I followed the steps described in the post above. Am I doing something wrong ?


## Post 3 by @Yann.Sagon (2025-06-06T07:02:31.285Z)

First thing to try is to remove the conda initialization from your `.bashrc` in your home directory. Then logout and login and try again.


## Post 4 by @walter.jauch (2025-06-07T12:38:04.651Z)

Hi,
I just tried removing the conda initialization from the .baschrc  and login in again but it seems that didn’t change the outcome. I still get the two errors I mentioned in my last reply…


## Post 5 by @Yann.Sagon (2025-06-10T11:57:21.828Z)

Hi @walter.jauch[@walter.jauch](https://hpc-community.unige.ch/u/walter.jauch)
it seems my instruction had a typo, probably a bad cut/paste sorry about that. Here is the correct line:
```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DPHANTOM_ENABLE_PYTHON_BINDING=ON
```


## Post 6 by @walter.jauch (2025-06-10T19:31:55.553Z)

Hi !
That worked !! I can now run the examples with no errors ! Thanks a lot !!
