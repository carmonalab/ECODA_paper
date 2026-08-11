# NETCDF library installation

- Source: https://hpc-community.unige.ch/t/3705

- Created: 2024-10-29T14:10:19.011Z

- Posts: 10

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Lyudmila.Utkina (2024-10-29T14:10:19.045Z)

## Primary informations
Username: utkina1
Cluster: Baobab
## Description
Can’t install the netcdf library in the created environment within Baobab. Gives the following error:
```
 pip install python=3.7/home/users/u/utkina1/baobab_python_env/bin/python3.7: error while loading shared libraries: libpython3.7m.so.1.0: cannot open shared object file: No such file or directory(baobab_python_env) 
```
Thanks for help
admin edit: we received this screen shot from you by email, adding it to the post
image
image1323×121 14.3 KB
[image1323×121 14.3 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/da73b69e89892d54725c9074576912ea93c0e119.png)


## Post 2 by @Yann.Sagon (2024-10-31T10:16:35.416Z)

Dear @Lyudmila.Utkina[@Lyudmila.Utkina](https://hpc-community.unige.ch/u/lyudmila.utkina)
we are already providing some variant of netCDF4-python:
```
(baobab)-[sagon@login1 ~]$ ml spider netcdf4-python

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  netcdf4-python:
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Description:
      Python/numpy interface to netCDF.

     Versions:
        netcdf4-python/1.4.2-Python-3.7.2
        netcdf4-python/1.5.7-Python-3.9.5
        netcdf4-python/1.5.7
```
Then you can load it:
```
(baobab)-[sagon@login1 ~]$ ml  GCC/10.3.0  OpenMPI/4.1.1 netcdf4-python/1.5.7
```
And use it:
```
(baobab)-[sagon@login1 ~]$ python
Python 3.9.5 (default, Aug 31 2021, 16:56:58)
[GCC 10.3.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> from netCDF4 import Dataset
```
Let me know if it is working for you.
Best


## Post 3 by @Lyudmila.Utkina (2024-10-31T10:46:24.878Z)

Yann.Sagon:
> `GCC/10.3.0  OpenMPI/4.1.1 netcdf4-python/1.5.7`
Hi And thanks for your reply
However I am facing the same problem - when I am sending my script it returns the error
image
image1108×296 30.3 KB
[image1108×296 30.3 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/b3df28ece59d3b567c87e5114c4ad8fd910c32fd.png)
What to do in this case?
Thank you again,
L.


## Post 4 by @Lyudmila.Utkina (2024-10-31T10:47:09.861Z)

(to precise that I have checked “module spider” for both packages and loaded everything (seemingly) )


## Post 5 by @Yann.Sagon (2024-10-31T10:53:46.318Z)

Lyudmila.Utkina:
> However I am facing the same problem - when I am sending my script it returns the error
Well, the error message has nothing in common with your old message, so you mean the same problem is because it isn’t working, right?
Please show your script and I’ll try to help.


## Post 6 by @Lyudmila.Utkina (2024-10-31T11:07:43.424Z)

Yes, because it isnt working.
Thank you.
Here is the script:
```
import netcdf4 as nc

import xarray as xr

import matplotlib.pyplot as plt

file_path = "/home/users/u/utkina1/scratch/mon_200101-200512.nc"

dataset = xr.open_dataset(file_path)

mean_precip = dataset['pr'].mean(dim=['rlat', 'rlon'])

median_precip = dataset['pr'].median(dim=['rlat', 'rlon'])

import matplotlib.pyplot as plt

plt.figure(figsize=(10, 6))

mean_precip.plot(label='Mean Precipitation', color='blue')

plt.xlabel('Time')

plt.ylabel('Precipitation (mm/day)')

plt.title('Mean Precipitation 1991-2000 monthly (Europe)')

plt.legend()

plt.grid(True)

plt.show()

plt.figure(figsize=(10, 6))

median_precip.plot(label='Median Precipitation', color='orange')

plt.xlabel('Time')

plt.ylabel('Precipitation (mm/day)')

plt.title('Median Precipitation 1991-2000 monthly (Europe)')

plt.legend()

plt.grid(True)

plt.show()

plt.figure(figsize=(10, 6))

mean_precip.plot(label='Mean Precipitation', color='blue')

median_precip.plot(label='Median Precipitation', color='orange')

plt.xlabel('Time')

plt.ylabel('Precipitation (mm/day)')

plt.title('Mean vs Median Precipitation 1991-2000 (Europe)')

plt.legend()

plt.grid(True)

from scipy.stats import wilcoxon

mean_values = mean_precip.values

median_values = median_precip.values

stat, p_value = wilcoxon(mean_values, median_values)

print(f'Wilcoxon signed-rank test statistic: {stat}')

print(f'P-value: {p_value}')
```


## Post 7 by @Lyudmila.Utkina (2024-10-31T11:09:24.026Z)

And here is the bash script:
```
#!/bin/bash
#SBATCH --job-name=python      # Job name
#SBATCH --output=python_output_%j.log        # Output file (log of your job)
#SBATCH --error=python_error_%j.log          # Error file (log of errors)
#SBATCH --partition=private-gap-cpu         # Partition (change if needed)
#SBATCH --nodes=1                     # Number of nodes
#SBATCH --ntasks=1                    # Number of tasks (1 task for serial jobs)
#SBATCH --cpus-per-task=4             # Number of CPU cores per task
#SBATCH --time=00:10:00               # Maximum run time (hh:mm:ss)
#SBATCH --mem=64G                      # Memory per node
#SBATCH --mail-type=END,FAIL          # Send mail on job end/fail
#SBATCH --mail-user=lyudmila.utkina@etu.unige.ch # Email address to send notifications

# Load any necessary modules (adjust to your cluster's environment)

# Activate your Python virtual environment if required
source /home/users/u/utkina1/baobab_python_env/bin/activate  # Adjust with your virtualenv path

module purge  # Clear any previously loaded modules
module load GCCcore/13.2.0 Python/3.11.5  # Example of loading a Python module
ml load GCC/11.3.0  OpenMPI/4.1.4

module load pandas
module load numpy
module load xarray
module load netcdf4-python/1.5.7 

# Run your Python script
python3 /home/users/u/utkina1/scratch/test.py
```


## Post 8 by @Yann.Sagon (2024-10-31T12:43:35.742Z)

Hi,
please see comments below about your sbatch script.
Line6: you can remove it (`#SBATCH --nodes=1`) only needed if you want more than one compute node (for a good reason)
Line10: be sure to only request this option if default memory isn’t enough. If you don’t know, you can remove it.
Line18: if you use our version of netCDF, why do you still need your virtualenv? As you are activating your virtualenv without any Python module loaded, you’ll get the default Python from system which is quite old. You probably need to remove the line.
Line21: you are loading GCC 13.2
Line22: you are loading GCC 11.3 => this will unload some of the modules loades from your previous load!
Line27: you are loading netcdf4-python but without the required dependencies, thus it fails.
Rules of thumb: you need to load all the modules for a single GCC version, i.e. choose one and stick to it.
`Pandas` is now provided by the module `SciPy-bundle`. You can see what is provided by this module like that:
```
(baobab)-[sagon@login1 ~]$ ml spider SciPy-bundle/2023.07

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  SciPy-bundle: SciPy-bundle/2023.07
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Description:
      Bundle of Python packages for scientific software

    You will need to load all module(s) on any one of the lines below before the "SciPy-bundle/2023.07" module is available to load.

      GCC/12.3.0

    Help:
      Description
      ===========
      Bundle of Python packages for scientific software

      More information
      ================
       - Homepage: https://python.org/

      Included extensions
      ===================
      beniget-0.4.1, Bottleneck-1.3.7, deap-1.4.0, gast-0.5.4, mpmath-1.3.0,
      numexpr-2.8.4, numpy-1.25.1, pandas-2.0.3, ply-3.11, pythran-0.13.1,
      scipy-1.11.1, tzdata-2023.3, versioneer-0.29
```
Numpy is as well included in the same module. And SciPy-bundle is a dependency for netcdf4-python/1.5.7, so no need to load it, it is automatically loaded. You just need to add `xarray/0.19.0` and that’s all.
New sbatch:
```
#!/bin/bash
#SBATCH --job-name=python               # Job name
#SBATCH --output=python_output_%j.log   # Output file (log of your job)
#SBATCH --error=python_error_%j.log     # Error file (log of errors)
#SBATCH --partition=private-gap-cpu     # Partition (change if needed)
#SBATCH --ntasks=1                      # Number of tasks (1 task for serial jobs)
#SBATCH --cpus-per-task=4               # Number of CPU cores per task
#SBATCH --time=00:10:00                 # Maximum run time (hh:mm:ss)

# Load any necessary modules (adjust to your cluster's environment)

module purge  # Clear any previously loaded modules
# load Numpy, xarray, netcdf4-python, pandas
ml GCC/10.3.0 OpenMPI/4.1.1 netcdf4-python/1.5.7 xarray/0.19.0

# Run your Python script
srun python3 /home/users/u/utkina1/scratch/test.py
```


## Post 9 by @Lyudmila.Utkina (2024-12-10T10:22:51.551Z)

Hi.
I didn’t need hpc over the last month but now when I am back to it I still have this issue with netcdf4.
The error that is shown after running the code:
```
Traceback (most recent call last):
  File "/home/users/u/utkina1/scratch/Histogram_quantiles.py", line 1, in <module>
srun: error: cpu224: task 0: Exited with exit code 1
    import netcdf4
ModuleNotFoundError: No module named 'netcdf4'
```
My bash script:
```
#!/bin/bash
#SBATCH --job-name=python               # Job name
#SBATCH --output=python_output_%j.log   # Output file (log of your job)
#SBATCH --error=python_error_%j.log     # Error file (log of errors)
#SBATCH --partition=private-gap-cpu     # Partition (change if needed)
#SBATCH --ntasks=1                      # Number of tasks (1 task for serial jobs)
#SBATCH --cpus-per-task=4               # Number of CPU cores per task
#SBATCH --time=00:10:00                 # Maximum run time (hh:mm:ss)

# Load any necessary modules (adjust to your cluster's environment)

module purge  # Clear any previously loaded modules
# load Numpy, xarray, netcdf4-python, pandas
module load GCC/10.3.0 OpenMPI/4.1.1
module load netcdf4-python/1.5.7-Python-3.9.5 
module load xarray/0.19.0 matplotlib/3.4.2

# Run your Python script
srun python3 /home/users/u/utkina1/scratch/Histogram_quantiles.py
```
The start of my Python script:
```
import netcdf4

import matplotlib.pyplot as plt

import numpy as np

# Load the NetCDF file

dataset = netcdf4.Dataset("/home/users/u/utkina1/scratch/mon_200101-200512.nc")

# Extract variables

time = dataset.variables['time'][:]

pr = dataset.variables['pr'][:, :, :]
```
I have installed all the modules and can’t comprehend what could be going wrong.
Please help, thanks.
Lucy


## Post 10 by @Yann.Sagon (2024-12-10T16:07:32.374Z)

Hi,
the python modules are case-sensitive and the correct name (at least in our installation) is `netCDF4`
```
(baobab)-[sagon@login1 ~]$ python3
Python 3.9.5 (default, Aug 31 2021, 16:56:58)
[GCC 10.3.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import netCDF4
```
Best
