# Specific module in python unavailable [job 41146649]

- Source: https://hpc-community.unige.ch/t/3983

- Created: 2025-06-17T11:45:20.520Z

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Louis.Perrin (2025-06-17T11:45:20.543Z)

Dear HPC staff,
First of all, thank you very much for all your help throughout the year!
I’m contacting you today regarding a python script I’d like to run on an array of jobs on HPC because it’s very time consuming to just run it locally (a couple of weeks, which I managed to get down to a couple of days thanks to python multiprocessing, and I’d like to convert it to an array of job to make it faster).
My current issue is the list of modules in python I’m using (list of imports hereunder):
from fenics import *
import time
import numpy as np
import sys
from scipy.interpolate import RegularGridInterpolator
import time
import os
Ive looked in the documentation and the forum, but didn’t find any mention of the fenics package on the HPC, which is why my code doesn’t run at the moment.
Would it be possible to install that package somehow? I unfortunately don’t know any equivalent for what I’m trying to do…
Here is the link to the package in case it is of any use: https://fenics.readthedocs.io/en/latest/installation.html[https://fenics.readthedocs.io/en/latest/installation.html](https://fenics.readthedocs.io/en/latest/installation.html)
Best regards,
Louis Perrin


## Post 2 by @Louis.Perrin (2025-06-18T11:20:34.666Z)




## Post 3 by @Yann.Sagon (2025-06-19T12:18:49.616Z)




## Post 4 by @Yann.Sagon (2025-06-20T08:31:47.630Z)

Dear @Louis.Perrin[@Louis.Perrin](https://hpc-community.unige.ch/u/louis.perrin)
I’ve tried to use the Docker image without luck as it is very outdated (latest release in 2020).
I tried then to install it as user:
```
(baobab)-[sagon@login1 ~]$ pip install fenics --user
Collecting fenics
  Downloading fenics-2019.1.0-py3-none-any.whl (1.5 kB)
Collecting fenics-ffc<2019.2,>=2019.1.0
  Downloading fenics_ffc-2019.1.0.post0-py3-none-any.whl (362 kB)
     |████████████████████████████████| 362 kB 4.4 MB/s
Collecting fenics-dijitso<2019.2,>=2019.1.0
  Downloading fenics_dijitso-2019.1.0-py3-none-any.whl (46 kB)
     |████████████████████████████████| 46 kB 1.8 MB/s
Collecting fenics-ufl<2019.2,>=2019.1.0
  Downloading fenics_ufl-2019.1.0-py3-none-any.whl (282 kB)
     |████████████████████████████████| 282 kB 121.9 MB/s
Collecting fenics-fiat<2019.2,>=2019.1.0
  Downloading fenics_fiat-2019.1.0-py3-none-any.whl (112 kB)
     |████████████████████████████████| 112 kB 133.2 MB/s
Collecting numpy
  Downloading numpy-2.0.2-cp39-cp39-manylinux_2_17_x86_64.manylinux2014_x86_64.whl (19.5 MB)
     |████████████████████████████████| 19.5 MB 56.5 MB/s
Collecting sympy
  Downloading sympy-1.14.0-py3-none-any.whl (6.3 MB)
     |████████████████████████████████| 6.3 MB 126.4 MB/s
Collecting mpmath<1.4,>=1.1.0
  Downloading mpmath-1.3.0-py3-none-any.whl (536 kB)
     |████████████████████████████████| 536 kB 136.2 MB/s
Installing collected packages: mpmath, sympy, numpy, fenics-ufl, fenics-fiat, fenics-dijitso, fenics-ffc, fenics
  WARNING: The script isympy is installed in '/home/sagon/.local/bin' which is not on PATH.
  Consider adding this directory to PATH or, if you prefer to suppress this warning, use --no-warn-script-location.
  WARNING: The scripts f2py and numpy-config are installed in '/home/sagon/.local/bin' which is not on PATH.
  Consider adding this directory to PATH or, if you prefer to suppress this warning, use --no-warn-script-location.
  WARNING: The script dijitso is installed in '/home/sagon/.local/bin' which is not on PATH.
  Consider adding this directory to PATH or, if you prefer to suppress this warning, use --no-warn-script-location.
  WARNING: The scripts ffc and ffc-3 are installed in '/home/sagon/.local/bin' which is not on PATH.
  Consider adding this directory to PATH or, if you prefer to suppress this warning, use --no-warn-script-location.
Successfully installed fenics-2019.1.0 fenics-dijitso-2019.1.0 fenics-ffc-2019.1.0.post0 fenics-fiat-2019.1.0 fenics-ufl-2019.1.0 mpmath-1.3.0 numpy-2.0.2 sympy-1.14.0
```
Test:
```
(baobab)-[sagon@login1 ~]$ python
Python 3.9.21 (main, Dec  5 2024, 00:00:00)
[GCC 11.5.0 20240719 (Red Hat 11.5.0-2)] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import ffc
>>>
```
Is this working for you like that?
