# Install sklearn module on Python

- Source: https://hpc-community.unige.ch/t/3374

- Created: 2024-03-14T07:11:29.393Z

- Posts: 15

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Helena.Bach (2024-03-14T07:11:29.453Z)

We are running a stata code on baobab that uses a the command “ddml” which calls a python package :
// now check that sklearn is available
cap python which sklearn
if _rc > 0 {
local dostata = 1
`qui’ di as res “Python module sklearn not available; using Stata -nl-…”
}
ad we get this error in the log file in baobab:
File “”, line 1, in
ModuleNotFoundError: No module named ‘sklearn’
(176 lines skipped)
(error occurred while loading _ddml_nnls.ado)
r(7102);
Can you provide any help?
Thank you


## Post 2 by @Yann.Sagon (2024-03-19T08:10:53.791Z)

Hi,
please try to load the scikit-learn module before:
```
(baobab)-[sagon@login2 ~]$ ml GCC/12.3.0 scikit-learn/1.3.1 Stata/18
```
Best
Yann


## Post 3 by @Helena.Bach (2024-08-08T14:14:21.533Z)

Dear Yann, I am encountering the same problem with a different code. I tried to implement your suggestion but it gives me this error: Please check the spelling or version number. Also try “module spider …”
It is also possible your cache file is out-of-date; it may help to try:
$ module --ignore_cache load “scikit-learn/1.3.1/Stata/18”
## The stata code calls some python packages. Through the stata terminal i could check what packages are missing and I am missing the package sklearn and scikit-learn.
## python query
```
Python Settings
  set python_exec      /usr/bin/python3
  set python_userpath  

Python system information
  initialized          no
  version              3.6.8
  architecture         64-bit
  library path         /usr/lib64/libpython3.so
```
. python which numpy
<module ‘numpy’ from ‘/usr/lib64/python3.6/site-packages/numpy/init.py’>
. python which sklearn
Python module sklearn not found
r(601);


## Post 4 by @Adrien.Albert (2024-08-08T16:21:34.528Z)

Dear @Helena.Bach[@Helena.Bach](https://hpc-community.unige.ch/u/helena.bach)
There is an error in your load command:
```
(baobab)-[alberta@login1 ~]$ module load scikit-learn/1.3.1/Stata/18
Lmod has detected the following error:  The following module(s) are unknown: "scikit-learn/1.3.1/Stata/18"

Please check the spelling or version number. Also try "module spider ..."
It is also possible your cache file is out-of-date; it may help to try:
  $ module --ignore_cache load "scikit-learn/1.3.1/Stata/18"

Also make sure that all modulefiles written in TCL start with the string #%Module
```
It says the module is unknown (not installed or input error on name):
> Lmod has detected the following error:  The following module(s) are unknown: “scikit-learn/1.3.1/Stata/18”
```
-module scikit-learn/1.3.1/Stata/18
+module scikit-learn/1.3.1 Stata/18
                          ^ here the missing space
```
After correcting it Then you will have another warning:
```
(baobab)-[alberta@login1 ~]$ module load scikit-learn/1.3.1 Stata/18
Lmod has detected the following error:  These module(s) or extension(s) exist but cannot be loaded as requested: "scikit-learn/1.3.1"
   Try: "module spider scikit-learn/1.3.1" to see how to load the module(s).
```
Module output says:
> Try: “module spider scikit-learn/1.3.1” to see how to load the module(s).
So I do it:
```
(baobab)-[alberta@login1 ~]$ ml spider scikit-learn/1.3.1

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  scikit-learn: scikit-learn/1.3.1
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Description:
      Scikit-learn integrates machine learning algorithms in the tightly-knit scientific Python world, building upon numpy, scipy, and matplotlib. As a machine-learning module, it provides versatile tools
      for data mining and analysis in any field of science and engineering. It strives to be simple and efficient, accessible to everybody, and reusable in various contexts.

    You will need to load all module(s) on any one of the lines below before the "scikit-learn/1.3.1" module is available to load.

      GCC/12.3.0
 
    Help:
      Description
      ===========
      Scikit-learn integrates machine learning algorithms in the tightly-knit scientific Python world,
      building upon numpy, scipy, and matplotlib. As a machine-learning module,
      it provides versatile tools for data mining and analysis in any field of science and engineering.
      It strives to be simple and efficient, accessible to everybody, and reusable in various contexts.
      
      
      More information
      ================
       - Homepage: https://scikit-learn.org/stable/index.html
      
      
      Included extensions
      ===================
      scikit-learn-1.3.1
```
Module output says:
> You will need to load all module(s) on any one of the lines below before the “scikit-learn/1.3.1” module is available to load.
> **      GCC/12.3.0**
And Finally with all the tips merged I  run:
```
(baobab)-[alberta@login1 ~]$ ml GCC/12.3.0 scikit-learn/1.3.1 Stata/18
```
And it works :slight_smile:


## Post 5 by @Helena.Bach (2024-08-08T21:43:47.203Z)

## Thank you! It seems I can run the command ml GCC/12.3.0 scikit-learn/1.3.1 Stata/18 and it works. However, when checking if this package exists in stata it still doesn’t show up:
## . python query
```
Python Settings
  set python_exec      /usr/bin/python3
  set python_userpath  

Python system information
  initialized          no
  version              3.6.8
  architecture         64-bit
  library path         /usr/lib64/libpython3.so
```
python which sklearn
Python module sklearn not found
r(601);
end of do-file
r(601);


## Post 6 by @Adrien.Albert (2024-08-09T08:37:29.434Z)

Hi @Helena.Bach[@Helena.Bach](https://hpc-community.unige.ch/u/helena.bach)
Could you try formatting your message using executed data, quotation marks, etc.? It’s quite difficult to differentiate between blocks of code/output, file content and simple text.
Thank you


## Post 7 by @Helena.Bach (2024-08-09T08:42:24.930Z)

Dear Adrien,
Yes sorry for this!
I managed to get this code running on the terminal: `ml GCC/12.3.0 scikit-learn/1.3.1 Stata/18`. However, when I run the stata code, it still doesn’t recognize the command:  `Cross-fitting fold 1 unrecognized command r(199);` .
If I check whether Stata recognizes the sklearn package it doesn’t seem to be installed:
```
 `Python Settings
  set python_exec      /usr/bin/python3
  set python_userpath  

Python system information
  initialized          no
  version              3.6.8
  architecture         64-bit
  library path         /usr/lib64/libpython3.so
python which sklearn
Python module sklearn not found
r(601);
end of do-file
r(601);`
```


## Post 8 by @Adrien.Albert (2024-08-09T08:46:29.838Z)

Could you give the do file  and the sbatch file to reproduce the error ?


## Post 9 by @Helena.Bach (2024-08-09T08:48:42.595Z)

Do file: ddml_baobab.do
Sbatch file: a.sh


## Post 10 by @Adrien.Albert (2024-08-09T08:49:01.621Z)

For now, I can import the sklearn module on python:
```
(baobab)-[alberta@login1 ~]$ ml GCC/12.3.0 scikit-learn/1.3.1 Stata/18

(baobab)-[alberta@login1 ~]$ python
Python 3.11.3 (main, Jun 24 2024, 15:34:31) [GCC 12.3.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import sklearn
>>>
```


## Post 11 by @Adrien.Albert (2024-08-09T09:18:27.688Z)

@Helena.Bach[@Helena.Bach](https://hpc-community.unige.ch/u/helena.bach)
Make a little effort, I’m not in your head, nor are the other readers.
Here the code block interest us:
Stata code:
```
#/home/users/b/bachh/ddml_baobab.do
ssc install ddml, replace
ssc install pystacked, replace

*python query
*python which numpy
*python which sklearn
[...]
```
sbatch:
```
#/home/users/b/bachh/a.sh
#!/bin/sh

#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --mem=120G
#SBATCH -J job
#SBATCH -e job-error.e%j
#SBATCH -o job-out.o%j
#SBATCH --mail-user=xxxx
#SBATCH --mail-type=ALL
#SBATCH --partition=shared-bigmem

module load Stata/18

srun stata-mp ddml_baobab.do
```
---
## My answer:
From what I can see, your sbatch script doesn’t load the required modules that are mentioned in your code. Here’s an example of my setup that works:
Stata Code:
```
(baobab)-[alberta@login1 stata]$ cat bachh.do 
ssc install ddml, replace
ssc install pystacked, replace

*python query
*python which numpy
*python which sklearn
```
Sbatch script:
```
(baobab)-[alberta@login1 stata]$ cat bachh.sh 
#!/bin/sh

#SBATCH --time=00:15:00
#SBATCH --ntasks=1
#SBATCH --mem=120G
#SBATCH -J job
#SBATCH -e job-error.e%j
#SBATCH -o job-out.o%j
#SBATCH --partition=shared-bigmem

ml GCC/12.3.0 scikit-learn/1.3.1 Stata/18

srun stata-mp bachh.do
```
Here’s the output from my job (everything is working correctly):
```
(baobab)-[alberta@login1 stata]$ cat job-out.o11895037

  ___  ____  ____  ____  ____ ®
 /__    /   ____/   /   ____/      18.0
___/   /   /___/   /   /___/       MP—Parallel Edition

 Statistics and Data Science       Copyright 1985-2023 StataCorp LLC
                                   StataCorp
                                   4905 Lakeway Drive
                                   College Station, Texas 77845 USA
                                   800-STATA-PC        https://www.stata.com
                                   979-696-4600        stata@stata.com

Stata license: 2-user 32-core network perpetual
Serial number: 501806302009
  Licensed to: University of Geneva
               Geneva

Notes:
      1. Stata is running in batch mode.
      2. Unicode is supported; see help unicode_advice.
      3. More than 2 billion observations are allowed; see help obs_advice.
      4. Maximum number of variables is set to 5,000 but can be increased;
          see help set_maxvar.

. do "bachh.do" 

. ssc install ddml, replace
checking ddml consistency and verifying not already installed...
installing into /home/users/a/alberta/ado/plus/...
installation complete.

. ssc install pystacked, replace
checking pystacked consistency and verifying not already installed...
installing into /home/users/a/alberta/ado/plus/...
installation complete.

. 
. *python query
. *python which numpy
. *python which sklearn
. 
end of do-file
```
You need to make sure that you load all the necessary modules into your sbatch script. This is why your job doesn’t run as expected.
Let me know how it goes with this correction.


## Post 12 by @Helena.Bach (2024-08-09T09:36:27.533Z)

Thank you for your detailed answer. It still doesn’t work.
The code code that fails is the following:
```
ddml init partial, kfolds(2)
			ddml E[Y|X]: reg `var'_`y'  $controls 
			ddml E[Y|X]: pystacked `var'_`y'  $controls , type(reg) method(ols lassocv)
			ddml E[D|X]: reg treatment $controls
			ddml E[D|X]: pystacked treatment $controls , type(reg) method(ols lassocv)
			ddml desc
			ddml crossfit 
			ddml estimate, robust
```
This code uses “pystacked” which works via  scikit-learn[scikit-learn](https://scikit-learn.org/stable/). Since the package is not installed, it doesn’t recognize the command.
After installing the scikit-learn package in the sbatch file, as you mentioned, I check whether Stata is able to recognise this package. I do this at the beginning of the dofile using the following code:
```
. do "ddml_baobab.do" 

. ssc install ddml, replace
checking ddml consistency and verifying not already installed...
all files already exist and are up to date.

. ssc install pystacked, replace
checking pystacked consistency and verifying not already installed...
all files already exist and are up to date.

. 
. python query
-------------------------------------------------------------------------------
    Python Settings
      set python_exec      /usr/bin/python3
      set python_userpath  

    Python system information
      initialized          no
      version              3.6.8
      architecture         64-bit
      library path         /usr/lib64/libpython3.so

. python
----------------------------------------------- python (type end to exit) -----
>>> end
-------------------------------------------------------------------------------

. python query
-------------------------------------------------------------------------------
    Python Settings
      set python_exec      /usr/bin/python3
      set python_userpath  

    Python system information
      initialized          yes
      version              3.6.8
      architecture         64-bit
      library path         /usr/lib64/libpython3.so

. *python which numpy
. python which sklearn
Python module sklearn not found
r(601);
```
I first initialize python, and then check whether it can identify sklearn. And it says “module sklearn not found”.
In the internet I’ve found the following guidelines to solve the " `crossfitting fold 1 unrecognized command`"  problem:
1、“ssc install ddml, replace” in Stata (which is already done in the dofile)
2、“ssc install pystacked, replace” in Stata (which is already done in the dofile)
3、“pip install scikit-learn” in cmd, to configure the python environment (
4、Try to use do-file in Stata to run the codes


## Post 13 by @Adrien.Albert (2024-08-09T09:54:10.776Z)

Helena,
I am sorry but I do not understand what is going wrong.
I can reproduce the error:
Output:
```
[..]
. python which sklearn
Python module sklearn not found
r(601);
```
But I correct it by loading (not installing) the module in my sbatch.
I just sent you an email for quick zoom meeting.


## Post 14 by @Adrien.Albert (2024-08-09T11:31:40.409Z)

After checking,
The problem comes from stata which uses the Python system instead of the Python loaded with module . Despite forcing the Python version, there is no change.
Output:
```
.  set python_exec /opt/ebsofts/Python/3.11.3-GCCcore-12.3.0/bin/python

.  python search
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 Python environments found:  
 /usr/bin/python2
 /usr/bin/python3
 /usr/bin/python3.6
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

.  python query
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Python Settings
      set python_exec      /usr/bin/python3
      set python_userpath  

    Python system information
      initialized          no
      version              3.6.8
      architecture         64-bit
      library path         /usr/lib64/libpython3.so
```
BUT with black magic :crystal_ball: it works, and now python comes from easybuild:
Beetween the two test I tried on OpenOnDemand with graphical Desktop (xstata)seeking some debug mode, I loaded the module by hand and now it’s working.
image
image1916×847 171 KB
[image1916×847 171 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/f185c2aa2ba8fce05327f13e3dd8ec3eaffd2b4b.jpeg)
And now with sbatch:
```
(baobab)-[alberta@cpu001 stata]$ stata-mp bachh.do

  ___  ____  ____  ____  ____ ®
 /__    /   ____/   /   ____/      18.0
___/   /   /___/   /   /___/       MP—Parallel Edition

 Statistics and Data Science       Copyright 1985-2023 StataCorp LLC
                                   StataCorp
                                   4905 Lakeway Drive
                                   College Station, Texas 77845 USA
                                   800-STATA-PC        https://www.stata.com
                                   979-696-4600        stata@stata.com

Stata license: 2-user 32-core network perpetual
Serial number: 501806302009
  Licensed to: University of Geneva
               Geneva

Notes:
      1. Unicode is supported; see help unicode_advice.
      2. More than 2 billion observations are allowed; see help obs_advice.
      3. Maximum number of variables is set to 5,000 but can be increased; see help set_maxvar.

. do "bachh.do" 

. *ssc install ddml, replace
. * ssc install pystacked, replace
. 
.  set python_exec /opt/ebsofts/Python/3.11.3-GCCcore-12.3.0/bin/python

.  python query
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Python Settings
      set python_exec      /opt/ebsofts/Python/3.11.3-GCCcore-12.3.0/bin/python
      set python_userpath  

    Python system information
      initialized          no
      version              3.11.3
      architecture         64-bit
      library path         /opt/ebsofts/Python/3.11.3-GCCcore-12.3.0/lib/libpython3.11.so

.  python which numpy
<module 'numpy' from '/opt/ebsofts/SciPy-bundle/2023.07-gfbf-2023a/lib/python3.11/site-packages/numpy/__init__.py'>

.  python which sklearn
<module 'sklearn' from '/opt/ebsofts/scikit-learn/1.3.1-gfbf-2023a/lib/python3.11/site-packages/sklearn/__init__.py'>

. 
end of do-file
```
It works with all slurm allocation methods (salloc, sbatch, srun and OpenOnDemand). I have no explanation for this resolution…


## Post 15 by @Helena.Bach (2024-08-09T11:44:27.117Z)

Thank you so much!
It is working indeed
