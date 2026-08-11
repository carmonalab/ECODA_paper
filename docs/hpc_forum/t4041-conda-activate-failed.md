# Conda activate failed

- Source: https://hpc-community.unige.ch/t/4041

- Created: 2025-08-11T20:23:43.753Z

- Posts: 13

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Camille.Christe (2025-08-11T20:23:43.795Z)

Dear HPC team,
I moved my data to bamboo and I tried to run the program Hybpiper on it.
For this I loaded Anaconda to install Hybpiper via conda
```
 ml Anaconda3/2024.02-1
conda create -c bioconda -c cona-forge -n hybpiper hybpiper
```
It seems to work properly, I could see the environment via conda env list
But when I tried to activate the environment, I got this message:
```
CondaError: Run 'conda init' before 'conda activate'
```
Then I am not able to write on the terminal.
Is there something that I missed to activate the environment I loaded?
Best regards,
Camille


## Post 2 by @Gael.Rossignol (2025-08-12T09:08:46.981Z)

Dear Camille,
Could you please test to launch on login1 command “conda init”?
Best regards,


## Post 3 by @Camille.Christe (2025-08-14T05:21:38.495Z)

Yes, here is what I get:
```
conda init

We trust you have received the usual lecture from the local System
Administrator. It usually boils down to these three things:

    #1) Respect the privacy of others.
    #2) Think before you type.
    #3) With great power comes great responsibility.

[sudo] password for christec: 
# >>>>>>>>>>>>>>>>>>>>>> ERROR REPORT <<<<<<<<<<<<<<<<<<<<<<

    Traceback (most recent call last):
      File "/opt/ebsofts/Anaconda3/2024.02-1/lib/python3.11/site-packages/conda/exception_handler.py", line 17, in __call__
        return func(*args, **kwargs)
               ^^^^^^^^^^^^^^^^^^^^^
      File "/opt/ebsofts/Anaconda3/2024.02-1/lib/python3.11/site-packages/conda/cli/main.py", line 83, in main_subshell
        exit_code = do_call(args, parser)
                    ^^^^^^^^^^^^^^^^^^^^^
      File "/opt/ebsofts/Anaconda3/2024.02-1/lib/python3.11/site-packages/conda/cli/conda_argparse.py", line 196, in do_call
        result = getattr(module, func_name)(args, parser)
                 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      File "/opt/ebsofts/Anaconda3/2024.02-1/lib/python3.11/site-packages/conda/cli/main_init.py", line 155, in execute
        return initialize(
               ^^^^^^^^^^^
      File "/opt/ebsofts/Anaconda3/2024.02-1/lib/python3.11/site-packages/conda/core/initialize.py", line 139, in initialize
        run_plan_elevated(plan2)
      File "/opt/ebsofts/Anaconda3/2024.02-1/lib/python3.11/site-packages/conda/core/initialize.py", line 856, in run_plan_elevated
        result = subprocess_call(
                 ^^^^^^^^^^^^^^^^
      File "/opt/ebsofts/Anaconda3/2024.02-1/lib/python3.11/site-packages/conda/gateways/subprocess.py", line 102, in subprocess_call
        stdout, stderr = process.communicate(input=stdin)
                         ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      File "/opt/ebsofts/Anaconda3/2024.02-1/lib/python3.11/subprocess.py", line 1209, in communicate
        stdout, stderr = self._communicate(input, endtime, timeout)
                         ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      File "/opt/ebsofts/Anaconda3/2024.02-1/lib/python3.11/subprocess.py", line 2088, in _communicate
        input_view = memoryview(self._input)
                     ^^^^^^^^^^^^^^^^^^^^^^^
    TypeError: memoryview: a bytes-like object is required, not 'str'

`$ /opt/ebsofts/Anaconda3/2024.02-1/bin/conda init`

  environment variables:
                 CIO_TEST=<not set>
               CONDA_ROOT=/opt/ebsofts/Anaconda3/2024.02-1
           CURL_CA_BUNDLE=<not set>
                    FPATH=/usr/share/lmod/lmod/init/ksh_funcs
               LD_PRELOAD=<not set>
                  MANPATH=/opt/ebsofts/Anaconda3/2024.02-
                          1/share/man:/opt/ebsofts/Anaconda3/2024.02-
                          1/man:/usr/share/lmod/lmod/share/man::
               MODULEPATH=/etc/modulefiles:/usr/share/modulefiles:/opt/modulefiles/Linux:/usr/sh
                          are/lmod/lmod/modulefiles/Core:/opt/ebmodules/all/Core
                     PATH=/opt/ebsofts/Anaconda3/2024.02-1:/opt/ebsofts/Anaconda3/2024.02-
                          1/sbin:/opt/ebsofts/Anaconda3/2024.02-
                          1/bin:/home/christec/.local/bin:/home/christec/bin:/usr/local/bin:/usr
                          /bin:/usr/local/sbin:/usr/sbin
          PKG_CONFIG_PATH=/opt/ebsofts/Anaconda3/2024.02-1/lib/pkgconfig
       REQUESTS_CA_BUNDLE=<not set>
            SSL_CERT_FILE=<not set>
 __LMOD_REF_COUNT_MANPATH=/opt/ebsofts/Anaconda3/2024.02-
                          1/share/man:1;/opt/ebsofts/Anaconda3/2024.02-
                          1/man:1;/usr/share/lmod/lmod/share/man:1;:1
__LMOD_REF_COUNT_MODULEPATH=/etc/modulefiles:1;/usr/share/modulefiles:1;/opt/modulefiles/Linux:1;/
                          usr/share/lmod/lmod/modulefiles/Core:1;/opt/ebmodules/all/Core:1
    __LMOD_REF_COUNT_PATH=/opt/ebsofts/Anaconda3/2024.02-1:1;/opt/ebsofts/Anaconda3/2024.02-
                          1/sbin:1;/opt/ebsofts/Anaconda3/2024.02-
                          1/bin:1;/home/christec/.local/bin:1;/home/christec/bin:1;/usr/local/bi
                          n:1;/usr/bin:1;/usr/local/sbin:1;/usr/sbin:1
__LMOD_REF_COUNT_PKG_CONFIG_PATH=/opt/ebsofts/Anaconda3/2024.02-1/lib/pkgconfig:1

     active environment : None
       user config file : /home/christec/.condarc
 populated config files : 
          conda version : 24.1.2
    conda-build version : 24.1.2
         python version : 3.11.7.final.0
                 solver : libmamba (default)
       virtual packages : __archspec=1=zen2
                          __conda=24.1.2=0
                          __glibc=2.34=0
                          __linux=5.14.0=0
                          __unix=0=0
       base environment : /opt/ebsofts/Anaconda3/2024.02-1  (read only)
      conda av data dir : /opt/ebsofts/Anaconda3/2024.02-1/etc/conda
  conda av metadata url : None
           channel URLs : https://repo.anaconda.com/pkgs/main/linux-64
                          https://repo.anaconda.com/pkgs/main/noarch
                          https://repo.anaconda.com/pkgs/r/linux-64
                          https://repo.anaconda.com/pkgs/r/noarch
          package cache : /opt/ebsofts/Anaconda3/2024.02-1/pkgs
                          /home/christec/.conda/pkgs
       envs directories : /home/christec/.conda/envs
                          /opt/ebsofts/Anaconda3/2024.02-1/envs
               platform : linux-64
             user-agent : conda/24.1.2 requests/2.31.0 CPython/3.11.7 Linux/5.14.0-503.40.1.el9_5.x86_64 rocky/9.5 glibc/2.34 solver/libmamba conda-libmamba-solver/24.1.0 libmambapy/1.5.6
                UID:GID : 209112:1000
             netrc file : None
           offline mode : False

An unexpected error has occurred. Conda has prepared the above report.
If you suspect this error is being caused by a malfunctioning plugin,
consider using the --no-plugins option to turn off plugins.

Example: conda --no-plugins install <package>

Alternatively, you can set the CONDA_NO_PLUGINS environment variable on
the command line to run the command without plugins enabled.

Example: CONDA_NO_PLUGINS=true conda install <package>

If submitted, this report will be used by core maintainers to improve
future releases of conda.
Would you like conda to send this report to the core maintainers? [y/N]: 
Sorry, try again.
[sudo] password for christec: 
No report sent. To permanently opt-out, use

    $ conda config --set report_errors false
```


## Post 4 by @Yann.Sagon (2025-08-15T08:32:23.953Z)

Hello @Camille.Christe[@Camille.Christe](https://hpc-community.unige.ch/u/camille.christe) it seems you used `sudo`? this isn’t permitted on the cluster. I’ve tried as user and it works for me:
Create the environment
```
(baobab)-[sagon@admin1 ~]$ conda create -c bioconda -c conda-forge -n hybpiper hybpiper
```
Try to activate (failed because activate wasn’t run previously
```
(baobab)-[sagon@admin1 ~]$ conda activate hybpiper

CondaError: Run 'conda init' before 'conda activate'
```
```
(baobab)-[sagon@admin1 ~]$ conda init
no change     /opt/ebsofts/Anaconda3/2024.02-1/condabin/conda
no change     /opt/ebsofts/Anaconda3/2024.02-1/bin/conda
no change     /opt/ebsofts/Anaconda3/2024.02-1/bin/conda-env
no change     /opt/ebsofts/Anaconda3/2024.02-1/bin/activate
no change     /opt/ebsofts/Anaconda3/2024.02-1/bin/deactivate
no change     /opt/ebsofts/Anaconda3/2024.02-1/etc/profile.d/conda.sh
no change     /opt/ebsofts/Anaconda3/2024.02-1/etc/fish/conf.d/conda.fish
no change     /opt/ebsofts/Anaconda3/2024.02-1/shell/condabin/Conda.psm1
no change     /opt/ebsofts/Anaconda3/2024.02-1/shell/condabin/conda-hook.ps1
no change     /opt/ebsofts/Anaconda3/2024.02-1/lib/python3.11/site-packages/xontrib/conda.xsh
no change     /opt/ebsofts/Anaconda3/2024.02-1/etc/profile.d/conda.csh
modified      /home/sagon/.bashrc

==> For changes to take effect, close and re-open your current shell. <==
```
I’ve checked, your `.bashrc` is already modified by conda.
By the way, our recommended solution when dealing with `conda` is to use cotainr[cotainr](https://doc.eresearch.unige.ch/hpc/applications_and_libraries#how_to_create_a_conda_environment_in_a_container) to encapsulate the environment in a container.. It is straightforward and you’ll get better performance.


## Post 5 by @Camille.Christe (2025-08-15T11:07:41.726Z)

I have tried the same way as you but after conda init I get the same message as before. I am not using sudo, even if the password is asked during the error message. I do not see it as a prompt but just in the error message. To use containr I would like to activate conda before in order to get the environment yml file.


## Post 7 by @Gael.Rossignol (2025-08-15T12:10:32.313Z)

Camille.Christe:
> `conda create -c bioconda -c cona-forge -n hybpiper hybpiper`
I think you have a typo in your command, try :
```
conda create -c bioconda -c conda-forge -n hybpiper hybpiper
```
And it will work fine.
Best regards,


## Post 8 by @Camille.Christe (2025-08-15T12:11:16.524Z)

Dear Gael,
I am sorry but I really do not understand. Here is what I have done juste before
“”"
ml Anaconda3/2024.02-1
(bamboo)-[christec@login1 ~]$ conda env list
# conda environments:
hybpiper                 /home/christec/.conda/envs/hybpiper
base                     /opt/ebsofts/Anaconda3/2024.02-1
(bamboo)-[christec@login1 ~]$ conda activate hybpiper
CondaError: Run ‘conda init’ before ‘conda activate’
(bamboo)-[christec@login1 ~]$ conda init
We trust you have received the usual lecture from the local System
Administrator. It usually boils down to these three things:
```
#1) Respect the privacy of others.
#2) Think before you type.
#3) With great power comes great responsibility.
```
[sudo] password for christec:
# >>>>>>>>>>>>>>>>>>>>>> ERROR REPORT <<<<<<<<<<<<<<<<<<<<<<
```
Traceback (most recent call last):
  File "/opt/ebsofts/Anaconda3/2024.02-1/lib/python3.11/site-packages/conda/exception_handler.py", line 17, in __call__
    return func(*args, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^
  File "/opt/ebsofts/Anaconda3/2024.02-1/lib/python3.11/site-packages/conda/cli/main.py", line 83, in main_subshell
    exit_code = do_call(args, parser)
                ^^^^^^^^^^^^^^^^^^^^^
  File "/opt/ebsofts/Anaconda3/2024.02-1/lib/python3.11/site-packages/conda/cli/conda_argparse.py", line 196, in do_call
    result = getattr(module, func_name)(args, parser)
             ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/opt/ebsofts/Anaconda3/2024.02-1/lib/python3.11/site-packages/conda/cli/main_init.py", line 155, in execute
    return initialize(
           ^^^^^^^^^^^
  File "/opt/ebsofts/Anaconda3/2024.02-1/lib/python3.11/site-packages/conda/core/initialize.py", line 139, in initialize
    run_plan_elevated(plan2)
  File "/opt/ebsofts/Anaconda3/2024.02-1/lib/python3.11/site-packages/conda/core/initialize.py", line 856, in run_plan_elevated
    result = subprocess_call(
             ^^^^^^^^^^^^^^^^
  File "/opt/ebsofts/Anaconda3/2024.02-1/lib/python3.11/site-packages/conda/gateways/subprocess.py", line 102, in subprocess_call
    stdout, stderr = process.communicate(input=stdin)
                     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/opt/ebsofts/Anaconda3/2024.02-1/lib/python3.11/subprocess.py", line 1209, in communicate
    stdout, stderr = self._communicate(input, endtime, timeout)
                     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/opt/ebsofts/Anaconda3/2024.02-1/lib/python3.11/subprocess.py", line 2088, in _communicate
    input_view = memoryview(self._input)
                 ^^^^^^^^^^^^^^^^^^^^^^^
TypeError: memoryview: a bytes-like object is required, not 'str'
```
`$ /opt/ebsofts/Anaconda3/2024.02-1/bin/conda init`
environment variables:
CIO_TEST=
CONDA_ROOT=/opt/ebsofts/Anaconda3/2024.02-1
CURL_CA_BUNDLE=
FPATH=/usr/share/lmod/lmod/init/ksh_funcs
LD_PRELOAD=
MANPATH=/opt/ebsofts/Anaconda3/2024.02-
1/share/man:/opt/ebsofts/Anaconda3/2024.02-
1/man:/usr/share/lmod/lmod/share/man::
MODULEPATH=/etc/modulefiles:/usr/share/modulefiles:/opt/modulefiles/Linux:/usr/sh
are/lmod/lmod/modulefiles/Core:/opt/ebmodules/all/Core
PATH=/opt/ebsofts/Anaconda3/2024.02-1:/opt/ebsofts/Anaconda3/2024.02-
1/sbin:/opt/ebsofts/Anaconda3/2024.02-
1/bin:/home/christec/.local/bin:/home/christec/bin:/usr/local/bin:/usr
/bin:/usr/local/sbin:/usr/sbin
PKG_CONFIG_PATH=/opt/ebsofts/Anaconda3/2024.02-1/lib/pkgconfig
REQUESTS_CA_BUNDLE=
SSL_CERT_FILE=
__LMOD_REF_COUNT_MANPATH=/opt/ebsofts/Anaconda3/2024.02-
1/share/man:1;/opt/ebsofts/Anaconda3/2024.02-
1/man:1;/usr/share/lmod/lmod/share/man:1;:1
__LMOD_REF_COUNT_MODULEPATH=/etc/modulefiles:1;/usr/share/modulefiles:1;/opt/modulefiles/Linux:1;/
usr/share/lmod/lmod/modulefiles/Core:1;/opt/ebmodules/all/Core:1
__LMOD_REF_COUNT_PATH=/opt/ebsofts/Anaconda3/2024.02-1:1;/opt/ebsofts/Anaconda3/2024.02-
1/sbin:1;/opt/ebsofts/Anaconda3/2024.02-
1/bin:1;/home/christec/.local/bin:1;/home/christec/bin:1;/usr/local/bi
n:1;/usr/bin:1;/usr/local/sbin:1;/usr/sbin:1
__LMOD_REF_COUNT_PKG_CONFIG_PATH=/opt/ebsofts/Anaconda3/2024.02-1/lib/pkgconfig:1
```
 active environment : None
   user config file : /home/christec/.condarc
```
populated config files :
conda version : 24.1.2
conda-build version : 24.1.2
python version : 3.11.7.final.0
solver : libmamba (default)
virtual packages : __archspec=1=zen2
__conda=24.1.2=0
__glibc=2.34=0
__linux=5.14.0=0
__unix=0=0
base environment : /opt/ebsofts/Anaconda3/2024.02-1  (read only)
conda av data dir : /opt/ebsofts/Anaconda3/2024.02-1/etc/conda
conda av metadata url : None
channel URLs : main/linux-64[main/linux-64](https://repo.anaconda.com/pkgs/main/linux-64)
main/noarch[main/noarch](https://repo.anaconda.com/pkgs/main/noarch)
r/linux-64[r/linux-64](https://repo.anaconda.com/pkgs/r/linux-64)
r/noarch[r/noarch](https://repo.anaconda.com/pkgs/r/noarch)
package cache : /opt/ebsofts/Anaconda3/2024.02-1/pkgs
/home/christec/.conda/pkgs
envs directories : /home/christec/.conda/envs
/opt/ebsofts/Anaconda3/2024.02-1/envs
platform : linux-64
user-agent : conda/24.1.2 requests/2.31.0 CPython/3.11.7 Linux/5.14.0-503.40.1.el9_5.x86_64 rocky/9.5 glibc/2.34 solver/libmamba conda-libmamba-solver/24.1.0 libmambapy/1.5.6
UID:GID : 209112:1000
netrc file : None
offline mode : False
An unexpected error has occurred. Conda has prepared the above report.
If you suspect this error is being caused by a malfunctioning plugin,
consider using the --no-plugins option to turn off plugins.
Example: conda --no-plugins install
Alternatively, you can set the CONDA_NO_PLUGINS environment variable on
the command line to run the command without plugins enabled.
Example: CONDA_NO_PLUGINS=true conda install
If submitted, this report will be used by core maintainers to improve
future releases of conda.
Would you like conda to send this report to the core maintainers? [y/N]:
Timeout reached. No report sent.
“”"


## Post 9 by @Gael.Rossignol (2025-08-15T12:13:19.655Z)

Ok, maybe you can remove folder “/home/christec/.conda” and try again.


## Post 10 by @Camille.Christe (2025-08-15T12:28:44.994Z)

Here is what I did just now:
rm -r .conda
ls -lah
ml Anaconda3/2024.02-1
conda create -c bioconda -c conda-forge -n hybpiper hybpiper
exit
ssh christec@login1.bamboo.hpc.unige.ch
ml Anaconda3/2024.02-1
conda env list
conda activate hybpiper
conda init
And I get the same error message.


## Post 11 by @Yann.Sagon (2025-08-15T14:51:07.493Z)

Workaround:
Add this to your .bashrc file:
```
# >>> conda initialize >>>

. /opt/ebsofts/Anaconda3/2024.02-1/etc/profile.d/conda.sh

export PATH="/opt/ebsofts/Anaconda3/2024.02-1/bin:$PATH"

# <<< conda initialize <<<
```
logout/login
and now `conda activate` should work.


## Post 12 by @Camille.Christe (2025-08-15T17:14:29.260Z)

Yann.Sagon:
```
# >>> conda initialize >>>

. /opt/ebsofts/Anaconda3/2024.02-1/etc/profile.d/conda.sh

export PATH="/opt/ebsofts/Anaconda3/2024.02-1/bin:$PATH"

# <<< conda initialize <<<
```
The problem is maybe here. I cannot write the .bashrc and change the permission


## Post 13 by @Camille.Christe (2025-08-17T13:44:28.324Z)

I erased the .bashrc and made a new one and now it is work.
Thank you very much for your help!!!


## Post 14 by @Yann.Sagon (2025-08-18T07:00:59.790Z)

Camille.Christe:
> The problem is maybe here. I cannot write the .bashrc and change the permission
Ah!!! good catch! Indeed we had several users with wrong permissions on files such as `.bashrc`.
Good to hear it is working now!
