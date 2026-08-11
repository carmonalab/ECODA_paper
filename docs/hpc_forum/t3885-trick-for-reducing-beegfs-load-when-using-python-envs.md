# Trick for reducing BeeGFS load when using Python envs

- Source: https://hpc-community.unige.ch/t/3885

- Created: 2025-03-25T09:11:15.157Z

- Posts: 3

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Berk.Gercek (2025-03-25T09:11:15.199Z)

# Using archives for faster python loading
The core issue I used to face was that when loading a complex python environment the import statements would take up to several minutes to run. This was because of the thousands of file access to each python package file.
This solution have been mentioned before, but I’ve been using this workaround to avoid loading custom Python environments file-by-file from BeeGFS network directories. For large environments it can speed up the initial imports to run a script by an order of magnitude. There are shortcut bash functions at the end of this post as well to make it simpler.
The trick is to essentially build a virtual environment in the `/tmp` directory which uses local storage, and then package the environment directory into an archive in the home directory. You then unpack that archive (one file access, fast) to `/tmp` as the prelude to any job instead of loading the environment from the network (many file accesses, slow).
## TL;DR
- Build your environments in fast local storage at `/tmp`
- Use `7z` or another tool to archive your environment into a single file that’s fast to access
- Whenever a job runs unpack that archive to the node’s `/tmp` so you don’t need to access network files for every `import mypackage` statement
- Speeds up complex environments (>1000 python files) tremendously, reduces load on BeeGFS metadata servers by making only one request instead of >1000.
- Shortcut functions at the end for one environment or many
- Key caveats: Cannot port archives between clusters (don’t know why, path/permission issues?), archiving environments with large files like CUDA toolkit libs can make this take up lots of home directory space
## Step 1: Building a virtual environment
I’ve only tested this method with `venv` and `uv venv` (highly recommend the latter for faster package management, see Astral UV for details[see Astral UV for details](https://docs.astral.sh/uv/)) but it should in principle work with anaconda provided you specify the paths to be the same.
Replace `uv venv` with `venv` if you don’t have UV. If you use conda, you will need to create a new conda installation in `/tmp` at this stage.
You start by creating a virtual environment in the `/tmp` directory which lives in local disk storage on the cluster:
```
$ cd /tmp
$ uv venv MY_ENV_NAME -p 3.10  # This creates a python 3.10 venv, adjust as needed. Make sure not to use the default .venv name!
$ cd MY_ENV_NAME && uv pip install PKG1 PKG2  # Install as usual until your environment is  ready
```
Verify that you have all the packages you need for your cluster jobs in that environment before step 2.
## Step 2: Archiving
I use 7zip since it’s already on the cluster to archive the environment folder. This would also work with tarballs, but I know 7z better and so used that. Also for very large environments such as those using CUDA compression is nice and 7z natively supports parallel compression unlike the default bzip.
```
$ 7z a -mx0 $HOME/MY_ENV_NAME.7z /tmp/MY_ENV_NAME
```
The `a` indicator tells 7zip to append new files to the archive if it exists, and `-mx0` indicates no compression (copy mode).
## Step 3: Unpacking
Because the `/tmp` local storage directory is not shared across nodes, the archive will need to be unpacked for each job. This can be done as follows:
```
$ 7z x $HOME/MY_ENV_NAME -o /tmp/
```
After this using `source /tmp/MY_ENV_NAME/bin/activate` as usual will set up the environment correctly.
## Step 4: New packages
Installing new packages and modifying the environment is very simple: If you are only appending new packages, the command from step 2 will work again after unpack the environment to `/tmp` and installing the new packages. If you remove packages it might be necessary to delete the existing archive and build from scratch.
## Handy alias and shortcut commands: Only one environment
You can use a set of simple aliases `.bashrc` to speed up the process of doing the above for virtual environments. These work if you have a single environment that you use for everything on the cluster. If you have multiple environments you’ll need the functions in the next section. Make sure to change the environment and file names!
```
alias vpack="7z a -mx0 $HOME/MY_ENV_NAME.7z /tmp/MY_ENV_NAME/"
alias vunpack="7z x $HOME/MY_ENV_NAME.7z -o/tmp/"
alias vact="source /tmp/MY_ENV_NAME/bin/activate"
```
For Anaconda you can replace the env name with with a conda installation in `/tmp`. Be warned that if you have CUDAtoolkit installations in your conda environments using 0 compression is not a good idea! The archives will be several GB. Set `-mx` to be 1 or greater in that case, and use `-mmt` to set multiple cores for compression/extraction for speed. Even then, it still may not be an advantage to use this trick.
## Shortcuts: Multiple environments
These functions can be added to the bashrc to handle multiple environments as arguments. If you use anaconda this is not possible since all environments live in the same folder structure. In the conda case the above aliases are enough.
```
function vpack {
    if [ -n "$1" ]; then
        if [ -d "/tmp/$1" ]; then
            7z a -mx0 "$HOME/$1.7z" "/tmp/$1"
        else
            echo "Error: Directory does not exist in /tmp, cannot create archive"
            return 1
        fi
    else
        echo "Error: No virtual environment name specified."
        return 1
    fi
}

function vunpack {
    if [ -n "$1" ]; then
        if [ -f "$HOME/$1.7z" ]; then
            7z x "$HOME/$1.7z" -o/tmp/
        else
            echo "Error: Environment archive does not exist in $HOME. Cannot unpack archive."
            return 1
        fi
    else
        echo "Error: No virtual environment name specified."
        return 1
    fi
}

function vact {
    if [ -n "$1" ]; then
        if [ -f "/tmp/$1/bin/activate" ]; then
            source "/tmp/$1/bin/activate"
        else
            echo "Error: Environment does not exist in /tmp. Cannot activate. Did you unpack first?"
            return 1
        fi
    else
        echo "Error: No virtual environment name specified."
        return 1
    fi
}
```
## Using this in a job script
Below is an example script and output using the above
```
#!/bin/bash
#SBATCH --job-name=pyenv-test
#SBATCH --time=00:01:00
#SBATCH --partition=debug-cpu
#SBATCH --output=testpyenv.out

source $HOME/.bashrc  # Necessary for commands to work: script doesn't source anything normally
vunpack MY_ENV_NAME
vact MY_ENV_NAME
echo "Using python at $(which python)"
python --version
```
Hope this helps! Cluster admins, please let me know if some part of this is bad practice.


## Post 2 by @Gael.Rossignol (2025-03-27T09:37:15.876Z)

Dear Berk,
Thanks for this tutorial, I would like to suggest it’s preferable to do actions for create env on a node (add “salloc” before), to avoid to have the login node with a full “/tmp” folder.
Best regard,


## Post 3 by @Berk.Gercek (2025-03-27T11:33:17.157Z)

Very good point, thank you. A related question: How much data can/should we put in /tmp on a node? Is there a recommended upper limit, and a hard upper limit?
