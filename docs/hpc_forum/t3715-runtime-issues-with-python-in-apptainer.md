# Runtime issues with Python in Apptainer

- Source: https://hpc-community.unige.ch/t/3715

- Created: 2024-11-08T14:15:29.606Z

- Tags: baobab

- Posts: 4

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Alexander.Froch (2024-11-08T14:15:29.650Z)

## Primary informations
Username: froch
Cluster: baobab
## Description
I’m currently running into some runtime issues when executing one of my scripts within one of my apptainer containers. When trying out the same script with a virtual env, it works but inside the container, it doesn’t.
## Steps to Reproduce
I’m running my container using the following command:
source /home/users/f/froch/software/scripts/Filespaths.txt
apptainer exec --contain --pwd ${PWD} -B /home -B /tmp -B /srv
docker://gitlab-registry.cern.ch/alfroch/root_ml:latest bash
Afterwards, I start my script with simply running:
python get_matching_eff.py -e ttbar -n 100
## Expected Result
The script should output two files and it seems to be unable to load the root file from scratch (based on some debug statements I put).
## Actual Result
I don’t get any output. The script tries to load the root file but I can’t. When I try to cancel the script, also nothing happens and I need to shut it down via another terminal and htop.


## Post 2 by @Gael.Rossignol (2024-11-15T16:09:30.728Z)

Dear Alexander,
There are some missing informations in your post to help me find a way to resolve issue.
I reproduce your command in my home directory by updating of course “Filespaths.txt” like that
```
PWD="/home/users/r/rossigng/"
PACKAGES="${PWD}/software/"
SCRIPTSDIR="${PACKAGES}/scripts"
SUBMISSION_SCRIPTS="${PACKAGES}/scripts"
UMAMIDIR="${PACKAGES}/umami"
SALTDIR="${PACKAGES}/salt"
PUMADIR="${PACKAGES}/puma"
UPPDIR="${PACKAGES}/umami-preprocessing"
FTAGTOOLSDIR="${PACKAGES}/atlas-ftag-tools"
```
I have done a copy of your software folder in my home previously and launched command
```
source /home/users/r/rossigng/software/scripts/Filespaths.txt
apptainer exec --contain --pwd ${PWD} -B /home -B /tmp -B /srv docker://gitlab-registry.cern.ch/alfroch/root_ml:latest bash
INFO:    Using cached SIF image
Apptainer>
```
First , script is launched on the login node and not on compute but you may have a sbatch file you didn’t provide.
Anyway, when I am in the apptainer I have done some checks
```
Apptainer> python --version
Python 3.10.12
Apptainer> find . -name get_matching_eff.py
Apptainer>
```
Then python is available on the container but the script “get_matching_eff.py” is not found.
Do you know where is this script?
Best regards,


## Post 3 by @Alexander.Froch (2024-11-25T11:14:54.552Z)

Hi Gael,
sorry for the late respone.
The script is usually run in an sbatch job, that’s true. I didn’t provide it, because the problem arose also in interactive jobs (where I usually test/run short stuff).
I traced down the issue a bit further. It doesn’t seem to be a direct issue with python, but more with accessing files. I’m using Uproot to load a root file from my scratch storage. While in a virtual env everything works well, inside the container the script dies when trying to load the root file.
The script btw. can be found here[here](https://gitlab.cern.ch/alfroch/matching_efficiency/-/blob/master/scripts/get_matching_eff.py?ref_type=heads)


## Post 4 by @Victor.Ferat (2024-12-06T10:54:01.906Z)

The scratch directory is a symlink and is not mounted when mouting the /home directory. You might need to mount the symlink location in addition to your current mount commands:
```
 --bind /srv/beegfs/scratch/users/${USER:0:1}/${USER}:/srv/beegfs/scratch/users/${USER:0:1}/${USER}
```
