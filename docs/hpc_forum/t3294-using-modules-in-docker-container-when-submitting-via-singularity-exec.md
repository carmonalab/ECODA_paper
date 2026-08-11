# Using modules in docker container when submitting via `singularity exec`?

- Source: https://hpc-community.unige.ch/t/3294

- Created: 2024-02-07T17:08:41.724Z

- Posts: 3

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Imahn.Shekhzadeh1 (2024-02-07T17:08:41.781Z)

Hi HPC Team,
So far, I have run a Python script via a bash submit file, where I loaded a conda environment. Here the relevant snippet:
```
module purge && module load GCC/11.3.0 OpenMPI/4.1.4 geopsy/3.4.2 Anaconda3/2022.05
source activate <conda-env>
python -B run.py
```
In `run.py`, the `geopsy` module is needed.
Now I wanted to run this in a docker container, so I’m doing
```
module purge && module load GCC/11.3.0 OpenMPI/4.1.4 geopsy/3.4.2 Anaconda3/2022.05
srun singularity exec --nv $HOME/docker/migrate_v002.sif python run.py
```
and the script doesn’t run through with the following error:
```
FileNotFoundError: [Errno 2] No such file or directory: '/opt/ebsofts/geopsy/3.4.2-foss-2022a/bin/gpdc'
```
My suspicion is that the `gpdc` module isn’t loaded (even though I wrote the `module load` commands before `srun singularity exec`). How can `gpdc` be made available in the docker container? :slightly_smiling_face:
The path to the full submit file is here:
```
/home/users/s/shekhza2/ant-migrate/ml4ant/ant/fwd_model/submit_files/dgp/07p02p24/17p11/run_dgp_17p11p41_190.sh
```
Best regards,
Imahn


## Post 2 by @Adrien.Albert (2024-02-08T14:50:43.928Z)

Hi @Imahn.Shekhzadeh1[@Imahn.Shekhzadeh1](https://hpc-community.unige.ch/u/imahn.shekhzadeh1)
/opt/ebsoft is not “mounted” in your container your need to bind it.
https://apptainer.org/user-docs/master/bind_paths_and_mounts.html[https://apptainer.org/user-docs/master/bind_paths_and_mounts.html](https://apptainer.org/user-docs/master/bind_paths_and_mounts.html)


## Post 3 by @Imahn.Shekhzadeh1 (2024-02-15T15:48:19.821Z)

Hi @Adrien.Albert[@Adrien.Albert](https://hpc-community.unige.ch/u/Adrien.Albert),
If I can come back to this question: While binding `/opt/ebsoft` in the docker container was indeed the solution, it turns out I’m not able to run `geopsy` in the docker container (even though I am with a `conda` environment on the cluster). I can tell that `geopsy` is not run, since the module takes an input file and spits out an output file. However, the output files are empty (when running this in the conda environment, they are not).
In my scripts, I use the path to `geopsy` as follows:
```
path__to_gpdc = "/opt/ebsofts/geopsy/3.4.2-foss-2022a/bin/gpdc"
proc = sp.Popen(
[
       path__to_gpdc,
       "-group",
       "-f",
       [...],
       input_vel_model,
],
       stdout=out_file,
       stdin=sp.PIPE,
       stderr=sp.PIPE,
       text=True,
)
```
I ran the singularity container on `yggdrasil` as follows:
```
srun singularity exec --bind /opt/ebsofts:/opt/ebsofts --env PATH="/opt/gpdc/bin:$PATH" --nv $HOME/docker/migrate_v002.sif python -m ant.fwd_model.generate_data [...]
```
Any hints on what might be going wrong?
Thanks in advance.
