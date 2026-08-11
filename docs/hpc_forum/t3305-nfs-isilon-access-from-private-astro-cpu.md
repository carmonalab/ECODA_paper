# NFS-Isilon access from private-astro-cpu

- Source: https://hpc-community.unige.ch/t/3305

- Created: 2024-02-12T12:09:08.207Z

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Max.Briel (2024-02-12T12:09:08.267Z)

## Primary information
Username: briel
Cluster: yggdrasil
## Description
Submitting a job array that writes to `nfs-isilon.astro.unige.ch:/ifs/astro/projects/posydon` results in the jobs immediately FAILING with an exit code 0:53, when ran from the private-astro-cpu partition.
I’m still able to access and create files on this filesystem from the login node.
Moreover, it’s also still possible to write to `nfs-isilon.astro.unige.ch:/ifs/astro/projects/posydon from the `debug-cpu`.
These jobs ran from the private-astro-cpu partition used to work fine, but now suddenly fail.
## Steps to Reproduce
- Create a slurm submission file with the `private-astro-cpu` partition on `nfs-isilon.astro.unige.ch:/ifs/astro/projects/posydon` and submit it from there.
- The submitted run fails.
Example
```
#!/bin/bash
#SBATCH -N 1
#SBATCH --array=0
#SBATCH --partition=private-astro-cpu
#SBATCH --ntasks-per-node 1
#SBATCH --time=0-00:05:00
#SBATCH --job-name="mesa_grid_\${SLURM_ARRAY_TASK_ID}"
#SBATCH --output=mesa_grid.%A_%a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=max.briel@unige.ch

echo "test"
```
## Expected Result
I expect the test submit to create a file with `mesa_grid_SLURM_ARRAY_TASK_ID.out` with ‘test’ written in it.
## Actual Result
I get no `.out` file and the job nearly immediately fails with a 0:53 error code.


## Post 2 by @Matthias.Kruckow (2024-02-14T13:15:10.952Z)

It does work again after the nfs-isilon got remounted on the nodes by Remy.
