# Baobab OpenOnDemand can not stat RStudio at the moment

- Source: https://hpc-community.unige.ch/t/3565

- Created: 2024-07-24T09:38:42.067Z

- Tags: baobab

- Posts: 5

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Fedor.Bezrukov (2024-07-24T09:38:42.166Z)

## Primary informations
Username: bezrukov
Cluster: baobab
## Description
The RStudio interactive app in OpenOnDemand does not work at the moment, using any Rstudio image suggested. Digging into the problem (in the output.log for hte session) I can see, that the reason is that `/acanas/m-BioinfoSupport` i snot mounted at the moment on the compute nodes.
- First, I guess that the problme will clean itself in a little while (at the login node this directory disappeared, but now I get it back mounted)
- Second, I would guess that probably adding the bind mount of `/acanas/m-BioinfoSupport` to the singularity command should be done only if the requested image resides there? Or if this directory is present and mounted. Then the app whould happily start if /acanas has problems (if I use docker based images, or my own singularity images)
## Steps to Reproduce
Start RStudio interactive App in baobab openondemand. The process completes nearly immediately and does not start the server.
## Expected Result
Expecting the app to start and the link to connect to RStudio to appear.
## Actual Result
You get the report of the form
image
image1353×471 35.6 KB
[image1353×471 35.6 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/aae2c40cf2c401c8d3461aa49829e69eaaad835b.png)
Digging into the debug information into output.log you can find that the reason is:
`FATAL:   container creation failed: mount hook function failure: mount /acanas/m-BioinfoSupport->/acanas/m-BioinfoSupport error: while mounting /acanas/m-BioinfoSupport: mount source /acanas/m-BioinfoSupport doesn't exist`


## Post 2 by @Fedor.Bezrukov (2024-07-24T09:43:52.374Z)

Oh, found the workaround :slight_smile:
Running
`srun -p public-interactive-cpu ls /acanas/m-BioinfoSupport`
trigered the automounter to mount the required directory on the compute node I was using.
However I would guess either attempting to trigger the automounter from the RStudio start script for openondemand (just by doing `ls /acanas/m-BioinfoSupport` before starting singularity) or mounting this directory optionally would improve stability of the service.
Best!


## Post 3 by @Adrien.Albert (2024-07-24T10:00:51.317Z)

Hi @Fedor.Bezrukov[@Fedor.Bezrukov](https://hpc-community.unige.ch/u/fedor.bezrukov)
The storage team performed a migration on this NASAC share. We had to unmount it on every HPC server, resulting in the unavailability of `/acanas/m-BioinfoSupport`.
Since it’s a private share, I have only contacted the owners. (My bad)
Fedor.Bezrukov:
> `srun -p public-interactive-cpu ls /acanas/m-BioinfoSupport`
It’s not a workaround; your `srun` command worked because the share has been re-mounted in the same time. :four_leaf_clover: :sparkles:
In this case, You may try to use another R image in OpenOnDemand as described in the example value.
image
image614×529 19.4 KB
[image614×529 19.4 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/0243fc7c91840e2b2d6cfc1fb2c63644ae2e7566.png)
As the share is back again, you can use /acanas/m-BioinfoSupport with R-OpenOnDemand.


## Post 4 by @Fedor.Bezrukov (2024-07-24T10:42:54.811Z)

Hi @Adrien.Albert[@Adrien.Albert](https://hpc-community.unige.ch/u/Adrien.Albert),
Thanks!
One note – while the mount was absent, the Rstudio app failed to work for any images, both those on /acanas/m-BioinfoSupport, and docker based, or stored in my home dir. The reason is that there singularity command in the app is invoked with `--bind /acanas/m-BioinfoSupport` opton, independent on what image is chosen, and not finding this directory it failed.


## Post 5 by @Adrien.Albert (2024-07-24T14:12:46.098Z)

Thanks for noticing, I will update it to make it more generic
