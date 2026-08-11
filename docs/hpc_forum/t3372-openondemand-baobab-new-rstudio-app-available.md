# [OpenOnDemand][Baobab] New Rstudio app available

- Source: https://hpc-community.unige.ch/t/3372

- Created: 2024-03-13T16:53:33.420Z

- Posts: 3

- Category: 6

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Adrien.Albert (2024-03-13T16:53:33.512Z)

Dear HPC User,
:tada: Guess what? You won’t R’egret it! RStudio is now available on OpenOnDemand[OpenOnDemand](https://openondemand.baobab.hpc.unige.ch/pun/sys/dashboard/)!
As this is a new version, some features may be missing. We will therefore update and improve them according to your suggestions.
In the near future, we plan to make official container images available in a shared directory to avoid data replication and facilitate the creation of interactive Rstudio sessions.
What do you think?


## Post 2 by @Lucille.Delisle1 (2025-03-12T16:30:21.524Z)

Hi,
I love the RStudio OpenOnDemand. There is a new version on the BioinfoSupport repo[BioinfoSupport repo](https://github.com/BioinfoSupport/ood-rstudio-baobab) would be good to have it.
I think would be really cool to have frequently used RStudio images in a shared directory.
I don’t use official container images as I could not find an official image with all the packages I wanted but my images are versioned and documented. Would you agree to host them?
For example, this image[this image](https://hub.docker.com/repository/docker/lldelisle/verse_with_more_packages/tags/4.4.2_0/sha256-d725e481aa8893b61ac73f9a2f785a7af6f4f0234299c1272729767d123970af) documented with this Dockerfile[this Dockerfile](https://github.com/lldelisle/lldelisle-docker/blob/verse_with_more_packages_4.4.2/verse_with_more_packages/image/Dockerfile). If gitlab.unige.ch CI enables the generation and hosting of the docker image I am happy to switch.
Best,
Lucille


## Post 3 by @Adrien.Albert (2025-03-13T10:26:28.398Z)

Hi @Lucille.Delisle1[@Lucille.Delisle1](https://hpc-community.unige.ch/u/lucille.delisle1)
I am on it, I need to check the modification :slight_smile:
