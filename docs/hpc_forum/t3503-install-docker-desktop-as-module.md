# Install Docker Desktop as module

- Source: https://hpc-community.unige.ch/t/3503

- Created: 2024-06-21T09:15:38.169Z

- Posts: 3

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Njiva.Andrianarivelo (2024-06-21T09:15:38.208Z)

Hi,
Thanks for your amazing work maintaining the HPC.
I’m sorting tons of electrophysiology data using spikeinterface that allows to run different spike sorters at the same time on the same data.
https://spikeinterface.readthedocs.io/en/latest/modules/sorters.html#containerizedsorters[https://spikeinterface.readthedocs.io/en/latest/modules/sorters.html#containerizedsorters](https://spikeinterface.readthedocs.io/en/latest/modules/sorters.html#containerizedsorters)
To have an easy access to the spike sorters I would need to install Docker desktop.
Docker Desktop: The #1 Containerization Tool for Developers | Docker[Docker Desktop: The #1 Containerization Tool for Developers | Docker](https://www.docker.com/products/docker-desktop/)
Do you think it’s possible to have it as a module ?
Thanks,
Andry


## Post 2 by @Adrien.Albert (2024-06-21T10:35:23.589Z)

hi @Njiva.Andrianarivelo[@Njiva.Andrianarivelo](https://hpc-community.unige.ch/u/njiva.andrianarivelo)
Docker-desktop is not available with module and Package for RHEL/Rocky is in Early Access and seems to require  License for “business” use:
Docker Documentation – 13 Jun 24[Docker Documentation – 13 Jun 24](https://docs.docker.com/desktop/install/rhel/#install-docker-desktop)
### Install Docker Desktop on RHEL[Install Docker Desktop on RHEL](https://docs.docker.com/desktop/install/rhel/#install-docker-desktop)
Instructions for installing Docker Desktop on RHEL
image
image928×614 31.4 KB
[image928×614 31.4 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/6a7dfefac917557938c19bfcfb4816f3a3e7f630.png)
More information about this issue, on docker forum:
Docker Community Forums – 12 Jun 24[Docker Community Forums – 12 Jun 24](https://forums.docker.com/t/installation-instructions-for-docker-desktop/141939)
### Installation instructions for docker desktop[Installation instructions for docker desktop](https://forums.docker.com/t/installation-instructions-for-docker-desktop/141939)
Docker Desktop for Linux
The instructions here https://docs.docker.com/desktop/install/rhel/  say to Install Docker Desktop,    Set up Docker’s package repository as follows:  $ sudo dnf config-manager --add-repo https://download.docker.com/linux/rhel/docker-ce.repo   ...
At this point, there’s not much we can do on our side.


## Post 3 by @Njiva.Andrianarivelo (2024-06-21T11:04:56.494Z)

No problem, thanks for looking guys.
