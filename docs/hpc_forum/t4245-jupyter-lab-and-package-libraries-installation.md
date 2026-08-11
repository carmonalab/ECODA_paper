# Jupyter Lab and package/libraries installation

- Source: https://hpc-community.unige.ch/t/4245

- Created: 2026-03-11T11:21:17.082Z

- Tags: baobab, bamboo

- Posts: 3

- Category: 8

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Sajad.Hamzenejadi (2026-03-11T11:21:17.147Z)

Hello,
I’m new to the HPC and have a few questions. I’d really appreciate your guidance.
I currently have a federated learning model written in a Jupyter notebook. From reading the forum, my understanding is that I can use Baobab and Bamboo through OOD, is that correct? If so, do you have any recommendations on how I choose between them?
Also, my federated learning model requires the following packages:
```
pip install timm
pip install torchinfo
pip install -q "flwr[simulation]" flwr-datasets
pip install thop
```
Should I install these packages directly within the JupyterLab environment, or do I need to request that they be installed by the HPC administrators?
Thank you very much for your help!


## Post 2 by @Sajad.Hamzenejadi (2026-03-16T09:14:29.301Z)

Dear @Yann.Sagon
If you have any suggestion/instructions, it is really appreciated!
Many thanks!


## Post 3 by @Yann.Sagon (2026-03-16T11:06:41.596Z)

Sajad.Hamzenejadi:
> I can use Baobab and Bamboo through OOD, is that correct?
You can use the three clusters: hpc:how_to_use_openondemand [eResearch Doc][hpc:how_to_use_openondemand [eResearch Doc]](https://doc.eresearch.unige.ch/hpc/how_to_use_openondemand#how_to_connect)
Sajad.Hamzenejadi:
> do you have any recommendations on how I choose between them
Check here:hpc:hpc_clusters [eResearch Doc][hpc:hpc_clusters [eResearch Doc]](https://doc.eresearch.unige.ch/hpc/hpc_clusters#hpc_clusters_the_clusters) but the doc is not very helpful I must admit! My advice is to use Bamboo as this is the more recent cluster with faster storage.
I’ve posted a tutorial that should do what you want:
[Tutorial] Create a Custom Kernel for JupyterLab Using cotainr[[Tutorial] Create a Custom Kernel for JupyterLab Using cotainr](https://hpc-community.unige.ch/t/tutorial-create-a-custom-kernel-for-jupyterlab-using-cotainr/4251) HPC Technical[HPC Technical](https://hpc-community.unige.ch/c/hpc-technical/5)
> This tutorial explains how to build a GPU‑ready container image using cotainr[cotainr](https://doc.eresearch.unige.ch/hpc/applications_and_libraries#conda), install PyTorch and other Python packages inside it, and expose it as a custom Jupyter kernel in an HPC environment using Open OnDemand[Open OnDemand](https://doc.eresearch.unige.ch/hpc/how_to_use_openondemand). 

1. Create the Conda Environment File (env.yml)
Create a file named env.yml with the following content: 
name: torch-gpu
channels:
  - pytorch
  - nvidia
  - conda-forge
  - defaults

dependencies:
  - python=3.10
  - pip
  - pytorch
  - torchvision
  - torchaudio
  - pytorch-cuda=1…


## Post 4 by @Yann.Sagon (2026-04-21T08:14:05.923Z)

2 posts were merged into an existing topic: [Tutorial] Create a Custom Kernel for JupyterLab Using cotainr[[Tutorial] Create a Custom Kernel for JupyterLab Using cotainr](https://hpc-community.unige.ch/t/tutorial-create-a-custom-kernel-for-jupyterlab-using-cotainr/4251/2)
