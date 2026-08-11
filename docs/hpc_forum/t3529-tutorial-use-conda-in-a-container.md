# [tutorial] use conda in a container

- Source: https://hpc-community.unige.ch/t/3529

- Created: 2024-07-08T11:47:48.524Z

- Posts: 1

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yann.Sagon (2024-07-08T11:47:48.627Z)

Dear users, many HPC users use Conda[Conda](https://www.anaconda.com/) to install the required software on the cluster and/or on their laptop.
The advantage of using Conda for the end user is great flexibility and easy software installation.
The downside is that each user installs their own copy of the software in their home or scratch directory, and a typical Conda installation creates thousands of files and puts a lot of pressure on our shared storage.
The good news: it is quite easy to create a container, in our case with [Apptainer] (hpc:applications_and_libraries [eResearch Doc][hpc:applications_and_libraries [eResearch Doc]](https://doc.eresearch.unige.ch/hpc/applications_and_libraries#apptainer_was_singularity)), that contains all the software you need. The big advantage is that there is only one file, easier to maintain and co use on all the clusters.
In this tutorial, we’ll use a project called cotainr[cotainr](https://cotainr.readthedocs.io/en/stable) for this, because one of the use cases[use cases](https://cotainr.readthedocs.io/en/stable/user_guide/conda_env.html#conda-environments) of this script is to create a container from a yaml environment file that you would use to install your Conda stuff.
- Download cotainr:
```
wget https://github.com/DeiC-HPC/cotainr/archive/refs/tags/2023.11.0.tar.gz
```
- Extract the files `tar xf 2023.11.0.tar.gz`
- Make a symlink to your `~/bin` folder so you can execute it without specifying the full path
```
mkdir ~/bin
cd ~/bin
ln -s ../cotainr/cotainr-2023.11.0/bin/cotainr .
```
- To use cotainr, you need Python 3.8.0 or higher, so let’s load this Python version using module[module](https://doc.eresearch.unige.ch/hpc/applications_and_libraries#module_-_lmod): `ml GCCcore/13.2.0 Python/3.11.5`.
- You will need to create an environment file: `my_conda_env.yml` which describe what you want to be installed in the container
```
channels:
  - conda-forge
Dependencies:
  - python=3.11.0
  - numpy=1.23.5
```
Now you are ready to build your container:
```
cotainr build my_conda_env_container.sif --baseimage=docker://ubuntu:22.04 --accept-licences --conda-env=my_conda_env.yml
```
- You can now use the container:
```
apptainer exec my_conda_env_container.sif python3 -c "import numpy; print(numpy.__version__)"
1.23.5`
```
Even if `cotainr` doesn’t directly support building an image from a `pip` requirements file, you can add the pip packages you need to your conda env file.
```
channels:
  - conda-forge
dependencies:
  - python=3.11.0
  - numpy=1.23.5
  - pip
  - pip:
    - scipy==1.9.3
```
If you convert an existing Conda environment using cotainr, do not forget to erase the `.conda` files not needed anymore.
Your feedback is greatly appreciated!
