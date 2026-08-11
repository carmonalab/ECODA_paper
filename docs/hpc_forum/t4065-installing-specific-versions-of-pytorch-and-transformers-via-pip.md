# Installing specific versions of pytorch and transformers via pip

- Source: https://hpc-community.unige.ch/t/4065

- Created: 2025-08-26T11:50:51.241Z

- Tags: baobab

- Posts: 11

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Erica.Lastufka (2025-08-26T11:50:51.301Z)

Hello,
I am trying to use the models and checkpoints from this repository: GitHub - facebookresearch/dinov3: Reference PyTorch implementation and models for DINOv3[GitHub - facebookresearch/dinov3: Reference PyTorch implementation and models for DINOv3](https://github.com/facebookresearch/dinov3)
With the pyTorch versions available on baobab (>2.2) I am not able to load the models. Same with the Transformers version, which is recommended to install via pip install git+https://github.com/huggingface/transformers.git
I have experimented with all the pyTorch versions available and all of them give an error when I try either to import torch, transformers, or AutoModel. Searches have said that sometimes the PyTorch version needs to be >2.2. I have tried installing pyTorch via pip and it does not seem to take effect. Similarly, if I load transformers via module load Transformers/X.X.XX, the version number after import is that one and not the latest version from the repository.
How can I use custom versions of pyTorch and transformers so that I can finally use the models?
Best regards,
~ Erica


## Post 2 by @Erica.Lastufka (2025-08-27T08:05:34.633Z)

For reference, I am able to load the models on my laptop with the following packages:
torch                     2.8.0                    pypi_0    pypi
torchmetrics              1.8.1                    pypi_0    pypi
torchvision               0.23.0                   pypi_0    pypi
transformers              4.56.0.dev0


## Post 3 by @Adrien.Albert (2025-08-27T10:09:45.196Z)

Hi @Erica.Lastufka[@Erica.Lastufka](https://hpc-community.unige.ch/u/erica.lastufka),
PyTorch is quite tricky to install, and EasyBuild , the software we use to deploy all modules, has not yet published PyTorch 2.8 recipy.
So, I think it would be interesting to try using a container instead.
## Procedure: Using cotainr[Using cotainr](https://doc.eresearch.unige.ch/hpc/applications_and_libraries#how_to_create_a_conda_environment_in_a_container)
- Create a env file (conda env export or chatgpt)
```
(baobab)-[alberta@cpu200 pytorch]$ cat pytorch_2.8.yaml 
name: pytorch_env
channels:
  - pytorch
  - conda-forge
dependencies:
  - python=3.11
  - pip
  - pip:
      - torch==2.8.0
      - torchvision==0.23.0
      - torchmetrics==1.8.1
      - transformers==4.55.0
```
Build the image:
```
ml GCCcore/13.3.0 cotainr
cotainr build pytorch-2.8.sif --base-image=docker://rockylinux:9 --accept-licenses --conda-env=pytorch_2.8.yaml
```
I have tested:
```
(baobab)-[alberta@cpu200 pytorch]$ apptainer exec pytorch-2.8.sif python3 -c "import transformers; print(transformers.__version__)"
4.55.0

(baobab)-[alberta@cpu200 pytorch]$ apptainer exec pytorch-2.8.sif  pip list
Package                  Version
------------------------ -----------
alembic                  1.13.2
antlr4-python3-runtime   4.9.3
backports.tarfile        1.2.0
banal                    1.0.6
build                    1.2.2
CacheControl             0.14.0
cachetools               5.5.0
certifi                  2024.8.30
cffi                     1.17.1
charset-normalizer       3.3.2
cleo                     2.1.0
crashtest                0.4.1
cryptography             43.0.1
dataset                  1.6.2
distlib                  0.3.8
dulwich                  0.21.7
faiss-cpu                1.8.0
fastjsonschema           2.20.0
filelock                 3.16.1
fsspec                   2025.7.0
gitdb                    4.0.11
GitPython                3.1.43
google-api-core          2.20.0
google-auth              2.35.0
google-cloud-bigquery    3.25.0
google-cloud-core        2.4.1
google-crc32c            1.6.0
google-resumable-media   2.7.2
googleapis-common-protos 1.65.0
greenlet                 3.0.3
grpcio                   1.66.1
grpcio-status            1.66.1
hf-xet                   1.1.8
huggingface-hub          0.34.4
idna                     3.10
importlib_metadata       8.5.0
installer                0.7.0
jaraco.classes           3.4.0
jaraco.context           5.3.0
jeepney                  0.8.0
Jinja2                   3.1.6
joblib                   1.4.2
keyring                  24.3.1
keyrings.alt             5.0.1
lightning-utilities      0.15.2
Mako                     1.3.5
MarkupSafe               2.1.5
more-itertools           10.3.0
mpmath                   1.3.0
msgpack                  1.1.0
networkx                 3.5
numpy                    2.0.0
nvidia-cublas-cu12       12.8.4.1
nvidia-cuda-cupti-cu12   12.8.90
nvidia-cuda-nvrtc-cu12   12.8.93
nvidia-cuda-runtime-cu12 12.8.90
nvidia-cudnn-cu12        9.10.2.21
nvidia-cufft-cu12        11.3.3.83
nvidia-cufile-cu12       1.13.1.3
nvidia-curand-cu12       10.3.9.90
nvidia-cusolver-cu12     11.7.3.90
nvidia-cusparse-cu12     12.5.8.93
nvidia-cusparselt-cu12   0.7.1
nvidia-nccl-cu12         2.27.3
nvidia-nvjitlink-cu12    12.8.93
nvidia-nvtx-cu12         12.8.90
packaging                24.1
pexpect                  4.9.0
pillow                   11.3.0
pip                      25.2
pkginfo                  1.11.1
platformdirs             4.3.6
poetry                   1.8.3
poetry-core              1.9.0
poetry-plugin-export     1.8.0
proto-plus               1.24.0
protobuf                 5.28.2
ptyprocess               0.7.0
pyasn1                   0.6.1
pyasn1_modules           0.4.1
pycparser                2.22
pyproject_hooks          1.2.0
python-dateutil          2.9.0.post0
PyYAML                   6.0.2
RapidFuzz                3.10.0
regex                    2025.7.34
requests                 2.32.3
requests-toolbelt        1.0.0
rsa                      4.9
safetensors              0.6.2
scikit-learn             1.5.1
scipy                    1.14.0
SecretStorage            3.3.3
setuptools               80.9.0
shellingham              1.5.4
six                      1.16.0
smmap                    5.0.1
SQLAlchemy               1.4.53
sympy                    1.14.0
threadpoolctl            3.5.0
tokenizers               0.21.4
tomlkit                  0.13.2
torch                    2.8.0
torchmetrics             1.8.1
torchvision              0.23.0
tqdm                     4.67.1
transformers             4.55.0
triton                   3.4.0
trove-classifiers        2024.9.12
typing_extensions        4.12.2
urllib3                  2.2.3
virtualenv               20.26.6
wheel                    0.45.1
zipp                     3.20.2
```


## Post 4 by @Erica.Lastufka (2025-08-28T08:03:28.156Z)

Thank you, Adrien! I was able to install the container. Is there a way to use this environment for debugging via IPython or JupyterLab for example? My naive attempt at using “apptainer exec pytorch-2.8.sif ipython“ yielded a fatal error


## Post 5 by @Adrien.Albert (2025-08-28T09:24:37.271Z)

Hi @Erica.Lastufka[@Erica.Lastufka](https://hpc-community.unige.ch/u/erica.lastufka),
Yes, that’s possible. Here is the step-by-step approach:
- Add jupyterlab to the list of `pip install` packages.
- Rebuild the container with cotainr.
- On a COMPUTE node, run the container with:
```
apptainer exec pytorch-2.8.sif jupyter lab --ip=$(hostname) --port=8888 --no-browser
```
- Copy the line printed by Jupyter, (for me:):
```
http://cpu185.baobab:8888/lab?token=a23859b3cb14ac924b1b3459fb80fcd54d85a1d708ce0975
```
- From your local machine, connect to the cluster with a dynamic port:
```
ssh -D XXXX user@login1.baobab.hpc.unige.ch
```
- Start Firefox on your local machine with FoxyProxy configured for the same port `XXXX`.
- Open the Jupyter URL you saved in step 4.
---
This is the handmade/DIY solution.
However, Open OnDemand can provide the same functionality in a much cleaner way. Our repository is open to the community, and we would be very grateful if you could participate in its implementation:
https://gitlab.unige.ch/hpc/unige-openondemand[https://gitlab.unige.ch/hpc/unige-openondemand](https://gitlab.unige.ch/hpc/unige-openondemand)


## Post 6 by @Erica.Lastufka (2025-08-28T09:53:06.783Z)

Hi Adrien,
In fact I prefer the great service of Open OnDemand. Is there a way to add container as a kernel for jupyterlab? I have had very mixed success with this in the past with normal conda environments, and have only got about half of them to be usable kernels.
And for running batch jobs, do I simply replace ‘srun’ with ‘apptainer exec pytorch-2.8.sif python3‘ ?


## Post 7 by @Adrien.Albert (2025-08-28T12:20:25.141Z)

As I am not a jupyter user, I do not know all the possibilities, but  That’s a good idea, this is my tests:
- generate the kernel (to build the tree)
```
(baobab)-[alberta@cpu206 pytorch]$ apptainer exec pytorch-2.8.sif bash
bash: /usr/share/git-core/contrib/completion/git-prompt.sh: No such file or directory
(conda_container_env) 
(conda_container_env) python -m ipykernel install --user --name=kernel-pytorch-2.2 --display-name "PyTorch Container"
```
- However it acts like it is not in container, you need to update the kernel file:
```
(baobab)-[alberta@cpu206 pytorch]$ ls /home/users/a/alberta/.local/share/jupyter/kernels/kernel-pytorch-2.2
kernel.json  logo-32x32.png  logo-64x64.png  logo-svg.svg

# Generate by ipykernel - NOT WORKING:
(conda_container_env) cat /home/users/a/alberta/.local/share/jupyter/kernels/kernel-pytorch-2.2/kernel.json 
{
 "argv": [
  "/opt/conda/envs/conda_container_env/bin/python",
  "-Xfrozen_modules=off",
  "-m",
  "ipykernel_launcher",
  "-f",
  "{connection_file}"
 ],
 "display_name": "PyTorch Container",
 "language": "python",
 "metadata": {
  "debugger": true
 }

# Modification to use  the conda env inside the containr:

(baobab)-[alberta@cpu206 pytorch]$ cat /home/users/a/alberta/.local/share/jupyter/kernels/kernel-pytorch-2.2/kernel.json
{
  "argv": [
    "/usr/bin/apptainer",
    "exec",
    "/home/users/a/alberta/cotainr/pytorch/pytorch-2.8.sif",
    "/opt/conda/envs/conda_container_env/bin/python",
    "-m",
    "ipykernel_launcher",
    "-f",
    "{connection_file}"
  ],
  "display_name": "PyTorch Container",
  "language": "python"
}
```
Once done; run an OpenOnDemand Jupyterlab session and you should see your kernel:
image
image1319×737 49.7 KB
[image1319×737 49.7 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/1ce8dc4f7f7b764a2cbc1d1c0b8e8615ffc2a583.png)
And to confirm list the package transformers:
image
image1501×469 35.7 KB
[image1501×469 35.7 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/74b8b36a22b53fb115d8f39e6b26a031e8c46444.png)


## Post 8 by @Erica.Lastufka (2025-08-28T14:54:09.430Z)

for me /opt/conda does not exist.
In my other kernel.json files, I have: “/opt/ebsofts/Python/3.10.4-GCCcore-11.3.0/bin/python”, “/opt/ebsofts/Python/3.9.6-GCCcore-11.2.0/bin/python3.9”, or “python”.
I tried modifying the kernel.json file as instructed, and also with “/opt/ebsofts/Python/3.10.4-GCCcore-11.3.0/bin/python” . The kernel is unable to connect.


## Post 9 by @Adrien.Albert (2025-08-29T08:46:38.387Z)

Cotainr builds a Conda environment inside a container , so you must use the Conda environment inside the container.
List the Conda environments available inside the container:
```
(baobab)-[alberta@cpu330 ~]$ apptainer exec pytorch-2.8.sif conda env list

# conda environments:
#
scgpt_conda            /home/users/a/alberta/.conda/envs/scgpt_conda
scgpt_conda_3.11.5     /home/users/a/alberta/.conda/envs/scgpt_conda_3.11.5
base                   /opt/conda
conda_container_env  * /opt/conda/envs/conda_container_env    # <========= Use this one
```
Since you have installed all your packages inside the Conda environment within the container, they are available through this environment.
You must use:
```
/opt/conda/envs/conda_container_env/bin/python
```
instead of:
```
/opt/ebsofts/Python/3.10.4-GCCcore-11.3.0/bin/python
```
Tips:
- Disable automatic activation of Conda environments in your `~/.bashrc`.
- Run `module purge` to clear all loaded modules before using the container.


## Post 10 by @Erica.Lastufka (2025-09-01T07:51:40.833Z)

I have not yet gotten this to work. Should I run module purge over ssh before connecting to OOD? How OOD is related to the rest of the system is not clear to me.
Additionally, with the container I cannot access my filesystem. I get the following errors:
File “”, line 225, in makedirs
FileNotFoundError: [Errno 2] No such file or directory: ‘/home/users/l/lastufka/scratch/GalaxyMNIST’
for folders which certainly exist.


## Post 11 by @Adrien.Albert (2025-09-01T09:01:26.142Z)

You need to use `--bind` option to mount other FS:
```
(baobab)-[alberta@cpu335 ~]$ apptainer exec pytorch-2.8.sif ls /srv/beegfs/scratch/users/a/alberta/
ls: cannot access '/srv/beegfs/scratch/users/a/alberta/': No such file or directory

(baobab)-[alberta@cpu335 ~]$ apptainer exec --bind /srv/beegfs/scratch/users/a/alberta/:/srv/beegfs/scratch/users/a/alberta/ pytorch-2.8.sif ls /srv/beegfs/scratch/users/a/alberta/
bin data sbatch
```
