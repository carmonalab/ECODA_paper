# Support RAPIDS 24.00.00 for cuGRAPH

- Source: https://hpc-community.unige.ch/t/3856

- Created: 2025-03-11T14:04:47.603Z

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Dimitrios.Proios (2025-03-11T14:04:47.633Z)

Hi I have trouble using the most recent rapids ( which is needed for some of the latest tools i.e https://developer.nvidia.com/blog/nvidia-rapids-24-10-introduces-accelerated-networkx-with-zero-code-change-updates-for-umap-and-cudf-pandas/[https://developer.nvidia.com/blog/nvidia-rapids-24-10-introduces-accelerated-networkx-with-zero-code-change-updates-for-umap-and-cudf-pandas/](https://developer.nvidia.com/blog/nvidia-rapids-24-10-introduces-accelerated-networkx-with-zero-code-change-updates-for-umap-and-cudf-pandas/) requires >= 24.00 ) -
Thanks to @Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon)  I can run and import cuml
```
apptainer run --nv docker://rapidsai/rapidsai
python 
import cuml 
cuml.__version__
'23.08.00'
```
but not with the latest version any ideas ?
Furthermore when trying to build images locall with docker and converting them to sif and  building them directly in baobab using a gpu node fails.
However I am not that experienced with nvidia and particularities of apptainer vs docker.
If anyne manages to run the latest one, please let me know .


## Post 2 by @Yann.Sagon (2025-03-12T10:27:50.473Z)

Dear @Dimitrios.Proios[@Dimitrios.Proios](https://hpc-community.unige.ch/u/dimitrios.proios)
I have good news for you: the latest image from rapidai[the latest image from rapidai](https://docs.rapids.ai/install?_gl=1*4jfrnn*_ga*MTU5MjEyNzM5Ni4xNzQxNzczMzU1*_ga_RKXFW6CM42*MTc0MTc3MzM1NS4xLjAuMTc0MTc3MzM1NS42MC4wLjA.#selector) seems to work out of the box with and provides cuML 25!
As this image is quite recent, you can’t use older GPUs such as titan or pascal.
```
(baobab)-[sagon@login1 ~]$ salloc --gpus=1 --partition=shared-gpu --time=01:00:00 --constraint="COMPUTE_TYPE_AMPERE|COMPUTE_TYPE_TURING"
salloc: Pending job allocation 15433518
salloc: job 15433518 queued and waiting for resources
salloc: job 15433518 has been allocated resources
salloc: Granted job allocation 15433518
salloc: Nodes gpu012 are ready for job
```
Then I run the latest rapidsai image:
```
(baobab)-[sagon@gpu012 ~]$ apptainer run --nv docker://rapidsai/base:25.02-cuda12.8-py3.12
INFO:    Using cached SIF image
This container image and its contents are governed by the NVIDIA Deep Learning Container License.
By pulling and using the container, you accept the terms and conditions of this license:
https://developer.download.nvidia.com/licenses/NVIDIA_Deep_Learning_Container_License.pdf

Python 3.12.9 | packaged by conda-forge | (main, Feb 14 2025, 08:00:06) [GCC 13.3.0]
Type 'copyright', 'credits' or 'license' for more information
IPython 9.0.0 -- An enhanced Interactive Python. Type '?' for help.
Tip: Use `F2` or %edit with no arguments to open an empty editor with a temporary file.

In [1]: import cuml

In [2]: cuml.__version__
Out[2]: '25.02.01'
```
