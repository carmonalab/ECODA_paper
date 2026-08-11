# Enable CUDA in jupyterlab environment?

- Source: https://hpc-community.unige.ch/t/3733

- Created: 2024-11-19T15:04:42.276Z

- Tags: openondemand

- Posts: 4

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Chloe.Mayere (2024-11-19T15:04:42.312Z)

Hi Hpc team !
First thank you this great tool that is ondemand :slight_smile:
I am using jupyter lab to perform single cell RNA seq analysis, and there are some tools requiring gpu acceleration (scvi-tools or drvi).
I managed to install it in a venv, but I think it does some mix up with python libraries from GCCore.
More important: even if I start the session on the GPU partition with one GPU. The GPU is detected but in the end I get an error message saying that torch is not compiled with CUDA enabled.
Is there an option I am missing ? Is there a GCCore environment available were CUDA is enabled ?
Otherwise, I see that it is possible to run scvi tools with a docker image (https://hub.docker.com/r/scverse/scvi-tools[https://hub.docker.com/r/scverse/scvi-tools](https://hub.docker.com/r/scverse/scvi-tools)). Is it possible to have jupyterlab running through this (as for RStudio ?)
Thank you very much, hope it is not to confuse!
Best,


## Post 2 by @Yann.Sagon (2024-11-20T15:14:56.343Z)

Dear @Chloe.Mayere[@Chloe.Mayere](https://hpc-community.unige.ch/u/chloe.mayere)
I have a solution for you! I’ve added an extra field in the JupyterLab submission form to specify extra modules you want to load, see in the screenshot below.
firefox_IYgeGoiLne
firefox_IYgeGoiLne710×349 10.8 KB
[firefox_IYgeGoiLne710×349 10.8 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/cc099152eac4627bf32b73fb92bbaca08698a595.png)
Best
Yann


## Post 3 by @Chloe.Mayere (2024-11-22T09:23:35.823Z)

Thank you so much !
I’ll give it a try straight away :blush:
Best !


## Post 4 by @Chloe.Mayere (2024-11-22T10:07:19.314Z)

Another question, where can I see which modules are available?
:slightly_smiling_face:
Do I still have to install these in my env (or venv) ?
(Sorry for these naive questions ^^)
