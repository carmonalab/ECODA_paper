#  Feasibility of Deploying Large Language Models Locally on Bamboo Cluster

- Source: https://hpc-community.unige.ch/t/3989

- Created: 2025-06-25T07:08:55.817Z

- Tags: bamboo

- Posts: 2

- Category: 8

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Hanzhang.Zheng (2025-06-25T07:08:55.924Z)

Dear HPC Community,
I’m currently working on a research project involving large-scale legal text processing using court verdicts.
We are considering deploying a large language model (e.g., Qwen, DeepSeek, or a quantized version of LLaMA) locally on the Bamboo cluster to perform batch inference on a large corpus (over 2 million documents).
Before proceeding, I would like to ask:
Is it feasible to locally deploy a large model (e.g., 7B parameters, quantized) on the A100 GPU nodes in Bamboo?
Are there any storage, containerization (e.g., Singularity), or dependency limitations we should be aware of when preparing the environment?
Has anyone in the community successfully run similar workloads (LLMs, inference pipelines) on Bamboo?
Would the HPC team recommend any specific best practices for running such large models efficiently in this cluster?
Any advice, shared experience, or technical pointers would be highly appreciated!
Hanzhang


## Post 2 by @Jamil.Zaghir (2025-06-25T07:42:56.160Z)

Hello,
What tools are you planning to use to deploy your LLM?
Here is a post / tutorial on how to deploy Llama using Ollama on Baobab: Seeking Help with Running Ollama on HPC Clusters[Seeking Help with Running Ollama on HPC Clusters](https://hpc-community.unige.ch/t/seeking-help-with-running-ollama-on-hpc-clusters/3563)
I assume the methodology should be similar with Bamboo.
Best regards,
Jamil
