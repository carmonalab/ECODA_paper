# HPC accounting summary

- Source: https://hpc-community.unige.ch/t/4056

- Created: 2025-08-21T15:37:41.601Z

- Posts: 1

- Category: 6

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yann.Sagon (2025-08-21T15:37:41.712Z)

Dear users, this is a summary of the cost model we introduced since beginning of 2025
Please refer to the following link for details regarding our updated cost model:  HPC Clusters - Cost Model[HPC Clusters - Cost Model](https://doc.eresearch.unige.ch/hpc/hpc_clusters#cost_model).
Summary:
- Starting this year, you receive a “CPU hours credit” based on the hardware you own (if any) in the cluster (private partition).
- You can find instructions on how to check your annual credit here:  Resources Available for Research Groups[Resources Available for Research Groups](https://doc.eresearch.unige.ch/hpc/accounting#resources_available_for_research_group). If you know your research group has bought some compute nodes but your PI doesn’t appear in the report, please contact us.
- The credit calculation in the provided script assumes a 5-year hardware ownership period. However, if this policy was introduced after your compute nodes were purchased, we can extend the duration until the end of their production lifecycle. To adjust this, modify the `--reference-year` argument (e.g., set it to 2022) when running the script.
- To ensure flexibility and simplicity, we have standardized resource usage by converting CPU and GPU hours into CPU hours, using different conversion ratios depending on the GPU type. More details can be found here:  Resource Accounting Uniformization[Resource Accounting Uniformization](https://doc.eresearch.unige.ch/hpc/accounting#resource_accounting_uniformization).
- You can use your credit across all three clusters (Baobab, Yggdrasil, and Bamboo), not just on your private compute nodes. However, when using your own compute nodes, you will receive a higher priority.
- To check your group’s current resource usage, visit:  Report and Statistics with sreport[Report and Statistics with sreport](https://doc.eresearch.unige.ch/hpc/accounting#report_and_statistics_with_sreport).
