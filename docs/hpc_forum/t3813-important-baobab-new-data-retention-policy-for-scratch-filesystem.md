# [Important][Baobab] New Data Retention Policy for Scratch Filesystem

- Source: https://hpc-community.unige.ch/t/3813

- Created: 2025-02-06T09:11:34.139Z

- Posts: 1

- Category: 6

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Adrien.Albert (2025-02-06T09:11:34.245Z)

Dear HPC Users,
We would like to inform you about the implementation of a new data retention policy concerning the scratch filesystem on our HPC clusters.
### Reminder:
The scratch storage is not intended for long-term data storage. Its primary purpose is to handle the large volumes of data generated and processed during HPC computations.
### New Policy:
:triangular_flag_on_post:  EDIT:
```
- Starting from baobab-scheduled-maintenance-18-21-february-2025
+ Starting from March 18 2025, after  hpc-lunch-13th-of-march-2025/ dedicated to data retention policy and cost model
```
Starting from March 18 2025, after [HPC][ON-SITE session] hpc-lunch 13th of March 2025[[HPC][ON-SITE session] hpc-lunch 13th of March 2025](https://hpc-community.unige.ch/t/hpc-on-site-session-hpc-lunch-13th-of-march-2025/3826/1) dedicated to data retention policy and cost model , any data on the scratch filesystem that is older than one year will be automatically deleted. This policy will be progressively adjusted up to 3 months.
- Deletion is based on the last access (read or writte) date of each file.
### Deployment Plan:
This new policy will be initially deployed on Baobab and progressively rolled out to all other HPC clusters. We will keep you informed about the specific timelines for each cluster  and policy developments.
image
image1551×460 42.4 KB
[image1551×460 42.4 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/fccc26b5fb7e6069f2b7aaadf0511f4bafe1be0a.png)
### Why This Change?
- Freeing Up Space: This policy will help purge cold or inactive data, ensuring that new data can be stored without the risk of reaching storage limits.
- Improved Performance: Removing outdated files reduces metadata overhead, leading to faster file system operations and improved I/O performance.
- Efficient Resource Management: Helps maintain optimal disk usage, reducing the likelihood of storage bottlenecks during peak computational periods.
- Reduced Maintenance Efforts: Simplifies system administration, allowing for quicker troubleshooting and system upgrades.
## What You Should Do:
- Read the documentation[documentation](https://doc.eresearch.unige.ch/hpc/storage_on_hpc#data_retention_policy) about it.
- Regularly review and clean up your data.
- Move any important long-term data to appropriate storage solutions. Unige provides Storage Solutions[Storage Solutions](https://www.unige.ch/researchdata/en/store/storage-unige/).
- Automate your workflows to account for this new policy.
---
For any questions or concerns, do no hesitate to ask in HPC-community.
Thank you for your cooperation.
Best regards,


## Post 2 by @Adrien.Albert (2025-02-10T16:34:23.803Z)




## Post 3 by @Gael.Rossignol (2025-02-18T06:53:35.658Z)
