# Home storage quota reduction on the Baobab cluster

- Source: https://hpc-community.unige.ch/t/4290

- Created: 2026-05-04T08:27:47.892Z

- Tags: baobab

- Posts: 1

- Category: 6

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yann.Sagon (2026-05-04T08:27:48.060Z)

Dear users,
Due to the continuous growth in the number of users on the Baobab cluster and current storage constraints, we will adjust the home storage quota.
:warning: This change applies to the Baobab cluster only.
Other clusters are not affected.
### :wrench: Quota change
- Current quota: 1 TB
- New quota: 0.8 TB
- Effective date: Monday, May 11, 2026
The home filesystem on Baobab is currently full. This change is necessary to ensure the stability and proper operation of the service.
### :warning: Impact on users
- Users with more than 0.8 TB of data in their home directory will no longer be able to write to this space until they free up disk space.
- Please clean up your home directory as soon as possible by removing or moving unnecessary data.
### :bar_chart: Check your quota and disk usage
You can check your home and scratch disk usage by following the official documentation:
:backhand_index_pointing_right: hpc:storage_on_hpc [eResearch Doc][hpc:storage_on_hpc [eResearch Doc]](https://doc.eresearch.unige.ch/hpc/storage_on_hpc#check_disk_usage_on_home_and_scratch)
### :pushpin: Reminder
Temporary data (intermediate results, job outputs, etc.) should preferably be stored on the `scratch` filesystem, which is significantly larger and designed for this purpose.
Thank you for your cooperation.
For any questions, please contact the HPC support team or consult the documentation.
Best regards,
The HPC team
