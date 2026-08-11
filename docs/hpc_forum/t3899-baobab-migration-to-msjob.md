# Baobab migration to MSJob

- Source: https://hpc-community.unige.ch/t/3899

- Created: 2025-04-01T08:40:55.830Z

- Tags: baobab

- Posts: 2

- Category: 6

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yann.Sagon (2025-04-01T08:40:55.942Z)

Dear Users,
We are proud to announce that Baobab has successfully migrated to a corporate operating system: Microsoft Windows. As part of this transition, we have replaced Slurm with MSJob, our new job scheduler.
Please note the following important changes:
- No more Slurm – All job submissions should now be done via MSJob. Please refer to the updated documentation for usage guidelines.
- Line endings matter – When writing scripts, ensure that lines end with `\r` instead of `\n` to maintain compatibility with Windows.
- Executables must be Windows-compatible – Some software compiled for Linux may not work. If you encounter errors, please recompile your software for Windows.
- Files are now stored directly in the cloud – Make sure to adjust your workflows accordingly, especially when dealing with file access and storage performance.
We appreciate your cooperation during this transition. If you have any questions or encounter issues, feel free to reach out to our support team.
Some screenshot example:
Connection to the cluster using `msconnect`:
```
Microsoft Windows [Version 10.0.19045.2728]  
(c) Microsoft Corporation. All rights reserved.  

C:\Users\Admin> echo Welcome to our new HPC service  
Welcome to our new HPC service  

C:\Users\Admin>_
```
Here’s an example of how you might submit a job using `MSJob` in the Windows terminal:
```
C:\Users\Admin> msjob submit --name=my_job --script=run_script.bat --time=01:00:00 --nodes=2 --cpus-per-node=4
Error: Executable format not supported by Windows. Please recompile your software for compatibility.

C:\Users\Admin> msjob submit --name=my_job --script=run_script.exe --time=01:00:00 --nodes=2 --cpus-per-node=4
Job my_job submitted successfully with Job ID: 12345

C:\Users\Admin> msjob status 12345
Job ID: 12345
Status: Running
Nodes: 2
CPUs per Node: 4
Time Remaining: 00:45:12

C:\Users\Admin> msjob cancel 12345
Job 12345 has been canceled.

C:\Users\Admin>_
```
Best regards, you HPC Team


## Post 2 by @Bharathkumar.Radhakrishnan (2025-04-01T11:51:27.376Z)

I was just getting riled up to write an angry post asking why we were not warned of this beforehand. You guys got me, I should have checked the date :joy:
