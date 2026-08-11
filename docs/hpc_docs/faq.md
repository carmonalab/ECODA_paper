# Source: https://doc.eresearch.unige.ch/hpc/FAQ
# Snapshot: 2026-08-11
# Crawled: 2026-08-11T14:32:28Z

---

~~QNA~~

## FAQ: Frequently Asked Question

#### General
**Q:** I have a problem, what should I do?
**A:** Please follow these steps:

- Review this FAQ to see if your issue is addressed.
- Check the current issues on the cluster here: https://hpc-community.unige.ch/t/2024-current-issues-on-hpc-cluster/ (A new post is created each year for reference).
- Post in the [HPC-community](https://hpc-community.unige.ch/c/hpc-support/hpc-issues/) under the category **HPC issue > HPC support** using **the Template**.

**Q:** Which cluster should I use ?
**A:** You can use the three clusters, but see [this link](https://doc.eresearch.unige.ch/hpc/hpc_clusters#hpc_clusters_the_clusters) to help you choose the right cluster.

**Q:** Must I include citations and acknowledgments in my publication?
**A:** Yes, according the [terms of use](https://www.unige.ch/eresearch/en/services/hpc/terms-use/) you **must** include at least:

"The computations were performed at University of Geneva using Baobab HPC service."

**Q:** Why is the cluster running slowly ?

!!!There could be several reasons for the cluster to slow down. It’s important to figure out where the slowness is happening:

- **Login Node**:The login node is designed for light tasks such as file editing, job submission, and monitoring, not for running heavy computations. To ensure fair usage and maintain responsiveness, each user is limited to 2 CPU cores and 8 GB of RAM on the login node

- **Compute Nodes**: Slowness on the compute nodes might be due to high CPU usage, storage issues, or other factors, which could cause your jobs to run more slowly.

- **Storage (Home, Scratch, Other)**: If there’s a problem with storage (like home directories or scratch space), it can slow down the entire cluster and affect your job performance.

**What You Can Do**:
Make sure you’re not contributing to the slowdown. Use the `htop` command on the login node to check CPU usage. If you see that all the CPUs are in use, take a screenshot and send it to us at [https://doc.eresearch.unige.ch/hpc@unige.ch](https://doc.eresearch.unige.ch/hpc@unige.ch) so we can look into it.

#### Cost
**Q:** What do you mean by 'billing metric' when discussing HPC usage and cost?
**A:** Due to the significant variation in GPU models and RAM usage, we now apply a uniform resource accounting system by converting GPU hours and memory usage into CPU-hour equivalents. This metric is called [billing](https://doc.eresearch.unige.ch/accounting#resource_accounting_uniformization).

**Q:** I'm using Open XDMoD to check how many hours my group used this year, but the amount is quite different from what appears on the invoice or in `ug_slurm_usage_per_user.py`.
**A:** For invoicing, we use a Slurm metric called billing. This metric aggregates CPU hours, memory usage, and GPU usage. Unfortunately, Open XDMoD does not currently support this metric and therefore does not take into account GPU type or memory usage. We have added a warning about this in [Open XDMoD documentation](https://doc.eresearch.unige.ch/accounting#openxdmod)

**Q:** My group only used our private partition. I heard from a colleague that it is free of charge, yet I still received an invoice.
**A:** Because GPU models and RAM usage vary significantly, we now apply a uniform resource accounting system by converting GPU hours and memory usage into CPU-hour equivalents.
This metric is called billing. If you own a private partition, you receive a fixed annual allocation of billing units that can be used on any cluster. More details are available here: [Resources available for research groups](https://doc.eresearch.unige.ch/accounting#resources_available_for_research_group)

**Q:** Your documentation states the cost per hour is 0.0157 CHF, but my invoice shows 0.02 CHF.
**A:** The price per hour shown on the invoice is rounded for display purposes only. The actual calculation uses the full precision of the rate.

**Q:** I'm a PI and want to see usage details per user.
**A:** You can use the script `ug_slurm_usage_per_user.py` for this purpose. An example is provided in our documentation: Aggregate usage by all users of a given PI. [Aggregate usage by all users of a given PI](https://doc.eresearch.unige.ch/accounting#aggregate_usage_by_all_users_of_a_given_pi)

**Q:** I have no idea why I received your email about 'HPC billing'.
**A:** The message is about the fact that the high performance computing serice known as Baobab will become a paid service after a free quota has been used. We sent the announcement to two mailing lists:
- baobab-announce: which includes all users of the Baobab service.
- hpc-community: very low-traffic mailing list containing all PIs and people interested in the HPC community. It may happen that you belong to the two mailings.

**Q:** I'm not interested in receiving further information about HPC at UNIGE, can you please remove me from the hpc-community mailing list?
**A:** If you are a UNIGE member or have a [switcheduid](https://eduid.ch/switcheduid) account, you can unsubscribe from the "hpc-community" list on [sympa web interface](https://listes.unige.ch/sympa/signoff/hpc-community?previous_action=review).

An alternate method is to send an email to [https://doc.eresearch.unige.ch/sympa@listes.unige.ch](https://doc.eresearch.unige.ch/sympa@listes.unige.ch) with the following mail body "UNSUBSCRIBE hpc-community". This mail must be sent using the email you wish to unsubscribe from.

If you are not a UNIGE member or if none of the previous steps worked, please send a request to [https://doc.eresearch.unige.ch/hpc@unige.ch](https://doc.eresearch.unige.ch/hpc@unige.ch), subject: "please unsubscribe me from the hpc-community mailing list".

Please note that you can't unsubscribe from the "baobab-announce" list if you still have an account on the Baobab.

**Q:** I'm a PI, how do I know which users are associated with me on Baobab?
**A:** If you have access to one of the clusters, you can use the `sshare` command:
```
(baobab)-[root@admin1 ~]$ sshare  -a -A <your_isis_username>
Account                    User  RawShares  NormShares    RawUsage  EffectvUsage  FairShare
-------------------- ---------- ---------- ----------- ----------- ------------- ----------
isis_pi                                 41    0.014594    73169235      0.031775   0.221089
 isis_pi                  user1          1    0.000768      130935      0.000239   0.805648
 isis_pi                  user2          1    0.000768     5069653      0.000300   0.762562
 isis_pi                  user3          1    0.000768           0      0.000000   1.000000
 isis_pi                  user4          1    0.000768           0      0.000000   1.000000
 isis_pi                  user5          1    0.000768     1707102      0.000285   0.773432
 [...]

```

You can also use [OpenXDmoD](https://doc.eresearch.unige.ch/hpc/accounting#job_accounting) to check user usage. Note that the list may be incomplete: for example, if a registered user has never used the cluster in the time period you specify, they won't appear at all.

**Q:** I'm a faculty/group manager, how may I have a list of every PI of a given dept?
**A:** You can use sacctmgr for that purpose
```
sacctmgr show assoc where parent=<your_deptartment_name> cluster=baobab format=account

```
If you don't know the name of your departement as registered in our cluster, you can list them by faculty:
```
sacctmgr show assoc where parent=sciences cluster=baobab format=account
   Account
----------
     astro
      biad
     biani
     bicel
     [...]

```

**Q:** If multiple PIs jointly purchase computing resources on the Baobab cluster, who receives the invoice?

**A:** When several Principal Investigators (PIs) collaborate to acquire computing resources on the Baobab cluster, the invoice will be sent to the designated contact person of the group associated with the partition named private_xxx, where xxx is the name of the partition. This person acts as the reference for the group and is responsible for managing the billing and communication related to the shared resources.

It is up to the PIs involved to agree among themselves on how the computing hours are distributed. The Baobab team does not manage internal allocations within shared purchases.

**Q:** I'm a PI, I tried to use OpenXDmoD to see the past usage of my group without success
**A:** We have a [tutorial](https://hpc-community.unige.ch/t/tutorial-see-your-past-computation-usage-using-openxdmod/3130) which explain how to do that.

**Q:** With OpenXDmoD how can I check usage on more than one partition?
**A:** Unfortunately, it seems that you need to do this operation for each partition separately.

**Q:** I want to login to OpenXDmoD, what are the login details?
**A:** User authentication isn't available at the moment. You can access all metrics without authentication. In the future, you'll be able to connect using your [switcheduid](https://eduid.ch/switcheduid) credentials, with the benefit of being able to create custom dashboards.

**Q:** I'm a user and I've noticed that I'm connected to two PIs, how is this possible?
**A:** The PI must be seen as a project. You can be part of two projects, and when you submit a job to the cluster, you can specify which project to charge to using the `--account` flag.
You can check who is your PI and which is your default account with the following command:
```
sacctmgr show user <your_username> -s format=user,account,defaultaccount cluster=baobab
      User    Account   Def Acct
---------- ---------- ----------
  username acct1       xxx
  username acct2       xxx

```

**Q:** I'm organizing a course and we need some HPC resources for the students. Do we have to pay for it?
**A:** The Baobab service is free for courses as long as the usage is low and for a defined period of time. Check [How our clusters work](https://doc.eresearch.unige.ch/hpc/hpc_clusters#use_baobab_for_teaching).

#### Account
**Q:** Who can be registered as PI (Primary Investigator) at Baobab?

!!!

- Anyone with an ISIs account (even an external one) with a fairly long validity period: teachers rather than assistants.
- Someone who knows the users they are going to invite. All their users will be under their responsibility.
- Someone who knows the service(s) for which they will be inviting outsiders. HPC in this case
- Person responsible for teaching or research (therefore responsible for the data generated and able to know what to do with it when a user leaves).

???When does my account expire ?

**A:** * If you have a non student account (Phd, postdoc, researcher), your account will expire at the same time your contract expire at UNIGE. Right now, there is a grace period after the end of your contract of around 6 months.
- If you have an outsider account, you need to check the expiration date you received when you filled the invitation.
- If you have an UNIGE student account, you can check the expiration date with the `chage` command:
```
(baobab)-[yourusername@login2 ~]$ chage -l yourusername
Last password change                                    : Apr 01, 2022
Password expires                                        : never
Password inactive                                       : never
Account expires                                         : never
Minimum number of days between password change          : 0
Maximum number of days between password change          : 99999
Number of days of warning before password expires       : 7

```

**Q:** I'm leaving UNIGE, can I continue to use Baobab HPC service?
**A:** For UNIGE or external members, please first refer to the guidelines about account expiration and the grace period: https://plone.unige.ch/distic/pub/isis/comment-fonctionne-isis#prolongations.
Please note that expired accounts are eligible for deletion at any time. We strongly recommend that you carefully prepare for your departure or contract extension.

However, it is possible to extend your access as long as you maintain close collaboration with your former research group. Your PI must [invite](https://gestion-externe.unige.ch/main/outsider-requests) you as an [outsider](https://doc.eresearch.unige.ch/hpc/access_the_hpc_clusters#outsider_account).
For technical reasons, your account must be expired before making the invitation request.
Once the invitation is processed, we will reactivate your account, and you will retain access to your data.

#### Connection to Cluster

**Q:** When I type my password, no characters are printed. Why?
**A:** Unlike Windows systems, Linux and Unix systems **do not display** any characters (not even *) when you enter your password in a terminal. The field remains blank, and the cursor will not move.

Simply type your password and press Enter. Your connection should be successful.

Please be cautious not to mistype your password multiple times, as you may be temporarily blocked (see below).

**Q:** When I tried to connect to the cluster, there is no response.
**A:** We employ `fail2ban` on the clusters to prevent brute-force attacks.

If you enter the wrong password three times consecutively, you will be banned for 15 minutes (`fail2ban` will blacklist your IP address). After 15 minutes, you can attempt to connect again.

If you are still unable to connect after 15 minutes, please contact us with the following information:

- Your username
- Your IP address (you can find it using [this web service](https://whatismyipaddress.com)).
- The cluster you are attempting to connect to.

**Q:** SSH "Could not resolve hostname XXXX: Name or service not known"
**A:** It means the specified hostname cannot be found, either due to a typo or because the DNS can't resolve it.
- check the login node [hostname](https://doc.eresearch.unige.ch/hpc/access_the_hpc_clusters#login_nodes)

PS: Keep in mind that baobab2 has been decommissioned for 2 years.

**Q:** When I try to connect to Clusters using `ssh` or `sftp`, I see the message: Connection refused
```
Connection refused

```
**A:** This may occur because you attempted to connect multiple times with incorrect credentials (e.g., wrong username or password), causing your IP address to be blacklisted. Your IP address will be automatically unblocked after 15 minutes.

Please note that your Baobab/Yggdrasil password is the same as your ISIS password, which we do not manage. If you forgot your password or need to verify it, please use the following service: [mdp.unige.ch](https://mdp.unige.ch).

**Q:** I have forgotten my password, can the HPC team reset it?
**A:** No, your Baobab/Yggdrasil/Bamboo password is your ISIS password, and we do not manage it.

If you **forgot your password** or need to verify it, please use the following service:
- [mdp.unige.ch](https://mdp.unige.ch).

**Q:** How to check my SshPublicKey ?
!!!
- If you are a **collaborator/student/external** user Check on [my-account](https://my-account.unige.ch/main/home)
- If you are an **Outsider** user Check on [applicant](https://applicant.unige.ch/main/outsider-info/update/)

For more informations please refer to [ssh PublicKey](https://doc.eresearch.unige.ch/hpc/access_the_hpc_clusters#ssh_publickey) page.

**Q:** Is it possible to access my account from more than one SSH key?
**A:** Yes, please check [Access the clusters](https://doc.eresearch.unige.ch/hpc/access_the_hpc_clusters#multiple_ssh_key)

**Q:** I tried to connect without success.
**A:** There are three possible reasons why you may not be able to connect:
- **The cluster is under maintenance.** Maintenance occurs periodically. Please check your email (including junk/spam folders) or visit the [HPC-community](https://hpc-community.unige.ch/) for announcements.

- **Your network is blocking access to our clusters or the SSH protocol.** We use public IP addresses for the login nodes. If you cannot connect, please contact your local network administrator to determine if there are any restrictions on accessing `login1.baobab.hpc.unige.ch`, `login1.yggdrasil.hpc.unige.ch`, or `login1.bamboo.hpc.unige.ch`, or if port 22 is blocked. you can receive this message : `ssh: connect to username@login1.baobab.hpc.unige.ch port 22: Connection timed out`

- **The login node is down.** While unlikely, if this occurs, please wait a little or contact us if the issue persists beyond 15 minutes.

**Q:** I'm having trouble connecting with software like PuTTY, FileZilla, or X2Go. What should I do?
**A:** The most often reason is your is not up to date. Check update and try again.
The most common reason is that your software is outdated. Please check for updates and try again.
If the issue persists, refer to the FAQ section on [https://doc.eresearch.unige.ch/hpc/faq#connection_to_cluster](https://doc.eresearch.unige.ch/hpc/faq#connection_to_cluster) for more troubleshooting steps.

**Q:** I'm an "outsider" and I can't connect to open on demand
**A:** This should be fixed in the future, but in the meantime, we have a workaround:
- connect to https://openondemand.baobab.hpc.unige.ch using the same email you used when register as outsider. You'll get the following error: `Error -- failed to map user ()`
- Then go to [the session page](https://openondemand.baobab.hpc.unige.ch/Shibboleth.sso/Session) and send us a screenshot. We'll activate manually your account.

#### X2GO-Desktop
**Q:** Why I can't connect with x2go ?
**A:** We have already identified a number of common problems:

- Check the general FAQ:  [connection_to_cluster](https://doc.eresearch.unige.ch/hpc/faq/#connection_to_cluster)
- [Check your quota](https://doc.eresearch.unige.ch/hpc/storage_on_hpc#check_disk_usage_on_the_clusters); reaching the limit will prevent you from writing to your directory, which means X2Go won’t be able to initialize the necessary configurations.
- If you're using Anaconda/conda, try commenting out the conda block in your .bashrc file.
```
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/path/to/your/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/path/to/your/anaconda3/etc/profile.d/conda.sh" ]; then
        . "/path/to/your/anaconda3/etc/profile.d/conda.sh"
    else
        export PATH="/path/to/your/anaconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

```
- Make a backup(**steps by steps**) of the folowing files or directories and try to login again:
  - **~/.bashrc**
  - **~/.Xauthority**
  - **~/.x2go**
  - **~/.local/session**
  - **~/.config/xfce**
#### Storage

**Q:** I have a question about the storage !?
- Where should I store my files?
- What should I do if I delete something by mistake?
- Is there a backup?
- How can I restore a deleted file?
- How much storage space is available?
- My job creates lots of temporary small files, and everything is slow. What should I do?

**A:** For detailed information on all storage-related topics, please refer to our [Storage page](https://doc.eresearch.unige.ch/hpc/storage_on_hpc). This page provides comprehensive guidance on file storage, recovery, and managing storage space efficiently.

If you need to store a large amount of data, consider using the "Academic NAS" service, which you can find here: Academic NAS.

**Q:** How can I request for a shared directory
!!!A shared directory is a directory with access permissions granted to a specific group for sharing data.
Since June 2025, all new shared directories must use groups declared in the Active Directory of UNIGE.

To request a shared directory please fill the [HPC form](https://dw.unige.ch/openentry.html?tid=hpc) on DW.
(For Oustsider please refer to your PI/repondant)

**Q:** How can I access to a shared directory?
**A:** To access a **shared directory**, you need to be added to the appropriate group.

For shared directories using a group from Active Directory (groups beginning with GG or GL - example: GL_S_SCIENCES_POSY_LET, please contact your CI (Correspondant Informatique) or use the [ADaccess form](https://dw.unige.ch/openentry.html?tid=adaccess) on DW.

For old group (share_XXX, private_XXX), please send an email to [https://doc.eresearch.unige.ch/hpc@unige.ch](https://doc.eresearch.unige.ch/hpc@unige.ch) including relevant information (Uusername, Group, private_partion etc...) with the responsible person for the share or partition in CC. The responsible person **must** approve the modification.

**Q:** How can I copy data from one cluster to another one?
**A:** If you have a lot of data, the best way is to use rsync between both clusers, so you won't have to copy the data to your laptop first. [Transfer data from one cluster to another](https://doc.eresearch.unige.ch/hpc/best_practices#transfer_data_from_one_cluster_to_another)

#### Applications

**Q:** What applications are installed on Clusters ?

**A:** You can find information about available applications [here](https://doc.eresearch.unige.ch/hpc/applications_and_libraries#find_installed_applications_with_module)

**Q:** The software I need is not available on Clusters: what should I do ?
**A:** Please check [this documentation](https://doc.eresearch.unige.ch/hpc/applications_and_libraries#what_do_i_do_when_an_application_is_not_available_via_module).

**Q:** Can I use any Microsoft Windows software ?
**A:** Baobab is a GNU/Linux only machine, like the majority of academic clusters. If you have a windows software that could run on a Windows cluster, contact us at [https://doc.eresearch.unige.ch/hpc@unige.ch](https://doc.eresearch.unige.ch/hpc@unige.ch), perhaps we could find some solutions.

**Q:** Can I use a proprietary licensed software ?
**A:** Yes we can install it, but you should pay the required license. Send us a request at [https://doc.eresearch.unige.ch/hpc@unige.ch](https://doc.eresearch.unige.ch/hpc@unige.ch).

**Q:** I need a different Linux distributions/version, am I stuck ?

**A:** No, please check the [Apptainer](https://doc.eresearch.unige.ch/hpc/applications_and_libraries#apptainer_was_singularity) documentation.

**Q:** Illegal instruction

**A:** If you run a program and it crashes with an error `"Illegal instruction"` the reason is probably because
you have compiled your program on Baobab login node and your program is running on an older compute node
on which the CPU lacks some specialized functionality that were used during the compilation.

You have two possibilities:
- Recompile your program with less optimization, or compile on an older node. See [Advanced users](https://doc.eresearch.unige.ch/hpc/hpc_clusters#for_advanced_users)
- Only run your program on newer servers. See [Specify the CPU type you want](https://doc.eresearch.unige.ch/hpc/slurm#specify_the_cpu_type_you_want) and [Compute nodes](https://doc.eresearch.unige.ch/hpc/hpc_clusters#compute_nodes).

**Q:** How can I use another Python version ?
**A:** You need to distinguish between the system-installed Python package and the Python versions provided by `module` or `easybuild`. Since we support a variety of software needs for our users, we use module to manage different software versions, including multiple Python versions. To switch between them, you can use the module command to load the specific Python version you need.

```
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  Python:
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Description:
      Python is a programming language that lets you work more quickly and integrate your systems more effectively.

     Versions:
        Python/2.7.11
        [...]
        Python/3.11.5

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  For detailed information about a specific "Python" package (including how to load the modules) use the module's full name.
  Note that names that have a trailing (E) are extensions provided by other modules.
  For example:

     $ module spider Python/3.11.5
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

```

**Q:** Can I load two versions of the same software? How can I use two different software versions with different GCC dependencies?
**A:** No, you cannot load two versions of the same software simultaneously. Additionally, if two software packages depend on different GCC versions, you will not be able to load them at the same time.

In this case you need to check if there is another version available compatible with the toolchain (`GCC`, `foss` etc...) you want to use. If not, please refer to [The software I need is not available on Clusters: what should I do ?](https://doc.eresearch.unige.ch/hpc/faq#the_software_i_need_is_not_ava).
#### Slurm: job scheduler
**Q:** What is Slurm ?
**A:** Slurm is a job scheduling system used to manage and allocate resources in a computing cluster. It helps you submit, monitor, and control jobs (tasks) on the cluster.
Please take a moment to review this very important section: [Slurm and job management](https://doc.eresearch.unige.ch/hpc/slurm)

**As a reminder**: It is **forbidden** to run heavy compute jobs on the login nodes, you **must** use a compute node instead.
**Q:** I am already familiar with `torque/pbs/sge/lsf/...`, what are the equivalent concepts in slurm ?
**A:** Have a look at this scheduler "rosetta stone", available here:  http://slurm.schedmd.com/rosetta.pdf

**Q:** Can I run some small test runs in the login node ?
**A:** **No never**. You **must** use SLURM to run any test. The `debug` partition is dedicated to small tests.

**Q:** What partition should I choose ?
See our documentation about [Slurm Partitions](https://doc.eresearch.unige.ch/hpc/slurm#which_partition_for_my_job).

**Q:** Can I launch a job longer than 4 days ?
**A:** No. Unfortunately you can't. If we raised this limit, you will have to wait longer before having your pending jobs started. We think that the 4 days limit is a good trade-off.    However there could be two work-around if you experience an issue with this limit:

- Some software feature **checkpointing**. During runtime, the program will periodically save its current state on the disks. In that case, this snapshot may be used to resume the computation by another job.  Check if your program allows checkpointing. If you cannot find the information, try contacting the developer or ask us at [https://doc.eresearch.unige.ch/hpc@unige.ch](https://doc.eresearch.unige.ch/hpc@unige.ch).
- You could add private notes to Baobab. In that case the limit will be raised to 7 days or even higher. If you are interested, contact us.

**Q:** How are the priorities computed ?
**A:** See [here](https://doc.eresearch.unige.ch/hpc/slurm#how_is_the_priority_of_a_job_determined)

To get the priority calculation details of the jobs in the pending queue, you can use the command: `sprio -w`. You can also have a look at the weights, by typing `sprio -l`.

**Q:** Why My jobs stay a long time in the pending queue ?

**A:** See
- [Which partition for my job](https://doc.eresearch.unige.ch/hpc/slurm#which_partition_for_my_job)
- [Job priorities](https://doc.eresearch.unige.ch/hpc/slurm#job_priorities)
- [Stop wasting resources](https://doc.eresearch.unige.ch/hpc/best_practices#stop_wasting_resources)

**Q:** Can I run interactive tasks ?

!!!Yes, you can. But it is really awkward because you cannot be sure when your job will start.

See [Interactive jobs](https://doc.eresearch.unige.ch/hpc/slurm#interactive_jobs)

You may be interesting about [OpenOnDemand](https://doc.eresearch.unige.ch/hpc/how_to_use_openondemand) which provide a graphical to start Interactive session ( JupyterLab, MatLab, VScode, R etc...)

**Q:** I want to run several time the same job with different parameters
!!!In that case you can use the **job arrays** feature of SLURM. Please, have a look at the documentation [Job array](https://doc.eresearch.unige.ch/hpc/slurm#job_array)

**Q:** Is the nodes from my partition in use?
**A:** You can't just check with the partition name as the node may be in use by another partition such as shared-cpu. Here is an example to check the usage of your comput nodes:
```
squeue -w $(sinfo --noheader --partition <to be replaced by partition-name> --format="%n" | nodeset -f)

```

**Q:** Why I'm not able to use all the cores of a compute node ?
!!!Indeed, we are reserving two cores per node for system tasks such as data transfer, and os stuff.

```
(yggdrasil)-[root@admin1 ~]$ scontrol show node cpu001
NodeName=cpu001 Arch=x86_64 CoresPerSocket=18
   CPUAlloc=0 CPUEfctv=34 CPUTot=36 CPULoad=0.01
   AvailableFeatures=GOLD-6240,XEON_GOLD_6240,V9
   ActiveFeatures=GOLD-6240,XEON_GOLD_6240,V9
   Gres=(null)
   NodeAddr=cpu001 NodeHostName=cpu001 Version=23.02.1
   OS=Linux 4.18.0-477.10.1.el8_8.x86_64 #1 SMP Tue May 16 11:38:37 UTC 2023
   RealMemory=187000 AllocMem=0 FreeMem=185338 Sockets=2 Boards=1
   CoreSpecCount=2 CPUSpecList=17,35 <==================== this means we have two specialization cores <<<<
   State=IDLE ThreadsPerCore=1 TmpDisk=150000 Weight=10 Owner=N/A MCS_label=N/A
   Partitions=debug-cpu
   BootTime=2023-08-10T12:08:11 SlurmdStartTime=2023-08-10T12:09:00
   LastBusyTime=2023-08-11T10:06:42 ResumeAfterTime=None
   CfgTRES=cpu=34,mem=187000M,billing=34
   AllocTRES=
   CapWatts=n/a
   CurrentWatts=0 AveWatts=0
   ExtSensorsJoules=n/s ExtSensorsWatts=0 ExtSensorsTemp=n/s

```

If you really need to use all the cores of a compute node, you can override this parameter: `--core-spec=0`. This will implicitly lead to an exclusive allocation of the node.

ref: https://slurm.schedmd.com/core_spec.html

**Q:** How can I access to a private slurm partition?
**A:** To use a **private Slurm partition**, you need to be added to the appropriate group.

Please send an email to [https://doc.eresearch.unige.ch/hpc@unige.ch](https://doc.eresearch.unige.ch/hpc@unige.ch) including relevant information (Uusername, Group, private_partion etc...) with the responsible person for the share or partition in CC. The responsible person **must** approve the modification.

#### Issues

**Q:** I have a keyboard issue using a Mac.
**A:** Please refer to this [keymap-issues-with-nx-from-mac-os-x](https://stackoverflow.com/questions/7018775/keymap-issues-with-nx-from-mac-os-x-lion-to-ubuntu/42094562#42094562) for a potential solution.

**Q:** When I ssh, I get the message : "cannot change locale (UTF-8): No such file or directory"
```console
-bash: warning: setlocale: LC_CTYPE: cannot change locale (UTF-8): No such file or directory

```
**A:** You can resolve this issue by following Step #1 [here](https://www.cyberciti.biz/faq/os-x-terminal-bash-warning-setlocale-lc_ctype-cannot-change-locale/).

Please ensure that you close all open terminals on your Mac and relaunch them.

**Q:** When I try to connect to the cluster from a Mac using `ssh -Y` and I receive an error like:
```
Can't connect to X11

```
**A:** This issue likely arises because Xorg is no longer provided natively on macOS. You need to install XQuartz.

Refer to this solution: [macOS High Sierra and X11 Forwarding](https://stackoverflow.com/questions/50035949/macos-high-sierra-and-x11-forwarding/50182736#50182736).

#### Switch edu-ID Login Issues

**Q:** I get an error message from Switch edu-ID while trying to access:
- https://hpc-community.unige.ch/
- https://openondemand.baobab.hpc.unige.ch/
- https://openondemand.bamboo.hpc.unige.ch/

**A:** Please follow these links for support:
- [SWITCH edu-ID Account - Welcome](https://plone.unige.ch/distic/pub/compte-switch-edu-id/compte-switch-edu-id-accueil#EN)
- [How to Create or Verify Your SWITCH edu-ID Account and Link It to UNIGE](https://plone.unige.ch/distic/pub/compte-switch-edu-id/comment-creer-compte-switch-edu-id#EN)

Ensure that you are using the email address linked to your Switch edu-ID account.

Please also note that your ISIS (UNIGE) password and your Switch edu-ID password are not the same. Verify that you are using the correct password when logging in.

#### HPC community forum

**Q:** I don't find a way to receive email summary of new post

**A:** You can activate the email summary following those steps:
