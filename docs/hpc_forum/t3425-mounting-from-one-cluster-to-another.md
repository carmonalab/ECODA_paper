# Mounting from one cluster to another

- Source: https://hpc-community.unige.ch/t/3425

- Created: 2024-04-23T09:04:55.549Z

- Posts: 7

- Category: 8

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Van-Khoa.Nguyen (2024-04-23T09:04:55.603Z)

Dear HPC team,
I often switch between two clusters to be able to use available resources if one of the clusters is being used extensively.
I created a common folder on Yggdrasil to store files and mounted this folder to Baobab so that I do not need to move files between clusters frequently.
I successfully did the mounting with `sshfs`, and I can access all the files in the folder from the login node.
```
sshfs -o default_permissions user@login1.yggdrasil.hpc.unige.ch:path1 path2
```
However, on the compute nodes of Baobab, I can not see all files in the mounted folder.
Could you explain this problem, and would it be possible for me to access the mounted folder from the compute nodes ? Thanks !
Best regards,
Van Khoa


## Post 2 by @Adrien.Albert (2024-04-24T08:43:35.238Z)

Dear @Van-Khoa.Nguyen[@Van-Khoa.Nguyen](https://hpc-community.unige.ch/u/van-khoa.nguyen)
As you are only mounting the login, the compute node does not have access to this `sshfs` mount point.
However, `sshfs` is not installed on the compute node for a number of reasons:
- Bandwidth limitations, impacting data transfer efficiency.
- Latency dependency, impacting the responsiveness of computing tasks.
- Reliability issues, with potential disruptions to network connectivity.
- Impact on performance
Keep in mind that Yggdrasil is in Sauverny, Baobab in Uni dufour, and geographical distance has a significant impact on the network.
We are well aware of the need. The project for a common storage space for the clusters is currently under study, and we hope to find a solution compatible with the current IT infrastructure.
Best Regards


## Post 3 by @Gael.Rossignol (2024-04-24T11:53:13.533Z)

Dear Van-Khoa,
To complete my collegue answer, there are somme documentation about the sync fo data from both clusters using “rsync” :white_check_mark:
https://doc.eresearch.unige.ch/hpc/best_practices?s[]=rsync[https://doc.eresearch.unige.ch/hpc/best_practices?s[]=rsync](https://doc.eresearch.unige.ch/hpc/best_practices?s%5B%5D=rsync)
If you can provide more information regarding the path size and the frequency of the data changes, we can help you to find a way.
Best regards,


## Post 4 by @Van-Khoa.Nguyen (2024-04-26T07:37:42.885Z)

Dear @Adrien.Albert[@Adrien.Albert](https://hpc-community.unige.ch/u/Adrien.Albert) ,
I understand better the problem with `sshfs`, and thank you for your thoughtful explanation.
Best regards,
Van Khoa


## Post 5 by @Van-Khoa.Nguyen (2024-04-26T07:40:40.815Z)

Dear @Gael.Rossignol[@Gael.Rossignol](https://hpc-community.unige.ch/u/Gael.Rossignol),
Thank you for your answer and the proposed solution.
I will try with your instruction.
Best regards,
Van Khoa


## Post 6 by @Van-Khoa.Nguyen (2024-04-26T14:07:17.005Z)

Dear @Gael.Rossignol[@Gael.Rossignol](https://hpc-community.unige.ch/u/Gael.Rossignol)  and @Adrien.Albert[@Adrien.Albert](https://hpc-community.unige.ch/u/Adrien.Albert) ,
Can we log in from Yggdrasil to Baobab with passwordless authentication?
I tried to add Yggdrasil’s public key to Baobab by following this tutorial, but it seems not working.
thegeekstuff.com[thegeekstuff.com](https://www.thegeekstuff.com/2008/11/3-steps-to-perform-ssh-login-without-password-using-ssh-keygen-ssh-copy-id/)
### 3 Steps to Perform SSH Login Without Password Using ssh-keygen & ssh-copy-id[3 Steps to Perform SSH Login Without Password Using ssh-keygen & ssh-copy-id](https://www.thegeekstuff.com/2008/11/3-steps-to-perform-ssh-login-without-password-using-ssh-keygen-ssh-copy-id/)
You can login to a remote Linux server without entering password in 3 simple steps using ssky-keygen and ssh-copy-id as explained in this article. ssh-keygen creates the public and private keys. ssh-copy-id copies the local-host’s public key to the...
I try to get rid of the password request every time I use rsync. Thanks!
Best regards,
Van Khoa


## Post 7 by @Yann.Sagon (2024-04-29T09:56:13.894Z)

Hi, better use our doc :stuck_out_tongue_winking_eye: hpc:access_the_hpc_clusters [eResearch Doc][hpc:access_the_hpc_clusters [eResearch Doc]](https://doc.eresearch.unige.ch/hpc/access_the_hpc_clusters#ssh_key)
