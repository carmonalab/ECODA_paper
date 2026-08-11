# Issue with SSH key for Baobab

- Source: https://hpc-community.unige.ch/t/4337

- Created: 2026-07-07T18:32:48.415Z

- Tags: baobab

- Posts: 5

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Lilou.Dehondt (2026-07-07T18:32:48.488Z)

## Primary informations
Username: dehondt4
Cluster: Baobab
## Description
Hi, I already raised this issue in the comment of my previous post but never got an answer. I created a SSH key according to the tutorial shared in the forum, but still can’t access baobab. I think the issue is that my HPC username is “dehondt4”, linked to my student account (@etu.unige.ch). However, when I log in to my-account.unige.ch to register my SSH public key, it automatically registers under my @unige.ch account instead of my @etu.unige.ch account. (I have both accounts because I am a student but also work as a RA). Could you please help me with this issue?
I need to access Baobab quickly to be able to do some work, which I have not been able to make progress on in the last week without the help of the cluster. Thank you very much for your help and have a very nice day!
Lilou


## Post 2 by @Yann.Sagon (2026-07-08T14:18:06.684Z)

Dear @Lilou.Dehondt[@Lilou.Dehondt](https://hpc-community.unige.ch/u/lilou.dehondt)
Can you please check this post [tutorial] 🔐 SSH login troubleshooting (key mismatch)[[tutorial] 🔐 SSH login troubleshooting (key mismatch)](https://hpc-community.unige.ch/t/tutorial-ssh-login-troubleshooting-key-mismatch/4339)
And update this post if there is an issue?


## Post 3 by @Yann.Sagon (2026-07-08T14:42:57.156Z)

Hi again, I’ve checked: we don’t see any of your ssh key with your profile. Are you sure you uploaded a key in my-account?


## Post 4 by @Lilou.Dehondt (2026-07-08T15:29:42.871Z)

Hi,
I believe it’s because it adds automatically to dehondt and not dehondt4. On the “my account” website I can’t choose.
As you can see below, the my account website maps to both my student account, that I use with baobab, and my work acount. At the bottom the button “add a ssh key” does no allow me to pick between the two of the them. The SSH key ends with my @unige.ch and not my @etu.unige.ch email.
Thank you for taking time to look at this issue,
Best,
Lilou Dehondt


## Post 5 by @Yann.Sagon (2026-07-09T07:46:53.150Z)

Dear @Lilou.Dehondt[@Lilou.Dehondt](https://hpc-community.unige.ch/u/lilou.dehondt)
I think I incorrectly checked that your profile had an SSH key associated with it. I checked again this morning and found that there is one SSH key associated with both of your accounts which is correct. With two profiles on one account, as you have, you can’t choose which profile the SSH key applies to: it is associated with both profiles.
Then please do the checks [tutorial] 🔐 SSH login troubleshooting (key mismatch)[[tutorial] 🔐 SSH login troubleshooting (key mismatch)](https://hpc-community.unige.ch/t/tutorial-ssh-login-troubleshooting-key-mismatch/4339)
