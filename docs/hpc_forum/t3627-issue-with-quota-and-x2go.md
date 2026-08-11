# Issue with quota and x2go

- Source: https://hpc-community.unige.ch/t/3627

- Created: 2024-09-03T09:32:00.882Z

- Tags: yggdrasil

- Posts: 9

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Elliot.Jiwani-Brown (2024-09-03T09:32:00.882Z)

I appear to be having X Authorisation issues with the yggdrasil login via X2Go (also affects the ssh when using -X). This occurred after trying to modify my session preferences as my connection was noticeably lagging when in use. Since then I have been able to login using X2Go or run GUI services from an ssh login.
If you have any help or advice on this (perhaps just a refresh of my login is the solution), this would be greatly appreciated.
Thanks
Elliot


## Post 2 by @daniel.forerosanchez (2024-09-04T06:47:54.348Z)

It seems Yggdrasil is down yet again.
Freezes at
```
(base) daniel@daniel-nb:~$ ssh -vv dforeros@login1.yggdrasil.hpc.unige.ch
OpenSSH_8.9p1 Ubuntu-3ubuntu0.10, OpenSSL 3.0.2 15 Mar 2022
debug1: Reading configuration data /home/daniel/.ssh/config
debug1: /home/daniel/.ssh/config line 45: Applying options for login1.yggdrasil.hpc.unige.ch
debug1: Reading configuration data /etc/ssh/ssh_config
debug1: /etc/ssh/ssh_config line 19: include /etc/ssh/ssh_config.d/*.conf matched no files
debug1: /etc/ssh/ssh_config line 21: Applying options for *
debug2: resolving "login1.yggdrasil.hpc.unige.ch" port 22
debug1: Connecting to login1.yggdrasil.hpc.unige.ch [129.194.64.11] port 22.
```
UPDATE: It was an issue of eduroam apparently, switching wifi worked


## Post 3 by @Adrien.Albert (2024-09-04T10:51:26.196Z)

Hi,
@daniel.forerosanchez[@daniel.forerosanchez](https://hpc-community.unige.ch/u/daniel.forerosanchez)
I can connect without any issue.
@Elliot.Jiwani-Brown[@Elliot.Jiwani-Brown](https://hpc-community.unige.ch/u/elliot.jiwani-brown)
Could you take a look to this FAQ:
https://doc.eresearch.unige.ch/hpc/faq#x2go-desktop[https://doc.eresearch.unige.ch/hpc/faq#x2go-desktop](https://doc.eresearch.unige.ch/hpc/faq#x2go-desktop)
If you still get an error, could you share it?


## Post 4 by @Elliot.Jiwani-Brown (2024-09-04T12:01:39.990Z)

I looked at the link and found that when checking the quota, apparently I have reached my limit, however when I view the contents of my home (jiwanibr) directory using < du -hsc * | sort -hr > it shows that I have used only 30GB of my quota. Is there a trash/deleted folder somewhere that won’t show up using this command and may have the hidden files?
Also, in the instructions on the link you provide, it is not clear what the following files need to be renamed as. Or do I need to create these folders if they are not already made?


## Post 5 by @Adrien.Albert (2024-09-04T13:28:49.828Z)

Hi @Elliot.Jiwani-Brown[@Elliot.Jiwani-Brown](https://hpc-community.unige.ch/u/elliot.jiwani-brown)
Indeed you have reach your quota on Home filesystem.
```
(yggdrasil)-[root@admin1 ~]$ beegfs-get-quota-home-scratch.sh -u jiwanibr 
          user/group                 ||           size          ||    chunk files
  storage     |   name        |  id  ||    used    |    hard    ||  used   |  hard
  ----------------------------|------||------------|------------||---------|---------
home        |       jiwanibr  |337317||    1.02 TiB| 1024.00 GiB||   343047|unlimited
scratch     |       jiwanibr  |337317||    7.21 TiB|   unlimited||   453546| 10000000
```
Note that the quota includes the entire `/home` directory, and therefore the `/home/share/cdff` directory in which you have data.


## Post 6 by @Adrien.Albert (2024-09-04T13:38:40.615Z)

Elliot.Jiwani-Brown:
> Also, in the instructions on the link you provide, it is not clear what the following files need to be renamed as. Or do I need to create these folders if they are not already made?
I understand what you mean, in this context rename means “make a backup by renaming”. But in your case (quota issue) you don’t need to follow these steps.
image
Exemple: Rename `~.bashrc` in `~.bashrc.orig`
```
mv ~/.bashrc{,.orig}
```
I’ve updated the FAQ, let me know if it’s clearer.


## Post 7 by @Elliot.Jiwani-Brown (2024-09-04T14:38:04.178Z)

Interesting. This share folder is the share folder of our research group and should be based in the scratch folder.
Do you know why it is linked to my home directory, and is there a simple solution to unlink it and keep the same data in the scratch folder?
(update) This also may have been originally been made thinking it was connected to the scratch folder mistakenly. I will copy it to my scratch folder and hopefully this fixes the issue.
Thank you for the response and updates


## Post 8 by @Elliot.Jiwani-Brown (2024-09-05T10:48:22.965Z)

Unfortunately this has not fixed the issue. The memory of my home folder (/home/users/j/jiwanibr/) only appears to contain 39GB of data, however still shows I have 1.02TB of data stored somewhere.
I moved the data (using mv) from /home/share/cdff to the srv/beegfs/scratch/shares folder however it appears to have not removed itself from the original location.
Could this be that the /home/share/cdff folder is linked somehow to my home login? I have a shortcut to it on my Desktop so perhaps this is the issue? Otherwise I have no idea where else to look for this 960GB of data that does not appear in my home folder


## Post 9 by @Adrien.Albert (2024-09-10T07:57:28.607Z)

Hi @Elliot.Jiwani-Brown[@Elliot.Jiwani-Brown](https://hpc-community.unige.ch/u/elliot.jiwani-brown)
The symbolic link in your home directory pointing to your scratch directory has no impact on your Home quota.
You have other data in the sub-directory `/home/share/cdff/...`
I can see you have 1TB on home filesystem: `/home/share/cdff/indocat`
```
ncdu 1.20 ~ Use the arrow keys to navigate, press ? for help                                                                                                                                                       
--- /home/share/cdff/indocat --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  850.7 GiB [##############################] /mseed_SDS                                                                                                                                                            
  160.1 GiB [#####                         ] /CCF_ML165_rma_corr
    1.8 MiB [                              ]  ML165_stations_nodes.xml

*Total disk usage:   1.0 TiB   Apparent size:   1.0 TiB   Items: 11584  
```
