# Storage usage script

- Source: https://hpc-community.unige.ch/t/4258

- Created: 2026-03-18T17:01:33.341Z

- Posts: 3

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Marco.Froelich (2026-03-18T17:01:33.398Z)

cluster: bamboo
user: froelicm
Hi!
I was logged out from Bamboo due to reaching disk quota. I was surprised, but was able to remove some files in my home and get back access promptly. To understand where the disk quota error was coming from, I wanted to inspect my allowances.
The documentation points to the following script to get the current personal disk usage.
beegfs-get-quota-home-scratch.sh
First of all, I get the following error:
```
home dir: /home/users/f/froelicm
scratch dir: /srv/beegfs/scratch/users/f/froelicm

storage     |   name        |  id  ||    used    |    hard    ||  used   |  hard----------------------------|------||------------|------------||---------|---------

Error: Failed to load map file: /etc/beegfs/home.d/beegfs-client.conf

[BeeGFS Control Tool Version: 7.4.7Refer to the default config file (/etc/beegfs/beegfs-client.conf)or visit http://www.beegfs.com to find out about configuration options.]

home        |

Error: Failed to load map file: /etc/beegfs/scratch.d/beegfs-client.conf
```
On the forum, this issue[issue](https://hpc-community.unige.ch//hpc-community.unige.ch/t/howto-check-disk-usage-on-baobab/860) was posted to check usage in a different way. I get the same error for both:
```
Error: Failed to load map file: /etc/beegfs/scratch.d/beegfs-client.conf

[BeeGFS Control Tool Version: 7.4.7
Refer to the default config file (/etc/beegfs/beegfs-client.conf)
or visit http://www.beegfs.com to find out about configuration options.]
```
In addition, I am working within a group lead by my PI. Is there a script to retrieve group-level disk usage (similarly to the script used to obtain CPUhours counts, which gives information about the group for all members of the group) ?
Thanks,
Marco


## Post 2 by @Adrien.Albert (2026-03-19T11:48:13.424Z)

Hello;
This has been fixed:
```
(bamboo)-[alberta@login1 ~]$ beegfs-get-quota-home-scratch.sh 
home dir: /home/users/a/alberta
scratch dir: /srv/beegfs/scratch/users/a/alberta

          user/group                 ||           size          ||    chunk files
  storage     |   name        |  id  ||    used    |    hard    ||  used   |  hard
  ----------------------------|------||------------|------------||---------|---------
home        |        alberta|401775||  185.41 GiB| 1024.00 GiB||   213852|unlimited
scratch     |        alberta|401775||      0 Byte|   unlimited||        0| 10000000
```
Let me know everything is working on your side.


## Post 3 by @Marco.Froelich (2026-03-19T14:50:11.654Z)

Adrien.Albert:
> `beegfs-get-quota-home-scratch.sh `
All three scripts working, merci Adrien !
