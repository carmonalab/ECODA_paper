# Issue on /home on bamboo

- Source: https://hpc-community.unige.ch/t/4070

- Created: 2025-08-28T07:06:14.626Z

- Tags: bamboo

- Posts: 12

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Lucille.Delisle1 (2025-08-28T07:06:14.680Z)

Hi,
It seems that there is an error with the /home on bamboo. When I ssh I get:
`Could not chdir to home directory /home/users/d/delislel: Communication error on send `
`-bash: /home/users/d/delislel/.bash_profile: Communication error on send `
`(bamboo)-[delislel@login1 /]$ `
If I try using on demand:
# Home directory not found
Your home directory appears to be missing. The home directory mount may be unavailable, or your home directory may need to still be created. Please contact support for help and attempt to restart your web server by clicking below when the problem has been fixed.
Thanks


## Post 2 by @Lucille.Delisle1 (2025-08-28T07:10:24.766Z)

This is fixed now, I don’t know if you solved it or if it was a temporary issue.


## Post 3 by @nicolas.clairis (2025-08-28T07:33:19.850Z)

thanks for the post, I was just wondering why all my batches crashed around 9am :sweat_smile:


## Post 4 by @nicolas.clairis (2025-08-28T07:35:42.497Z)

Although the access to /home seems to be working, I still get an error message when trying to access the scratch via FileZilla though
`Commande :	cd “/home/users/c/clairis/scratch”`
`Erreur :	Directory /home/users/c/clairis/scratch: received failure with description ‘Failure’`
`Erreur :	Impossible de récupérer le contenu du dossier`
and same when trying via Putty:
`cd scratch`
`-bash: cd: scratch: Communication error on send`


## Post 5 by @Lucille.Delisle1 (2025-08-28T07:37:29.002Z)

Indeed, I have the same problem in command lines:
```
(bamboo)-[delislel@login1 ~]$ ls scratch
ls: cannot access 'scratch': Communication error on send 
(bamboo)-[delislel@login1 ~]$ ls /srv/beegfs/scratch/users
ls: cannot access '/srv/beegfs/scratch/users': Communication error on send 
```


## Post 6 by @Lucille.Delisle1 (2025-08-28T07:49:55.313Z)

This is now solved for me.


## Post 7 by @nicolas.clairis (2025-08-28T08:02:22.706Z)

thanks for warning! looks ok here as well indeed, weirdly they didn’t say anything yet so I hope it’s stable


## Post 8 by @Adrien.Albert (2025-08-28T08:34:38.544Z)

Dear Users,
You can follow all issue in the post [2025] Current issues on HPC Cluster - #21 by Yann.Sagon[[2025] Current issues on HPC Cluster - #21 by Yann.Sagon](https://hpc-community.unige.ch/t/2025-current-issues-on-hpc-cluster/3788/21)
Best regards


## Post 9 by @nicolas.clairis (2025-08-28T08:36:17.604Z)

yes I know and I was watching, but it was marked as solved while we still had the issue with the access to scratch :sweat_smile:  Is it solved now?


## Post 10 by @Adrien.Albert (2025-08-28T08:37:37.777Z)

Yes it’s should work :wink:


## Post 11 by @Yann.Sagon (2025-08-28T08:53:12.902Z)

nicolas.clairis:
> it was marked as solved while we still had the issue
sorry about that, we figured out after the reboot of the server that a service didn’t started as expected, but yes, I confirm it is now running as expected. We contacted the vendor to have some insight to try to fix the root cause.


## Post 12 by @nicolas.clairis (2025-08-28T09:54:53.494Z)

ok great thanks for the update and confirmation that it works! :slight_smile: fingers crossed for no more crashing then
