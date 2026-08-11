# .conda ghost files

- Source: https://hpc-community.unige.ch/t/3495

- Created: 2024-06-17T17:29:11.927Z

- Posts: 3

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Max.Briel (2024-06-17T17:29:11.976Z)

Hi,
cluster: yggdrasil
user: briel
I’m experiencing an issue with installing and removing conda environments.
I ran into issues with conda being unable to remove its files using `conda clean --all` and being unable to install some new packages in a new environments (such as jupyter) due to files being unable to be removed. I suspect my .conda folder got corrupted during installation or removal of part of an environment.
As such, I tried resetting my conda environment. However, there are currently “ghost” files in my .conda folder that cannot be deleted by me. rm -rf does not work.
Can these please be deleted (or the complete .conda folder in my home directory)?
Or could it be reverted back to an earlier working date (ie yesterday)?
example:
`/home/users/b/briel/.conda/pkgs/webencodings-0.5.1-py311h06a4308_1/lib/python3.11/site-packages/webencodings-0.5.1-py3.11.egg-info/PKG-INFO`
`ls` returns:
```
ls: cannot access 'PKG-INFO': No such file or directory
ls: cannot access 'SOURCES.txt': No such file or directory
ls: cannot access 'dependency_links.txt': No such file or directory
ls: cannot access 'top_level.txt': No such file or directory
total 0
-????????? ? ? ? ?            ? dependency_links.txt
-????????? ? ? ? ?            ? PKG-INFO
-????????? ? ? ? ?            ? SOURCES.txt
-????????? ? ? ? ?            ? top_level.txt
```


## Post 2 by @Adrien.Albert (2024-06-21T15:11:07.994Z)

Hi @Max.Briel[@Max.Briel](https://hpc-community.unige.ch/u/max.briel)
Indeed, it seems that some metadata was out of sync, I forced the file system to check itself :
```
(yggdrasil)-[root@login1 ~]$ ls /home/users/b/briel/.conda/pkgs/webencodings-0.5.1-py311h06a4308_1/lib/python3.11/site-packages/webencodings-0.5.1-py3.11.egg-info/ -l
total 0
(yggdrasil)-[root@login1 ~]$ ls /home/users/b/briel/.conda/pkgs/webencodings-0.5.1-py311h06a4308_1/lib/python3.11/site-packages/ -l
total 1
drwxr-sr-x 3 briel hpc_users 1 Jun 21 02:26 webencodings
drwxr-sr-x 2 briel hpc_users 0 Jun 21 02:26 webencodings-0.5.1-py3.11.egg-info
```
I think it’s all good, could you confirm it ?


## Post 3 by @Max.Briel (2024-06-27T12:20:35.153Z)

Thank you! This is now fixed! :smiley:
