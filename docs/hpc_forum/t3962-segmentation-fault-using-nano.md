# Segmentation fault using nano

- Source: https://hpc-community.unige.ch/t/3962

- Created: 2025-05-30T07:27:58.499Z

- Tags: yggdrasil

- Posts: 4

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Charles.Pouchon (2025-05-30T07:27:58.549Z)

Dear @Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon),
I installed Aster using Conda this morning, and now I’m encountering an issue with `nano`.
```
(base) (yggdrasil)-[pouchon@login1 ~]$ nano  
Segmentation fault (core dumped)
```
I can’t use it anymore. Do you have any idea what might be causing this?
Thanks in advance,
Best regards,
Charles


## Post 2 by @Yann.Sagon (2025-05-30T08:53:31.053Z)

Hello,
you can see that since you installed/used Conda, ano is provided by Conda:
```
(base) (yggdrasil)-[pouchon@login1 ~]$ which nano
~/miniforge3/bin/nano
```
Indeed, this appears in login1 logs:
```
[Fri May 30 08:29:18 2025] traps: nano[3312930] general protection fault ip:14f27ed59a3f sp:7fffd97387b8 error:0 in libncursesw.so.6.5[14f27ed45000+25000]
[Fri May 30 08:29:29 2025] traps: nano[3313050] general protection fault ip:15414c952a3f sp:7ffdaee4d928 error:0 in libncursesw.so.6.5[15414c93e000+25000]
```
But anyway, I’m able to launch this version of nano on login1.yggdrasil as user pouchon without issue, so I don’t understand what is different in your case. Did you logged out and login again?
If you want to disable conda, you can edit your `.bashrc` file and remove all the stuff that conda wrote to it.
```
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/home/users/p/pouchon/miniforge3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/users/p/pouchon/miniforge3/etc/profile.d/conda.sh" ]; then
        . "/home/users/p/pouchon/miniforge3/etc/profile.d/conda.sh"
    else
        export PATH="/home/users/p/pouchon/miniforge3/bin:$PATH"
    fi
fi
unset __conda_setup

if [ -f "/home/users/p/pouchon/miniforge3/etc/profile.d/mamba.sh" ]; then
    . "/home/users/p/pouchon/miniforge3/etc/profile.d/mamba.sh"
fi
# <<< conda initialize <<<
```


## Post 3 by @Charles.Pouchon (2025-05-30T09:18:11.770Z)

Hi,
Thanks a lot for your reply.
I installed `nano` with Conda just after my previous message, and it’s now working correctly, as you mentioned:
```
(base) (yggdrasil)-[pouchon@login1 ~]$ which nano  
~/miniforge3/bin/nano
```
Previously, I was getting the error with the `nano` installed here:
```
(base) (yggdrasil)-[pouchon@login1 ~]$ which nano  
/usr/bin/nano
```
So, I guess I can safely use the new `nano` now.
Thanks again,
Charles


## Post 4 by @Yann.Sagon (2025-05-30T10:41:15.257Z)

Ok I see. `nano` needs `libncursesw.so.6` which is provided by `conda` in your case and isn’t probably  compatible with legacy `nano` from the OS due to version change.
