# Cannot connect to yggdrasil via X2Go, no error message

- Source: https://hpc-community.unige.ch/t/3847

- Created: 2025-03-04T11:27:50.237Z

- Tags: yggdrasil

- Posts: 8

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yuzheng.Kang (2025-03-04T11:27:50.281Z)

Hi SysAdmin,
Recently, I am trying to connect Yggdrasil via X2Go on Mac. It was working before. Turns out there is an empty error message. I could not debug it, could you help me to solve this issue? Here are the screenshot and settings I have.
Moreover, I tired to comment out the conda block in the .bashrc file. But it didn’t work.
Screenshot 2025-03-04 at 12.23.15
Screenshot 2025-03-04 at 12.23.15908×1018 77.3 KB
[Screenshot 2025-03-04 at 12.23.15908×1018 77.3 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/7b10f6d3b3c500a406bcc129de479658cab202c6.jpeg)
Screenshot 2025-03-04 at 12.23.42
Screenshot 2025-03-04 at 12.23.421297×860 98.3 KB
[Screenshot 2025-03-04 at 12.23.421297×860 98.3 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/385f8f4044c552c990964e3b047fc0867db70670.jpeg)


## Post 2 by @Gael.Rossignol (2025-03-06T08:54:47.624Z)

Dear Yuzheng,
Could you please try to desactivate your personal .bashrc and try again?
Best regard,


## Post 3 by @Yuzheng.Kang (2025-03-06T10:17:06.164Z)

Dear Gael,
Unfortunately, it didn’t work.


## Post 4 by @Gael.Rossignol (2025-03-06T15:30:19.859Z)

Gael.Rossignol:
> Dear Yuzheng,
> Could you please try to desactivate your personal .bashrc and try again?
> Best regard
Dear Yuzheng,
Could you please try again, I fixed something.
Best regard


## Post 5 by @Yuzheng.Kang (2025-03-06T15:36:58.275Z)

Hi,
It still not working. Just to make sure, I removed .bashrc file, is that you want me to try? I did it again, but not working.
Best,
Yuzheng


## Post 6 by @Gael.Rossignol (2025-03-06T15:47:06.407Z)

Yuzheng.Kang:
> gain, but not working
And you don’t have any error message? Just a question do you have updated yoiur computer to latest mac os?
Maybe reinstall x2go on your personal computer?
Thanks,


## Post 7 by @Yuzheng.Kang (2025-03-07T09:20:31.840Z)

Dear Gael,
Thanks for your help. In the end, problem solved by re-installing the x2go. So I guess the software was corrupted.
Best,
Yuzheng


## Post 8 by @Gael.Rossignol (2025-03-10T08:38:53.437Z)

Hello,
I’m happy you solved you problems.
Have a nice day,
