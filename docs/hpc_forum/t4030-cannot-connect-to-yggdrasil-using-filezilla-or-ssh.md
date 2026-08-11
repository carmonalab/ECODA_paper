# Cannot connect to yggdrasil using FileZilla, or ssh

- Source: https://hpc-community.unige.ch/t/4030

- Created: 2025-08-04T12:43:28.210Z

- Posts: 5

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @William.Ceva (2025-08-04T12:43:28.250Z)

Hello,
For some reason, I can only connect to yggdrasil via x2go at the moment.
I need to transfer some files off of yggdrasil and onto my own computer, but I cannot connect to yggdrasil at all using FileZilla.  I also cannot connect to yggdrasil via standard ssh…only x2go works…?
Thank you for any help you can provide,
Will Ceva


## Post 2 by @Adrien.Albert (2025-08-04T12:53:11.761Z)

Hi @William.Ceva[@William.Ceva](https://hpc-community.unige.ch/u/william.ceva)
Please post on HPC Support > HPC issues[HPC Support > HPC issues](https://hpc-community.unige.ch/c/hpc-support/hpc-issues/14) to use the template post.
Could you give me more information about your issue, do you have any errors messages ?


## Post 3 by @William.Ceva (2025-08-04T13:25:57.832Z)

Hello, here are the error messages I receive from FileZilla when I attempt to connect
Screenshot 2025-08-04 092524
Screenshot 2025-08-04 0925241912×566 67.8 KB
[Screenshot 2025-08-04 0925241912×566 67.8 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/439b88b9128d5967f9b568071f279056f187a5e0.png)


## Post 4 by @Adrien.Albert (2025-08-04T14:41:36.969Z)

Hi @William.Ceva[@William.Ceva](https://hpc-community.unige.ch/u/william.ceva)
I comment you conda block in your bashrc, could you try again ?
Let me know if it’s working ?
Also you may be interesting using conda in container:
doc.eresearch.unige.ch[doc.eresearch.unige.ch](https://doc.eresearch.unige.ch/hpc/applications_and_libraries#how_to_create_a_conda_environment_in_a_container)
### hpc:applications_and_libraries [eResearch Doc][hpc:applications_and_libraries [eResearch Doc]](https://doc.eresearch.unige.ch/hpc/applications_and_libraries#how_to_create_a_conda_environment_in_a_container)


## Post 5 by @William.Ceva (2025-08-04T15:54:26.244Z)

Hi Adrien,
Yes, I was able to establish an sftp connection with FileZilla and grab the files I needed.  Weird though, that the conda block in my bashrc was the issue (I have connected via FileZilla numerous times before over the years, even though the conda block has been in my bashrc all that time…strange that it is suddenly causing a problem…)
Regardless, thank you!
