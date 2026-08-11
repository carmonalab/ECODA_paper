# OpenOnDemand SandBox App Permission denied

- Source: https://hpc-community.unige.ch/t/3307

- Created: 2024-02-13T15:53:32.125Z

- Posts: 3

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Victor.Ferat (2024-02-13T15:53:32.185Z)

Hey,
I’m trying to create a SandBox App on the OpenOnDemand platform. However, when I start the App, the `output.log` show an `Permission denied` error:
```
Script starting...
Generating connection YAML file...
/var/spool/slurmd/job7404075/slurm_script: line 128: /home/users/f/ferat/ondemand/data/sys/dashboard/batch_connect/dev/fMRIprepOnDemand/output/8033afc5-4c61-41fb-91f4-0fb446d4d634/script.sh: Permission denied
Cleaning up...
```
You can find the full App definition at this address:
GitHub[GitHub](https://github.com/vferat/fMRIprepOnDemand)
### GitHub - vferat/fMRIprepOnDemand[GitHub - vferat/fMRIprepOnDemand](https://github.com/vferat/fMRIprepOnDemand)
Contribute to vferat/fMRIprepOnDemand development by creating an account on GitHub.
And an example of a failed run HERE[HERE](https://ondemand.baobab.hpc.unige.ch/pun/sys/dashboard/files/fs/home/users/f/ferat/ondemand/data/sys/dashboard/batch_connect/dev/fMRIprepOnDemand/output/8033afc5-4c61-41fb-91f4-0fb446d4d634)
Might it be due to the linux authorizations `644` on the file ? Do you know if it is something I need to configure from the App side or if it’s a server side configuration ?
Thanks for your help,


## Post 2 by @Victor.Ferat (2024-02-15T12:16:29.058Z)

I solve the issue by adding read,write and execution permissions to the `ondemand` folder.
```
chmod -R u+rwx ./ondemand
```
Maybe it is somehing that should be set by default to avoid users having issue starting SandBox Apps


## Post 3 by @Adrien.Albert (2024-02-21T08:41:31.030Z)

Dear @Victor.Ferat[@Victor.Ferat](https://hpc-community.unige.ch/u/victor.ferat)
I am glade to see your contribution on OpenOnDemand. Let us know if your app need to be shared with the community :slight_smile:
