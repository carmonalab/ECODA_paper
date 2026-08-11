# Unable to mount NASAC share on

- Source: https://hpc-community.unige.ch/t/3435

- Created: 2024-04-30T12:49:56.634Z

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Jade.Awada (2024-04-30T12:49:56.634Z)

Hey all!
I’m having the same problem.
I followed the steps on your documentation + the procedure you described here but it didn’t work.
When I mount a NasAc and try to access to .gvfs I obtain the following message: Transport endpoint is not connected.
Then I tried to unmount it with the following commands:
`gio mount -u smb://nasac-m2.isis.unige.ch/m-kaiser`, then
`pkill -u awada dbus, and rm -rf .gvfs`. I got this message: `rm: cannot remove '.gvfs': is a directory`.
Does anyone know what could be the problem?
Many thanks!!!
Best regards,
Jade


## Post 2 by @Yann.Sagon (2024-05-02T14:20:36.293Z)

Dear @Jade.Awada[@Jade.Awada](https://hpc-community.unige.ch/u/jade.awada)
on which cluster do you have the issue?
do you have the issue on the login node or on the compute node?
> I got this message: `rm: cannot remove '.gvfs': is a directory`.
This is weird. The command `rm -rf` is able to delete a directory. Please copy paste the exact command you type and the output.
