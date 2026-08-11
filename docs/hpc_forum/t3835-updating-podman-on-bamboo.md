# Updating Podman on Bamboo

- Source: https://hpc-community.unige.ch/t/3835

- Created: 2025-02-21T11:30:34.495Z

- Posts: 4

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Berk.Gercek (2025-02-21T11:30:34.539Z)

Hello all,
I’m happy to see Podman is available on the clusters, as I’d like to start building images on the cluster itself instead of on a separate machine (how it is necessary to do with Apptainer).
However on Bamboo it seems like Podman isn’t up to date and doesn’t work correctly, while on Yggdrasil a “hello world” container runs successfully.
Attempting to run a “hello world” container on bamboo yields the following error:
```
gercek@login1:~$ podman run quay.io/podman/hello
WARN[0000] The cgroupv2 manager is set to systemd but there is no systemd user session available
WARN[0000] For using systemd, you may need to log in using a user session
WARN[0000] Alternatively, you can enable lingering with: `loginctl enable-linger 334847` (possibly as root)
WARN[0000] Falling back to --cgroup-manager=cgroupfs
WARN[0000] Network file system detected as backing store.  Enforcing overlay option `force_mask="700"`.  Add it to storage.conf to silence this warning
WARN[0000] The cgroupv2 manager is set to systemd but there is no systemd user session available
WARN[0000] For using systemd, you may need to log in using a user session
WARN[0000] Alternatively, you can enable lingering with: `loginctl enable-linger 334847` (possibly as root)
WARN[0000] Falling back to --cgroup-manager=cgroupfs
Trying to pull quay.io/podman/hello:latest...
Getting image source signatures
Copying blob 81df7ff16254 done   |
Copying config 5dd467fce5 done   |
Writing manifest to image destination
ERRO[0002] Cleaning up container 9abb530e87ad03732041067c079bd53f71a0b283043c19159abe6d8efc96cd16: unmounting container 9abb530e87ad03732041067c079bd53f71a0b283043c19159abe6d8efc96cd16 storage: cleaning up container 9abb530e87ad03732041067c079bd53f71a0b283043c19159abe6d8efc96cd16 storage: unmounting container 9abb530e87ad03732041067c079bd53f71a0b283043c19159abe6d8efc96cd16root filesystem: replacing mount point "/home/gercek/.local/share/containers/storage/overlay/585848117a606fc5dbf711ee14e16e2c9f000750c97716e2487c2fbe4c916904/merged": file exists
Error: crun: open `/home/gercek/.local/share/containers/storage/overlay/585848117a606fc5dbf711ee14e16e2c9f000750c97716e2487c2fbe4c916904/merged/etc/resolv.conf`: No such file or directory: OCI runtime attempted to invoke a command that was not found
```
The podman version on Bamboo is
```
gercek@login1:~$ podman --version
podman version 4.9.4-rhel
```
while on Yggdrasil, where the container is able to run the expected output, it yields:
```
gercek@login1:~$ podman --version
podman version 5.2.2
```
Running the same container produces the expected output, but with the same set of warnings that occur on Bamboo:
```
gercek@login1:~$ podman run quay.io/podman/hello
WARN[0000] The cgroupv2 manager is set to systemd but there is no systemd user session available
WARN[0000] For using systemd, you may need to log in using a user session
WARN[0000] Alternatively, you can enable lingering with: `loginctl enable-linger 334847` (possibly as root)
WARN[0000] Falling back to --cgroup-manager=cgroupfs
WARN[0000] The cgroupv2 manager is set to systemd but there is no systemd user session available
WARN[0000] For using systemd, you may need to log in using a user session
WARN[0000] Alternatively, you can enable lingering with: `loginctl enable-linger 334847` (possibly as root)
WARN[0000] Falling back to --cgroup-manager=cgroupfs
!... Hello Podman World ...!

         .--"--.
       / -     - \
      / (O)   (O) \
   ~~~| -=(,Y,)=- |
    .---. /`  \   |~~
 ~/  o  o \~~~~.----. ~~
  | =(X)= |~  / (O (O) \
   ~~~~~~~  ~| =(Y_)=-  |
  ~~~~    ~~~|   U      |~~

Project:   https://github.com/containers/podman
Website:   https://podman.io
Desktop:   https://podman-desktop.io
Documents: https://docs.podman.io
YouTube:   https://youtube.com/@Podman
X/Twitter: @Podman_io
Mastodon:  @Podman_io@fosstodon.org
```
Would it be possible to update the runtime on Bamboo (and Baobab if it’s in the same version)? Also, it would be good to change the config for podman to something other than `cgroupv2=systemd` to suppress the warning that occurs on both platforms.
Thanks for the help!
Berk


## Post 2 by @Berk.Gercek (2025-03-05T15:03:26.332Z)

Just wanted to bump this thread to know if podman will be updated or whether I should stick to singularity on bamboo?


## Post 3 by @Gael.Rossignol (2025-03-06T08:47:38.942Z)

Berk.Gercek:
> cgroupv2=systemd
Dear Berk,
Sorry for delay on the answer, we actually update each clusters avery 3 month. Next maintenance of Bamboo will occurs :
Bamboo scheduled maintenance: 25-27 March 2025 - HPC Announce - HPC Community[Bamboo scheduled maintenance: 25-27 March 2025 - HPC Announce - HPC Community](https://hpc-community.unige.ch/t/bamboo-scheduled-maintenance-25-27-march-2025/3846)
So this will correct issue with podman and we will be up to date.
I will check if I can change podman config.
Best regards,


## Post 4 by @Berk.Gercek (2025-03-10T13:08:46.426Z)

Thanks Gael, no worries, I know I submitted this during a maintenance so it likely fell under the radar.
I’ll keep an eye out then!
