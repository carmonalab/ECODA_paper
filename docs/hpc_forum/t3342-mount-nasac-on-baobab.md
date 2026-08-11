# Mount nasac on baobab

- Source: https://hpc-community.unige.ch/t/3342

- Created: 2024-02-28T13:39:53.821Z

- Posts: 5

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Nada.Kojovic (2024-02-28T13:39:53.894Z)

Hi all
I am trying to mount nasac-m2 to baobab and I would need your help before making any mistake.
On Nasac we have a lot of heavy files that I would like to avoid loading to scratch and if possible process directly there. I see there is a mount for nasac on srv/nasac.unige.ch but how can I have my partition of interest there (m-AutismLab)? Thank you!
Nada


## Post 2 by @Yann.Sagon (2024-02-29T13:51:14.666Z)

Hi @Nada.Kojovic[@Nada.Kojovic](https://hpc-community.unige.ch/u/nada.kojovic), did you check here[here](https://doc.eresearch.unige.ch/hpc/storage_on_hpc#access_external_storage)?


## Post 3 by @Nada.Kojovic (2024-03-01T14:11:40.735Z)

Dear Yann
Thank you. I did actually but I have nothing in /run/user/your_uid/gvfs/my uid is 296556 if not wrong.
I do have something here it can navigate nasac directories but does not list files in them: /home/kojovic/.gvfs/smb-share:server=nasac-m2.unige.ch,share=m-autismlab
Many thanks for your help
Nada


## Post 4 by @Yann.Sagon (2024-03-04T09:48:42.487Z)

Are you sure the rights are set correctly in the NASAC? Are you able to access those files from windows with the same username? Maybe the user rights are set per group under windows and the group isn’t used in Baobab. Can you show a screenshot of advanced properties under windows of one of the file?
If unsure, please contact the storage team storage@unige.ch.


## Post 5 by @Nada.Kojovic (2024-03-19T13:35:04.624Z)

Hi Yann
Here is it:
image
image1156×224 83.5 KB
[image1156×224 83.5 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/1800527d242f082f6c2a3c7b3770ace099805010.png)
