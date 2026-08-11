# Email notifications are not sent anymore

- Source: https://hpc-community.unige.ch/t/3500

- Created: 2024-06-19T09:27:53.636Z

- Posts: 3

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Leonhard.Unterlerchner (2024-06-19T09:27:53.676Z)

Dear HPC team and community,
I got an issue with email notifications.
#Username: unterlel
#Cluster: Baobab
Email notifications are not sent anymore. It used to work perfectly until the 6 of May. To get notifiactions on jobs I added the following lines in the Bash shell:
```
#SBATCH --mail-user=baobabunter@gmail.com
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
```
I have checked the spam folder of my email client, nothing in there. Anyone having the same issue? How can I resolve it?
Best,
Leonhard.


## Post 2 by @Adrien.Albert (2024-06-19T12:32:52.077Z)

Hi @Leonhard.Unterlerchner[@Leonhard.Unterlerchner](https://hpc-community.unige.ch/u/leonhard.unterlerchner)
it should work again. However, I’m not sure your current sbatch will work as expected.
From slurm sbatch documentation:
https://slurm.schedmd.com/sbatch.html[https://slurm.schedmd.com/sbatch.html](https://slurm.schedmd.com/sbatch.html)
> Multiple type values may be specified in a comma separated list.


## Post 3 by @Leonhard.Unterlerchner (2024-06-19T12:44:22.840Z)

Thank you for your work!
It works again.Thank you for your remark on my script. My Sbatch already send the right notifications but it will be more concise with your suggestion.
Best,
Leonhard.
