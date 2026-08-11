# Engine crash Yggdrasil

- Source: https://hpc-community.unige.ch/t/3954

- Created: 2025-05-19T10:35:37.218Z

- Tags: yggdrasil

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Elisa.Denis (2025-05-19T10:35:37.266Z)

Username: $denisel0
Cluster: $Yggdrasil
I am trying to use ipyrad Welcome to ipyrad — ipyrad documentation[Welcome to ipyrad — ipyrad documentation](https://ipyrad.readthedocs.io/en/master/) but I get the following error. Prior to the Yggdrasil update, it worked well.
## ipyrad [v.0.9.105]
## Interactive assembly and analysis of RAD-seq data
Parallel connection |
Step 1: Demultiplexing fastq data to Samples
Encountered an Error.
Message: This operation requires engines. Try client.wait_for_engines(n) to wait for engines to register.
I also tried to debug as mentionned on the ipyrad FAQ with the following commands but the error remains the same.
$ echo “# Prevent ipyparallel engines from dying in a headless environment” >> ~/.bashrc
$ echo “export QT_QPA_PLATFORM=offscreen” >> ~/.bashrc
$ echo “export QT_QPA_PLATFORM=offscreen” >> ~/.profile
$ source ~/.bashrc
$ source ~/.profile
$ env | grep QT
These are the parameters I am using in my script:
#SBATCH --job-name ipyrad
#SBATCH --error %j.e.
#SBATCH --output %j.o.
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 12
#SBATCH --partition shared-bigmem
#SBATCH --time 10:00:00
Did anyone had similar issues with engines?
Many thanks
Elisa


## Post 2 by @Gael.Rossignol (2025-05-28T15:14:04.549Z)

Dear Elisa,
Could you please give me id of failled jobs to check what happens?
Thanks,
edit: and the path to your sbatch please.
