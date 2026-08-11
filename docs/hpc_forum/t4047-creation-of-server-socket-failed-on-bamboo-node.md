# Creation of server socket failed on Bamboo node

- Source: https://hpc-community.unige.ch/t/4047

- Created: 2025-08-13T07:57:28.836Z

- Posts: 7

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Simon.Hug (2025-08-13T07:57:28.872Z)

Dear all
I encountered an annoying problem when running an array of jobs on Bamboo in which R carries out a parallel bootstrap. While I have done this multiple times without problems, yesterday (August 12th) one of these jobs halted with the following (R) error message:
Error in serverSocket(port = port) :
creation of server socket failed: port 11835 cannot be opened
Calls: makeCluster → makePSOCKcluster → serverSocket
As all other jobs in the array are working perfectly fine, I’m wondering whether there is/was a problem on one of the Bamboo nodes (I use the public-cpu partition).
Thanks for any insights, simon


## Post 2 by @Gael.Rossignol (2025-08-13T08:17:56.452Z)

Dear Simon,
Could you please give me the job number to check the problem?
Best regards,


## Post 3 by @Simon.Hug (2025-08-13T11:54:34.293Z)

The job number was  2680491_7
best, simon


## Post 4 by @Gael.Rossignol (2025-08-14T13:18:36.631Z)




## Post 5 by @Yann.Sagon (2025-08-15T07:29:06.218Z)

Without knowing your sbatch and more details, it is hard to debug. Maybe you started two jobs on the same compute node and both tried to open the same port number? By the way, why do you need to open a port?


## Post 6 by @Simon.Hug (2025-08-20T07:16:09.265Z)

The sbatch looks as follows:
```
#SBATCH --time=95:00:00
#SBATCH --cpus-per-task=12
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=3000 # in MB
#SBATCH -o myjob-%A\_%a.out

module load GCC/11.3.0  OpenMPI/4.1.4 R/4.2.1
module load rgdal/1.6-6

OUTFILE=srcvr_mc9s$SLURM_ARRAY_TASK_ID

Rscript --slave srcvr_mc6s_n.r > $OUTFILE
```
The part of the code in the R-file srcvr_mc6s_n.r that led to the error is the following:
```
(cores<-as.numeric(Sys.getenv(“SLURM_CPUS_PER_TASK”)))
cl ← makeCluster(cores)
setDefaultCluster(cl)
```
To reiterate: this code has worked perfectly fine multiple times and just once it produced the error mentioned before.
best wishes, simon


## Post 7 by @Yann.Sagon (2025-08-20T08:04:31.958Z)

Thanks for sharing your code. I did some investigation and testing.
If I fix the port number in the `makeCluster` function and launch multiple jobs on the same compute node, I have the same error as you do, this is normal and my first intuition that you had a conflict because you were running multiple R jobs on a compute node is probably correct.
```
Error in serverSocket(port = port) :
  creation of server socket failed: port 11000 cannot be opened
Calls: makeCluster -> makePSOCKcluster -> serverSocket
Execution halted
srun: error: cpu329: task 0: Exited with exit code 1
```
It appears that R is trying to open random port to avoid such conflict but maybe the mechanism isn’t very reliable or it was a race condition as you started multiple jobs simultaneously using job arrays.
The easiest workaround is to use `FORK` backend for cluster instead of the default “PSOCK”. In this case no port is opened. If you are interested to spread your work in more than one compute node, you should check the `MPI` backend.
```
cl <- makeCluster(cores, "FORK")
```
my 2 cents: It seems `--slave` is obsolete and can be replaced by `--vanilla` to run R in your sbatch script.


## Post 8 by @Simon.Hug (2025-08-20T08:09:22.186Z)

Thanks for this investigating and proposed solution: I will check this out. best wishes, simon
