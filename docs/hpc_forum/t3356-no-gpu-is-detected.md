# No GPU is detected

- Source: https://hpc-community.unige.ch/t/3356

- Created: 2024-03-07T08:29:26.859Z

- Posts: 2

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Jingze.Duan (2024-03-07T08:29:26.916Z)

If you are asking for help, try to provide information that can help us solve your issue, such as :
My jobs running on gpu008 always failed after a few second with the same error (on yggdrasil).
```
Executable:   /opt/ebsofts/GROMACS/2023.1-foss-2022a-CUDA-11.7.0/bin/gmx_mpi
Data prefix:  /opt/ebsofts/GROMACS/2023.1-foss-2022a-CUDA-11.7.0
Working dir:  /srv/beegfs/scratch/users/d/duanjin3/GdmCl_0mV
Command line:
  gmx_mpi mdrun -deffnm Ez_pos -s Ez_pos.tpr -v -nsteps -1 -pin on -nb gpu -cpi Ez_pos.cpt -noappend

Fatal error:
Cannot run short-ranged nonbonded interactions on a GPU because no GPU is
detected.
```
Screenshot from 2024-03-07 09-23-49
Screenshot from 2024-03-07 09-23-49915×55 11.1 KB
[Screenshot from 2024-03-07 09-23-49915×55 11.1 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/dbc926738f63665f2935780b12c55d3b96b865d0.png)
What should I do or how can I escape gpu008?


## Post 2 by @Gael.Rossignol (2024-04-04T09:25:35.259Z)

Hello,
I think this issue is solved, please read this post :white_check_mark:
No GPU is detected on yggdrasil[No GPU is detected on yggdrasil](https://hpc-community.unige.ch/t/no-gpu-is-detected-on-yggdrasil/3291/4) HPC Technical[HPC Technical](https://hpc-community.unige.ch/c/hpc-technical/5)
> Dear @Jingze.Duan[@Jingze.Duan](https://hpc-community.unige.ch/u/jingze.duan) 
Gromac has been recompiled.
Best regards,
