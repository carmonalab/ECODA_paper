# thread-MPI ranks (-ntmpi) is required to resolve this issue

- Source: https://hpc-community.unige.ch/t/3790

- Created: 2025-01-21T08:47:45.082Z

- Posts: 2

- Category: 1

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Chakradhara.Uppugunduri (2025-01-21T08:47:45.128Z)

Dear HPC team,
I hope this email finds you well. I am writing to inform you about an issue I encountered while running a GROMACS simulation on the updated environment with CUDA 11.7.
After the recent change from CUDA/12.3.0 to CUDA/11.7.0, I attempted to execute the following command:
gmx mdrun -v -deffnm md -cpi md.cpt
Unfortunately, the simulation terminated with a fatal error message:
```
When using GPUs, setting the number of OpenMP threads without specifying the number of ranks can lead to conflicting demands. Please specify the number of thread-MPI ranks as well (option -ntmpi).
```
The system suggests that specifying the number of thread-MPI ranks (-ntmpi) is required to resolve this issue. Additionally, the log indicates that the compiled SIMD is AVX2_256, whereas AVX_512 might be more optimal for this host.
Could you please assist in addressing this issue? If adjustments to the GROMACS configuration or reinstallation with different options are required, I kindly request your guidance or intervention to ensure compatibility with the current hardware and CUDA setup.
Thank you for your attention and support. Please let me know if you require any additional details or logs from my end.
Looking forward to your assistance.


## Post 2 by @Adrien.Albert (2025-01-22T09:25:06.510Z)

Hi @Chakradhara.Uppugunduri[@Chakradhara.Uppugunduri](https://hpc-community.unige.ch/u/chakradhara.uppugunduri)
Could you please provides your sbatch and the block code to reproduce your issue ?
