# Request for ASHS Software Installation

- Source: https://hpc-community.unige.ch/t/3972

- Created: 2025-06-04T14:17:47.604Z

- Tags: bamboo

- Posts: 3

- Category: 10

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @aybuke.calik (2025-06-04T14:17:47.652Z)

Hello,
I would like to perform brain segmentation using structural MRI data, which requires installing the ASHS software (as outlined here: ASHS Documentation - Installation[ASHS Documentation - Installation](https://sites.google.com/view/ashs-dox/local-ashs/installation)). I’m aware that installing software directly in home directories is not permitted.
Would it be possible for the HPC team to install ASHS on Bamboo? Alternatively, is there another recommended way to access or use this software on the HPC system?
Thank you for your help in advance.
Best,
Aybüke Çalık


## Post 2 by @aybuke.calik (2025-06-05T08:12:03.699Z)

Moreover, I need installation of MATLAB 2019b runtime to be able to use FreeSurfer.
ERROR: cannot find Matlab 2019b runtime in location:
/opt/ebsofts/FreeSurfer/7.3.2-centos7_x86_64/MCRv97
It is looking for either:
- bin/glnxa64/libmwlaunchermain.so    (Linux 64b) or*
- bin/maci64/libmwlaunchermain.dylib (Mac 64b)*
The hippocampal/amygdala and brainstem modules require the (free) Matlab runtime.
You will need to download the Matlab Compiler Runtime (MCR) for Matlab 2019b.
To do so, please run the following command (you might need root permissions):
fs_install_mcr R2019b
For further details, please visit MatlabRuntime - Free Surfer Wiki[MatlabRuntime - Free Surfer Wiki](https://surfer.nmr.mgh.harvard.edu/fswiki/MatlabRuntime)


## Post 3 by @Gael.Rossignol (2025-06-11T14:56:29.672Z)

Hello,
I compile a new version of FreeSurfer :
FreeSurfer/7.4.1-centos8_x86_64
Could you test if it’s working with this one?
Best regards,
