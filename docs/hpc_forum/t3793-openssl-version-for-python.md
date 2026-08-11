# OpenSSL version for python

- Source: https://hpc-community.unige.ch/t/3793

- Created: 2025-01-22T13:46:05.675Z

- Tags: baobab

- Posts: 8

- Category: 10

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Maria.Reva (2025-01-22T13:46:05.728Z)

Hi,
I am trying to use the `bertopic` package in Python, but I’m facing issues because of a mismatch in OpenSSL versions. Here’s the situation:
- Python’s Built-In OpenSSL: Python 3.8.6 is built with OpenSSL 1.0.2o-fips (2018).
- Verified via: `python -c "import ssl; print(ssl.OPENSSL_VERSION)"`
- System OpenSSL: The system’s OpenSSL version is 1.1.1k.
- Verified via: `openssl version`
`bertopic` requires OpenSSL 1.1.1+, which is incompatible with Python’s current OpenSSL.
I would like to run Python with a version of OpenSSL compatible with the `bertopic` package (1.1.1+).
Thanks for your help,
Maria


## Post 2 by @Adrien.Albert (2025-01-22T16:39:15.756Z)

Hi Maria
For future request pleas post on HPC issues[HPC issues](https://hpc-community.unige.ch/c/hpc-support/hpc-issues/14)  catergorie and use the provided template.
Could you please give me the sbatch and the code to reproduce the error ?
How did you install `bertopic`


## Post 3 by @Maria.Reva (2025-01-23T14:41:29.323Z)

Hi,
Ive loaded the module and active my venv (the batch is in /home/users/r/reva/SFN_LNP/scripts/bert.sh) :
module load GCCcore/10.2.0 Python/3.8.6
. /home/users/r/reva//env_name/bin/activate
The error pops up when I simply do import, take aside running the code, just import fails : from bertopic import BERTopic
bertopic was installed directly in my venv via pip install.
Let me know if I can provide any other info.
Thanks,
Maria


## Post 4 by @Maria.Reva (2025-01-24T13:47:59.636Z)

Hi,
let me know if I can provide more information.
Thank you,
Maria


## Post 5 by @Adrien.Albert (2025-01-24T13:50:13.128Z)

I am sorry for the wait time, we are on urgent issue impacting serverals users :confused:


## Post 6 by @Maria.Reva (2025-01-28T07:31:43.888Z)

Hi Adrien,
No problems. Please let me know if you had time to look into it.
Thank you,
Maria


## Post 7 by @Yann.Sagon (2025-01-28T09:51:43.972Z)

@Maria.Reva[@Maria.Reva](https://hpc-community.unige.ch/u/maria.reva)
I’ve recompiled Python 3.8.6 with more recent OpenSSL version:
Before:
```
(baobab)-[sagon@login1 ~]$ python -c "import ssl; print(ssl.OPENSSL_VERSION)"
OpenSSL 1.0.2o-fips  27 Mar 2018
```
After:
```
(baobab)-[sagon@login1 ~]$ python -c "import ssl; print(ssl.OPENSSL_VERSION)"
OpenSSL 1.1.1k  FIPS 25 Mar 2021
```
By the way: is there a reason to use this old Python version instead of the latest we provide?
Best


## Post 8 by @Maria.Reva (2025-01-28T14:48:49.671Z)

Thanks a lot!
In my experience 3.8 worked well when prototyping on my personal machine, python version >3.9 resulted in the incompatibilities with others lib that we use along the way.
Thanks again,
Maria
