# [tutorial] run an old binary on the cluster

- Source: https://hpc-community.unige.ch/t/4278

- Created: 2026-04-16T12:55:29.182Z

- Posts: 1

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yann.Sagon (2026-04-16T12:55:29.315Z)

# Running Legacy 32-bit Binaries on Rocky Linux 9 with Apptainer
If you need to run an ancient 32-bit ELF binary (e.g., compiled for Linux Kernel 2.2.5) on a modern 64-bit system like Rocky Linux 9, the cleanest solution is to use Apptainer with a 32-bit CentOS 7 base.
## The Problem
Modern 64-bit distributions lack the 32-bit dynamic loader (/lib/ld-linux.so.2) and legacy versions of shared libraries (like libpng.so.2 or libXmu.so.6) required by software from the late 90s/early 2000s.
## Step 1: Create the Apptainer Definition File
Create a file named `xrec_env.def`. This configuration uses the i386 architecture and fixes the CentOS 7 repositories, which moved to the “Vault” archives after June 2024.
```
Bootstrap: docker
From: i386/centos:7

%post
    # 1. Point Yum to the CentOS Vault (CentOS 7 is EOL)
    cat <<EOF > /etc/yum.repos.d/CentOS-Base.repo
[base]
name=CentOS-7.9.2009 - Base - vault.centos.org
baseurl=https://vault.centos.org/altarch/7.9.2009/os/i386/
gpgcheck=1
gpgkey=file:///etc/pki/rpm-gpg/RPM-GPG-KEY-CentOS-7

[updates]
name=CentOS-7.9.2009 - Updates - vault.centos.org
baseurl=https://vault.centos.org/altarch/7.9.2009/updates/i386/
gpgcheck=1
gpgkey=file:///etc/pki/rpm-gpg/RPM-GPG-KEY-CentOS-7
EOF

    # 2. Install required 32-bit libraries
    yum clean all
   yum -y install glibc.i686 libgcc.i686 libstdc++.i686 libXmu.i686 libXext.i686 libXp.i686 libpng-1.5.13-8.el7.i686

    # 3. Fix legacy library versions (The Symlink Hack)
    # If the binary expects libpng.so.2, we link it to the available libpng15
    ln -s /usr/lib/libpng15.so.15 /usr/lib/libpng.so.2

%runscript
    exec ./xrec "$@"
```
# Step 2: Build the Image
Build the SIF (Apptainer Image Format) file.
```
apptainer build xrec_env.sif xrec_env.def
```
## Step 3: Execution and Verification
You can now run your legacy binary or inspect it within the container environment.
### Check dependencies:
```
apptainer exec xrec_env.sif ldd ./xrec
```
### Run the program:
```
apptainer run xrec_env.sif
```
## Summary
- Isolation: Avoids installing 32-bit compatibility layers directly on the Rocky 9 host.
- Architecture: Forcing i386/ ensures the entire container stack is 32-bit.
- Vault: Essential for CentOS 7 since its end-of-life in mid-2024.
