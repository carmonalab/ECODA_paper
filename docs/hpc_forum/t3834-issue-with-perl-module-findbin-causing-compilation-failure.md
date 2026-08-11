# Issue with Perl Module "FindBin" Causing Compilation Failure

- Source: https://hpc-community.unige.ch/t/3834

- Created: 2025-02-20T18:02:16.533Z

- Tags: yggdrasil

- Posts: 4

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Siddharth.Bhatnagar (2025-02-20T18:02:16.587Z)

Hello!
Username: `bhatnags`
Cluster: `yggdrasil`
I am writing to report an issue I encountered while trying to compile my model (a Global Climate Model, GCM):
Here are some relevant details:
- Experience: I have been using the model since 2022 on the cluster without any problem. The GCM uses modules that I’ve specified in my bash profile (GCC, NetCDF, etc.), which remain unchanged. The last time I compiled the model was in October 2024 → this was successful. But then, I tried doing it again a couple days ago and had this issue.
- Command Executed: I executed the following command from `/home/users/b/bhatnags/trunk_r3397/LMDZ.COMMON`:
- `./makelmdz_fcm -arch BAO_YGG -p std -d 64x48x30 newstart`
- Error Message: The compilation failed in 2-3 seconds with the following error message:
```
Can't locate FindBin.pm in @INC (you may need to install the FindBin module) (@INC contains: /home/users/b/bhatnags/FCM_V1.2/lib /usr/local/lib64/perl5/5.32 /usr/local/share/perl5/5.32 /usr/lib64/perl5/vendor_perl /usr/share/perl5/vendor_perl /usr/lib64/perl5 /usr/share/perl5) at /home/users/b/bhatnags/FCM_V1.2/lib/Fcm/Config.pm line 25.
BEGIN failed--compilation aborted at /home/users/b/bhatnags/FCM_V1.2/lib/Fcm/Config.pm line 25.
Compilation failed in require at /home/users/b/bhatnags/FCM_V1.2/bin/fcm line 32.
BEGIN failed--compilation aborted at /home/users/b/bhatnags/FCM_V1.2/bin/fcm line 32.
```
- Relevant Files and Paths: FCM executable: `/home/users/b/bhatnags/FCM_V1.2/bin/fcm` | Configuration file: `/home/users/b/bhatnags/FCM_V1.2/lib/Fcm/Config.pm`
- Additional Information: I suspect that a recent update or change in the system (possibly affecting Perl or environment variables such as `PERL5LIB`) might have altered the include paths, causing Perl to be unable to locate its core module FindBin.
I would appreciate any help in dealing with this. Please let me know if I can provide more details.
Thank you very much for the help in advance!
Cheers,
Siddharth


## Post 2 by @Siddharth.Bhatnagar (2025-02-24T17:26:33.328Z)

It turned out that indeed the `PERL5LIB` environment variable was unset (but I had never set it previously) and Perl was strangely not available in the `@INC` path anymore.
The fix that works for me is to manually set `PERL5LIB`:
`export PERL5LIB=/opt/ebsofts/Perl/5.32.1-GCCcore-10.3.0-minimal/lib64/perl5:$PERL5LIB`
You can add this line to your `.bash_profile` to make the change persistent.


## Post 3 by @Adrien.Albert (2025-02-25T15:19:16.222Z)

Hi @Siddharth.Bhatnagar[@Siddharth.Bhatnagar](https://hpc-community.unige.ch/u/siddharth.bhatnagar)
Did you load the module Perl/5.32.1-minimal before your compilation ?


## Post 4 by @Siddharth.Bhatnagar (2025-04-16T09:36:47.124Z)

Hi @Adrien.Albert[@Adrien.Albert](https://hpc-community.unige.ch/u/Adrien.Albert)
Sorry, I just saw this message.
But yes, manually loading Perl did the trick – only then could I run my compilation command successfully.
