# Perl error: missing xlocale.h

- Source: https://hpc-community.unige.ch/t/3805

- Created: 2025-02-01T10:09:57.822Z

- Tags: yggdrasil

- Posts: 9

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Dominique.Eckert (2025-02-01T10:09:57.869Z)

Dear colleagues,
I am attempting to compile a code that includes some Perl scripts. Config was successful. However, when building the software, I encountered the following error:
```
/opt/ebsofts/Perl/5.34.0-GCCcore-11.3.0/lib/perl5/5.34.0/x86_64-linux-thread-multi/CORE/perl.h:922:13: fatal error: xlocale.h: No such file or directory
  922 | #   include <xlocale.h>
      |             ^~~~~~~~~~~
compilation terminated.
```
I looked this up and it appears that the perl.h library linked above contains an include statement for xlocale.h, which doesn’t exist. Can you please tell me how to deal with this?
Note: other Perl versions installed on Yggdrasil also appear to contain this include statement to the non-existent xlocale.h.
Cheers,
Dominique


## Post 2 by @Yann.Sagon (2025-02-03T08:20:38.038Z)

Dear  @Dominique.Eckert[@Dominique.Eckert](https://hpc-community.unige.ch/u/dominique.eckert)
This is the block in `perl.h`
```
#ifdef I_XLOCALE
#   include <xlocale.h>
#endif
```
`xlocale.h` is optional. Maybe you have to disable something in your configure script?
Did you do the configure step from the beginning on the cluster? What are you trying to build?
Best
Yann


## Post 3 by @Dominique.Eckert (2025-02-04T08:35:00.538Z)

Hi @Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon) , thanks for the response. I am trying to build HEASOFT v6.34
https://heasarc.gsfc.nasa.gov/lheasoft/[https://heasarc.gsfc.nasa.gov/lheasoft/](https://heasarc.gsfc.nasa.gov/lheasoft/)
I am on the main node of the cluster, and ran the configure script directly from there, linking PERL to the Perl5.34.0 version


## Post 4 by @Dominique.Eckert (2025-02-06T19:36:18.662Z)

Hi @Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon) . I sent a request to the HEASOFT helpdesk and got the following reply:
> Hi Dominique, The CFITSIO Perl interface code includes “perl.h” which was found on your system in /opt/ebsofts/Perl/5.34.0-GCCcore-11.3.0/lib/perl5/5.34.0/x86_64-linux-thread-multi/CORE/, and it’s that system perl.h file that is trying to include xlocale.h.
> So it sounds like there’s a problem in the Perl installation, i.e. that it includes a perl.h which tries to include a file that isn’t there.  Does Rocky Linux use the Red Hat package manager ‘dnf’?  I’ve tried to discern which dnf package might include xlocale.h but without success so far.  My best guess is that it’s in the “perl-devel” package.  If that doesn’t work, my best suggestion is to try to install and use a different version of Perl.
Can you please have a look?


## Post 5 by @Gael.Rossignol (2025-02-07T08:13:56.849Z)

Dominique.Eckert:
> > So it sounds like there’s a problem in the Perl installation, i.e. that it includes a perl.h which tries to include a file that isn’t there. Does Rocky Linux use the Red Hat package manager ‘dnf’? I’ve tried to discern which dnf package might include xlocale.h but without success so far. My best guess is that it’s in the “perl-devel” package. If that doesn’t work, my best suggestion is to try to install and use a different version of Perl.
Dear Dominique,
I’m actually recompiling HEASOFT, I will keep you informed as soon as I finish. Which perl are you using during your compilation?
Best regards,


## Post 6 by @Dominique.Eckert (2025-02-07T08:49:13.511Z)

Hi @Gael.Rossignol[@Gael.Rossignol](https://hpc-community.unige.ch/u/Gael.Rossignol) . I was using Perl 5.34 but as far as I can see, any relatively recent version of Perl would do.
Thanks for looking into this


## Post 7 by @Gael.Rossignol (2025-02-10T13:10:06.753Z)

Dear Dominique,
I compile HEASoft with EasyBuild :
```
(yggdrasil)-[rossigng@login1 ~]$ module spider HEASoft

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  HEASoft: HEASoft/6.34
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Description:
      HEASOFT is the HEASARC's software suite that currently encompasses both new and legacy mission-independent or multi-mission FTOOLS, the XANADU suite (XSPEC, XRONOS, and XIMAGE),
      XSTAR, and numerous mission-specific packages.

    You will need to load all module(s) on any one of the lines below before the "HEASoft/6.34" module is available to load.

      GCC/13.3.0

    Help:
      Description
      ===========
      HEASOFT is the HEASARC's software suite that currently encompasses both new and legacy mission-independent or multi-mission FTOOLS, the XANADU suite (XSPEC, XRONOS, and XIMAGE), XSTAR,
 and numerous mission-specific packages.

      More information
      ================
       - Homepage: https://heasarc.gsfc.nasa.gov/lheasoft/
```
You can load with command :
```
(yggdrasil)-[rossigng@login1 ~]$ ml GCC/13.3.0 GCC/13.3.0
```
Could you please test if this version is working fine?
Have a nice day,


## Post 8 by @Dominique.Eckert (2025-02-10T14:23:03.825Z)

That appears to work, thanks very much @Gael.Rossignol[@Gael.Rossignol](https://hpc-community.unige.ch/u/Gael.Rossignol) . Clearly this is a better solution than me compiling it on my own.
I’ll keep testing it and will let you know if everything works.


## Post 9 by @Dominique.Eckert (2025-02-13T08:43:55.392Z)

Hi @Gael.Rossignol[@Gael.Rossignol](https://hpc-community.unige.ch/u/Gael.Rossignol) . There are files missing in your HEASoft installation:
```
Failed to read continuum data from 3.0.9_coco.fits
Failed to read 3.0.9_coco.fits
```
I searched for this file (which is called apec_v3.0.9_coco.fits) and other similar files in the directory where HEASoft was installed and could not find them. These files are absolutely needed for a functional installation of HEASoft.
