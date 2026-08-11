# Wolframscript linking issues on Baobab

- Source: https://hpc-community.unige.ch/t/3362

- Created: 2024-03-11T19:02:08.175Z

- Posts: 2

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Bharathkumar.Radhakrishnan (2024-03-11T19:02:08.239Z)

There is a slight weirdness about the installation of Mathematica 12.3.0 on the baobab cluster that I would like to point out.
If I am not mistaken, every distribution of Mathematica comes with the binary called wolframscript. Indeed, if you check `/unige/mathematica_11.1/bin` you have the symlink `wolframscript -> /unige/mathematica_11.1/SystemFiles/Kernel/Binaries/Linux-x86-64/wolframscript`. It exists in v12.3.0 as well. However you do not have a symlink to wolframscript in `/unige/mathematica_12.3/bin`, even though `/unige/mathematica_12.3/SystemFiles/Kernel/Binaries/Linux-x86-64/wolframscript` exists as a fully functional binary.
Now when using `ml load`, only `/unige/mathematica_12.3/bin` gets added to `PATH`, and thus `which wolframscript` would give the error that it could not find it, and one needs to manually append `PATH` with the exact binary to be able to use wolframscript from v12.3.0.
The reason why I believe it is important to point this out is that I have noticed a few undocumented instabilities with wolframscript in v11.1.1, which seems to have been solved in v12.3.0, and thus if I specifically want to use v12.3.0, the only way is to first load 12.3.0 and then manually add the path to the binary in `PATH`, which is a perfectly functional but non-elegant solution.


## Post 2 by @Gael.Rossignol (2024-04-03T13:38:06.237Z)

Hello,
It seems the issue has been resolved :white_check_mark:
> (baobab)-[rossigng@login2 ~]$ ml Mathematica/12.3.0
> (baobab)-[rossigng@login2 ~]$ which Mathematica
> /unige/mathematica_12.3/bin/Mathematica
> (baobab)-[rossigng@login2 ~]$ ll /unige/mathematica_12.3/bin/Mathematica
> lrwxrwxrwx 1 nobody nobody 47 Jun  4  2021 /unige/mathematica_12.3/bin/Mathematica → /unige/mathematica_12.3/Executables/Mathematica
> (baobab)-[rossigng@login2 ~]$ ll /unige/mathematica_12.3/bin/
> total 256
> lrwxrwxrwx 1 nobody nobody 40 Jun  4  2021 math → /unige/mathematica_12.3/Executables/math
> lrwxrwxrwx 1 nobody nobody 47 Jun  4  2021 mathematica → /unige/mathematica_12.3/Executables/mathematica
> lrwxrwxrwx 1 nobody nobody 47 Jun  4  2021 Mathematica → /unige/mathematica_12.3/Executables/Mathematica
> lrwxrwxrwx 1 nobody nobody 46 Jun  4  2021 MathKernel → /unige/mathematica_12.3/Executables/MathKernel
> lrwxrwxrwx 1 nobody nobody 39 Jun  4  2021 mcc → /unige/mathematica_12.3/Executables/mcc
> lrwxrwxrwx 1 nobody nobody 43 Jun  4  2021 wolfram → /unige/mathematica_12.3/Executables/wolfram
> lrwxrwxrwx 1 nobody nobody 49 Jun  4  2021 WolframKernel → /unige/mathematica_12.3/Executables/WolframKernel
> lrwxrwxrwx 1 nobody nobody 49 Mar 14 15:02 wolframscript → /unige/mathematica_12.3/Executables/wolframscript
Could you please confirm all is now ok?
Best regards,
