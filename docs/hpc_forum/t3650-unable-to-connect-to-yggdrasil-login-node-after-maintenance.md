# Unable to connect to Yggdrasil login node after maintenance

- Source: https://hpc-community.unige.ch/t/3650

- Created: 2024-07-17T14:07:25.286Z

- Tags: yggdrasil

- Posts: 7

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Elliot.Jiwani-Brown (2024-07-17T14:07:25.286Z)

Hi, I am unable to login to the yggdrasil server using my same X2Go session. Should I delete and make a new session or is there an update I need to make now that yggdrasil is updated?
Also, your documentation website that has all the information on how to use the HPC doesn’t appear to be working, do you know why this is?
Thanks
Elliot


## Post 2 by @Yann.Sagon (2024-07-17T14:22:21.018Z)

Dear  @Elliot.Jiwani-Brown[@Elliot.Jiwani-Brown](https://hpc-community.unige.ch/u/elliot.jiwani-brown) we completely reinstalled the login node, so yes for shure, your session was killed, try to start a new one.
We had issue with our doc, it was only reachable from unige network but this is solved, I tried right now. Do you still have the issue? If yes, can you share more details about that please?
Best


## Post 3 by @Elliot.Jiwani-Brown (2024-07-17T14:43:37.652Z)

Dear @Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon), it still appears not to be working. I get the error:
```
Warning: Protocol mismatch or no X authentification
Session: Terminating session at 'Wed Jul 17 16:40:32 2024'
Info: Your session was closed before reaching a useable state.
Info: This can be due to the local X server refusing access to the client.
Info: Please check authorization provided by the remote X application.
SessionL Session terminated at 'Wed Jul 17 16:40:32 2024'
```
I can log in with an ssh session but cannot open anything that may be graphic e.g. gedit (even using -X)
Thanks
Elliot


## Post 4 by @Elliot.Jiwani-Brown (2024-07-18T08:56:07.585Z)

@Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon) do you have any folllow ups on this? I am still experiencing the same issue.
Also when using ssh -X, I get the following error (which is also similar to the X2Go error):
`/usr/bin/xauth:  error in locking authority file /home/users/j/jiwanibr/.Xauthority`
For the documents link, I still could only access it using an unige VPN


## Post 5 by @Adrien.Albert (2024-07-18T15:18:08.321Z)

Hi @Elliot.Jiwani-Brown[@Elliot.Jiwani-Brown](https://hpc-community.unige.ch/u/elliot.jiwani-brown)
It seems that your quota on Yggdrasil has been reached, and is blocking the creation of files necessary for x2go to function.
I’ve increased your quota until July 27, 2024. Once you’ve logged in, clean out your personal directory to prevent the problem from recurring…
Best Regards


## Post 6 by @Elliot.Jiwani-Brown (2024-07-18T15:37:11.297Z)

@Adrien.Albert[@Adrien.Albert](https://hpc-community.unige.ch/u/Adrien.Albert) Thank you very much for allowing the quota to fix this. I have cleared my personal folder so hopefully this should not be a problem in future. X2Go is working well now.
I assume the quota is not extended to the scratch folder and here there is more memory allowed?


## Post 7 by @Adrien.Albert (2024-07-18T15:39:54.532Z)

Hi @Elliot.Jiwani-Brown[@Elliot.Jiwani-Brown](https://hpc-community.unige.ch/u/elliot.jiwani-brown)
Greate to know eveything is working
Elliot.Jiwani-Brown:
> I assume the quota is not extended to the scratch folder and here there is more memory allowed?
You can find more information about quota in our documentation:
https://doc.eresearch.unige.ch/hpc/storage_on_hpc?s[]=quota#cluster_storage[https://doc.eresearch.unige.ch/hpc/storage_on_hpc?s[]=quota#cluster_storage](https://doc.eresearch.unige.ch/hpc/storage_on_hpc?s%5B%5D=quota#cluster_storage)
For now, we only enable inode quota on `/srv/beegfs/scratch`.
