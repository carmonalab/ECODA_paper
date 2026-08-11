# Obscure MPI warning on a specific node

- Source: https://hpc-community.unige.ch/t/3377

- Created: 2024-03-16T17:58:05.972Z

- Posts: 3

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Sara.Ricossa (2024-03-16T17:58:06.036Z)

Hi everyone! I am using baobab for semidefinite optimization with a program called `sdpb`. The specifics don’t really matter to this question, other than the fact that it relies on `OpenMPI/4.1.4`. If however you need to see a script or any other information do let me know.
I just stumbled upon this UCX warning that crashes the software I’m using:
```
[1710610763.983411] [cpu319:2620412:0]           ib_md.c:1234 UCX  WARN  IB: ibv_fork_init() was disabled or failed, yet a fork() has been issued.
[1710610763.983419] [cpu319:2620412:0]           ib_md.c:1235 UCX  WARN  IB: data corruption might occur when using registered memory.
```
The warning is repeated many times until `sdpb` crashes with this error message:
```
Process 9 caught error message:
in allocate_blocks() at ../src/sdp_solve/Block_Info/allocate_blocks.cxx:72: 
  Assertion '!block_indices.empty()' failed:
    No SDP blocks were assigned to rank=9. node=9 node_rank=0
```
which I believe is saying that the node 319 (node number 9 out of the ten I was using) is not behaving well.
To verify this I tried running the same batch script with less nodes. When I am not assigned node 319, everything runs smoothly. As soon as node 319 is assigned I run into this error.
I have tried looking up what this warning means but it is still very unclear to me. I was hoping some of you might know what it means and how I could avoid running into it.
Secondary question is who should I contact in the future if I notice a certain node is not behaving well and might need rebooted?
Thank you very much for your time and have a great rest of your weekend. :smiling_face:
Edit: more nodes are now presenting the same issue…
Other edit: Added the version of OpenMPI I’m using. Unfortunately I cannot use an earlier version since it clashes with other dependencies. Also, the errors started appearing two days ago and got worse over the weekend.


## Post 2 by @Gael.Rossignol (2024-04-03T12:53:37.765Z)

Dear Sergio,
I see that cpu319 has a memory problem thanks for reporting it!
```
[Apr 1 15:47] core: [Hardware Error]: Machine check events logged
[  +0.006110] [Hardware Error]: Corrected error, no action required.
[  +0.006352] [Hardware Error]: CPU:3 (17:31:0) MC17_STATUS[Over|CE|MiscV|AddrV|-|-|SyndV|CECC|-|-|Scrub]: 0xdc2041000000011b
[  +0.011297] [Hardware Error]: Error Addr: 0x00000004010f4080
[  +0.005835] [Hardware Error]: IPID: 0x0000009600650f00, Syndrome: 0x112201010a800301
[  +0.007913] [Hardware Error]: Unified Memory Controller Ext. Error Code: 0, DRAM ECC error.
[  +0.010255] [Hardware Error]: cache level: L3/GEN, tx: GEN, mem-tx: RD
```
We will drain the node and replace the hardware.
Best regards,


## Post 3 by @Gael.Rossignol (2024-04-03T13:33:55.791Z)

Sara.Ricossa:
> `ibv_fork_init() was disabled or failed`
Regarding errors on other nodes, could you please share you sbatch to check environment used?
Best regards,
