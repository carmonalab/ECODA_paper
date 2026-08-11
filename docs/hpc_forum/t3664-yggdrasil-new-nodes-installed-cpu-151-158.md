# [Yggdrasil] New nodes installed cpu[151-158]

- Source: https://hpc-community.unige.ch/t/3664

- Created: 2024-10-01T13:24:52.410Z

- Posts: 1

- Category: 6

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @remy.ressegaire (2024-10-01T13:24:52.441Z)

Dear users,
We have installed 8 new nodes on Yggdrasil :
```
(yggdrasil)-[root@login1 ~]$ sinfo -N --format="%10N %.6D %20P %.4c %.8z %.6m %.8d %G %.57f" -n cpu[151-158]
NODELIST    NODES PARTITION            CPUS    S:C:T MEMORY TMP_DISK GRES                                            AVAIL_FEATURES
cpu151          1 shared-cpu            128   2:64:1 515500   150000 (null)                                EPYC-7763,AMD_EPYC_7763,V8
cpu151          1 private-astro-pairs   128   2:64:1 515500   150000 (null)                                EPYC-7763,AMD_EPYC_7763,V8
cpu152          1 shared-cpu            128   2:64:1 515500   150000 (null)                                EPYC-7763,AMD_EPYC_7763,V8
cpu152          1 private-astro-pairs   128   2:64:1 515500   150000 (null)                                EPYC-7763,AMD_EPYC_7763,V8
cpu153          1 shared-cpu            128   2:64:1 515500   150000 (null)                                EPYC-7763,AMD_EPYC_7763,V8
cpu153          1 private-astro-pairs   128   2:64:1 515500   150000 (null)                                EPYC-7763,AMD_EPYC_7763,V8
cpu154          1 shared-cpu            128   2:64:1 515500   150000 (null)                                EPYC-7763,AMD_EPYC_7763,V8
cpu154          1 private-astro-pairs   128   2:64:1 515500   150000 (null)                                EPYC-7763,AMD_EPYC_7763,V8
cpu155          1 shared-cpu            128   2:64:1 515500   150000 (null)                                EPYC-7763,AMD_EPYC_7763,V8
cpu155          1 private-astro-pairs   128   2:64:1 515500   150000 (null)                                EPYC-7763,AMD_EPYC_7763,V8
cpu156          1 shared-cpu            128   2:64:1 515500   150000 (null)                                EPYC-7763,AMD_EPYC_7763,V8
cpu156          1 private-astro-pairs   128   2:64:1 515500   150000 (null)                                EPYC-7763,AMD_EPYC_7763,V8
cpu157          1 shared-cpu            128   2:64:1 515500   150000 (null)                                EPYC-7763,AMD_EPYC_7763,V8
cpu157          1 private-astro-pairs   128   2:64:1 515500   150000 (null)                                EPYC-7763,AMD_EPYC_7763,V8
cpu158          1 shared-cpu            128   2:64:1 515500   150000 (null)                                EPYC-7763,AMD_EPYC_7763,V8
cpu158          1 private-astro-pairs   128   2:64:1 515500   150000 (null)                                EPYC-7763,AMD_EPYC_7763,V8
```
