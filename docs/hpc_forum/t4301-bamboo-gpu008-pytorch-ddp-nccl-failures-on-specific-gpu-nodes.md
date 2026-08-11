# [bamboo][gpu008] PyTorch DDP/NCCL failures on specific GPU nodes

- Source: https://hpc-community.unige.ch/t/4301

- Created: 2026-06-06T21:20:21.292Z

- Tags: bamboo

- Posts: 1

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @maciej.falkiewicz (2026-06-06T21:20:21.368Z)

Hi HPC Team,
I would like to report a recurring issue with PyTorch distributed jobs on gpu008 on Bamboo. The same training code and environment work correctly on other nodes, but on this problematic node, the job hangs/fails during DDP initialization, before actual training starts.
Observed behavior:
- The job is launched with `torchrun` using 4 local ranks.
- All ranks appear to start correctly and construct the model.
- Each rank prints the same model size: `number of parameters: 26.75M`.
- The failure happens when wrapping the model with `torch.nn.parallel.DistributedDataParallel`.
- The first NCCL collective times out:
```
Watchdog caught collective operation timeout:
WorkNCCL(SeqNum=1, OpType=ALLGATHER, NumelIn=1, NumelOut=4, Timeout(ms)=600000)
```
- The stack trace points to DDP’s parameter-shape verification:
```
_verify_param_shape_across_processes
torch/nn/parallel/distributed.py:1020
```
- There are also repeated c10d socket warnings before the NCCL timeout:
```
The client socket cannot be initialized to connect to [localhost]:...
errno: 97 - Address family not supported by protocol

The client socket cannot be initialized to connect to [gpu008.bamboo]:...
errno: 97 - Address family not supported by protocol
```
- The final DDP error reports inconsistent parameter counts across ranks, for example:
```
DDP expects same model across all ranks, but Rank 2 has 51 params,
while rank 0 has inconsistent 0 params.
```
Given that the same job works on other nodes, this looks more like a node-specific distributed runtime / networking / NCCL issue than a problem in the training code itself. The “inconsistent 0 params” message may be a consequence of the failed initial collective rather than a real model mismatch.
Could you please help diagnose the affected nodes? Possible things to check might include:
- node-local hostname / localhost resolution,
- IPv4 vs IPv6 configuration,
- NCCL interface selection,
- InfiniBand / RDMA status,
- GPU peer-to-peer / NVLink / PCIe topology,
- driver / CUDA / NCCL health,
- whether simple `nccl-tests` or a minimal PyTorch DDP all-reduce test fails on the same nodes.
I attached the relevant logs. The issue appears at the first NCCL `ALLGATHER` during DDP initialization and times out after 600 seconds.
Best,
Maciej
```
[...]
[rank2]:[E606 09:09:08.164985364 ProcessGroupNCCL.cpp:757] [Rank 2] Watchdog caught collective operation timeout: WorkNCCL(SeqNum=1, OpType=ALLGATHER, NumelIn=1, NumelOut=4, Timeout(ms)=600000) ran for 600004 milliseconds before timing out.
[rank3]:[E606 09:09:08.164985344 ProcessGroupNCCL.cpp:757] [Rank 3] Watchdog caught collective operation timeout: WorkNCCL(SeqNum=1, OpType=ALLGATHER, NumelIn=1, NumelOut=4, Timeout(ms)=600000) ran for 600004 milliseconds before timing out.
[rank1]:[E606 09:09:08.164985344 ProcessGroupNCCL.cpp:757] [Rank 1] Watchdog caught collective operation timeout: WorkNCCL(SeqNum=1, OpType=ALLGATHER, NumelIn=1, NumelOut=4, Timeout(ms)=600000) ran for 600004 milliseconds before timing out.
[rank0]:[E606 09:09:08.164985384 ProcessGroupNCCL.cpp:757] [Rank 0] Watchdog caught collective operation timeout: WorkNCCL(SeqNum=1, OpType=ALLGATHER, NumelIn=1, NumelOut=4, Timeout(ms)=600000) ran for 600004 milliseconds before timing out.
[rank0]:[E606 09:09:08.169193939 ProcessGroupNCCL.cpp:2395] [PG ID 0 PG GUID 0(default_pg) Rank 0]  failure detected by watchdog at work sequence id: 1 PG status: last enqueued work: 1, last completed work: -1
[rank3]:[E606 09:09:08.169193719 ProcessGroupNCCL.cpp:2395] [PG ID 0 PG GUID 0(default_pg) Rank 3]  failure detected by watchdog at work sequence id: 1 PG status: last enqueued work: 1, last completed work: -1
[rank1]:[E606 09:09:08.169193689 ProcessGroupNCCL.cpp:2395] [PG ID 0 PG GUID 0(default_pg) Rank 1]  failure detected by watchdog at work sequence id: 1 PG status: last enqueued work: 1, last completed work: -1
[rank2]:[E606 09:09:08.169332087 ProcessGroupNCCL.cpp:2395] [PG ID 0 PG GUID 0(default_pg) Rank 2]  failure detected by watchdog at work sequence id: 1 PG status: last enqueued work: 1, last completed work: -1
[rank1]:[E606 09:09:08.170090097 ProcessGroupNCCL.cpp:801] Stack trace of the failed collective:
#0 _verify_param_shape_across_processes from /opt/venv/lib/python3.13/site-packages/torch/distributed/utils.py:286
#1 __init__ from /opt/venv/lib/python3.13/site-packages/torch/nn/parallel/distributed.py:1020
#2 main from /workspace/experiments/scripts/train.py:692
#3 <module> from /workspace/experiments/scripts/train.py:1647

[rank3]:[E606 09:09:08.170103647 ProcessGroupNCCL.cpp:801] Stack trace of the failed collective:
#0 _verify_param_shape_across_processes from /opt/venv/lib/python3.13/site-packages/torch/distributed/utils.py:286
#1 __init__ from /opt/venv/lib/python3.13/site-packages/torch/nn/parallel/distributed.py:1020
#2 main from /workspace/experiments/scripts/train.py:692
#3 <module> from /workspace/experiments/scripts/train.py:1647

[rank0]:[E606 09:09:08.170147966 ProcessGroupNCCL.cpp:801] Stack trace of the failed collective:
#0 _verify_param_shape_across_processes from /opt/venv/lib/python3.13/site-packages/torch/distributed/utils.py:286
#1 __init__ from /opt/venv/lib/python3.13/site-packages/torch/nn/parallel/distributed.py:1020
#2 main from /workspace/experiments/scripts/train.py:692
#3 <module> from /workspace/experiments/scripts/train.py:1647

[rank1]:[E606 09:09:08.170157826 ProcessGroupNCCL.cpp:2732] [PG ID 0 PG GUID 0(default_pg) Rank 1] First PG on this rank to signal dumping.
[rank0]:[E606 09:09:08.170172106 ProcessGroupNCCL.cpp:2732] [PG ID 0 PG GUID 0(default_pg) Rank 0] First PG on this rank to signal dumping.
[rank3]:[E606 09:09:08.170204045 ProcessGroupNCCL.cpp:2732] [PG ID 0 PG GUID 0(default_pg) Rank 3] First PG on this rank to signal dumping.
[rank2]:[E606 09:09:08.170255835 ProcessGroupNCCL.cpp:801] Stack trace of the failed collective:
#0 _verify_param_shape_across_processes from /opt/venv/lib/python3.13/site-packages/torch/distributed/utils.py:286
#1 __init__ from /opt/venv/lib/python3.13/site-packages/torch/nn/parallel/distributed.py:1020
#2 main from /workspace/experiments/scripts/train.py:692
#3 <module> from /workspace/experiments/scripts/train.py:1647

[rank2]:[E606 09:09:08.170276834 ProcessGroupNCCL.cpp:2732] [PG ID 0 PG GUID 0(default_pg) Rank 2] First PG on this rank to signal dumping.
[rank3]:[E606 09:09:08.172519935 ProcessGroupNCCL.cpp:1987] [PG ID 0 PG GUID 0(default_pg) Rank 3] Received a dump signal due to a collective timeout from this local rank and we will try our best to dump the debug info. Last enqueued NCCL work: 1, last completed NCCL work: -1.This is most likely caused by incorrect usages of collectives, e.g., wrong sizes used across ranks, the order of collectives is not same for all ranks or the scheduled collective, for some reason, didn't run. Additionally, this can be caused by GIL deadlock or other reasons such as network errors or bugs in the communications library (e.g. NCCL), etc.
[rank1]:[E606 09:09:08.172524105 ProcessGroupNCCL.cpp:1987] [PG ID 0 PG GUID 0(default_pg) Rank 1] Received a dump signal due to a collective timeout from this local rank and we will try our best to dump the debug info. Last enqueued NCCL work: 1, last completed NCCL work: -1.This is most likely caused by incorrect usages of collectives, e.g., wrong sizes used across ranks, the order of collectives is not same for all ranks or the scheduled collective, for some reason, didn't run. Additionally, this can be caused by GIL deadlock or other reasons such as network errors or bugs in the communications library (e.g. NCCL), etc.
[rank2]:[E606 09:09:08.172543984 ProcessGroupNCCL.cpp:1987] [PG ID 0 PG GUID 0(default_pg) Rank 2] Received a dump signal due to a collective timeout from this local rank and we will try our best to dump the debug info. Last enqueued NCCL work: 1, last completed NCCL work: -1.This is most likely caused by incorrect usages of collectives, e.g., wrong sizes used across ranks, the order of collectives is not same for all ranks or the scheduled collective, for some reason, didn't run. Additionally, this can be caused by GIL deadlock or other reasons such as network errors or bugs in the communications library (e.g. NCCL), etc.
[rank0]:[E606 09:09:08.172551974 ProcessGroupNCCL.cpp:1987] [PG ID 0 PG GUID 0(default_pg) Rank 0] Received a dump signal due to a collective timeout from this local rank and we will try our best to dump the debug info. Last enqueued NCCL work: 1, last completed NCCL work: -1.This is most likely caused by incorrect usages of collectives, e.g., wrong sizes used across ranks, the order of collectives is not same for all ranks or the scheduled collective, for some reason, didn't run. Additionally, this can be caused by GIL deadlock or other reasons such as network errors or bugs in the communications library (e.g. NCCL), etc.
[rank3]:[E606 09:09:08.172565234 ProcessGroupNCCL.cpp:1700] [PG ID 0 PG GUID 0(default_pg) Rank 3] ProcessGroupNCCL preparing to dump debug info. Include stack trace: 1, only active collectives: 0
[rank1]:[E606 09:09:08.172565154 ProcessGroupNCCL.cpp:1700] [PG ID 0 PG GUID 0(default_pg) Rank 1] ProcessGroupNCCL preparing to dump debug info. Include stack trace: 1, only active collectives: 0
[rank2]:[E606 09:09:08.172583394 ProcessGroupNCCL.cpp:1700] [PG ID 0 PG GUID 0(default_pg) Rank 2] ProcessGroupNCCL preparing to dump debug info. Include stack trace: 1, only active collectives: 0
[rank0]:[E606 09:09:08.172592054 ProcessGroupNCCL.cpp:1700] [PG ID 0 PG GUID 0(default_pg) Rank 0] ProcessGroupNCCL preparing to dump debug info. Include stack trace: 1, only active collectives: 0
[rank2]: Traceback (most recent call last):
[rank2]:   File "/workspace/experiments/scripts/train.py", line 1647, in <module>
[rank2]:     main()
[rank2]:     ~~~~^^
[rank2]:   File "/workspace/experiments/scripts/train.py", line 692, in main
[rank2]:     model = DDP(
[rank2]:         model, device_ids=[ddp_local_rank] if device_type == "cuda" else None
[rank2]:     )
[rank2]:   File "/opt/venv/lib/python3.13/site-packages/torch/nn/parallel/distributed.py", line 1020, in __init__
[rank2]:     _verify_param_shape_across_processes(self.process_group, parameters)
[rank2]:     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
[rank2]:   File "/opt/venv/lib/python3.13/site-packages/torch/distributed/utils.py", line 286, in _verify_param_shape_across_processes
[rank2]:     return dist._verify_params_across_processes(process_group, tensors, logger)
[rank2]:            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
[rank2]: RuntimeError: DDP expects same model across all ranks, but Rank 2 has 51 params, while rank 0 has inconsistent 0 params.
[rank0]: Traceback (most recent call last):
[rank0]:   File "/workspace/experiments/scripts/train.py", line 1647, in <module>
[rank0]:     main()
[rank0]:     ~~~~^^
[rank0]:   File "/workspace/experiments/scripts/train.py", line 692, in main
[rank0]:     model = DDP(
[rank0]:         model, device_ids=[ddp_local_rank] if device_type == "cuda" else None
[rank0]:     )
[rank0]:   File "/opt/venv/lib/python3.13/site-packages/torch/nn/parallel/distributed.py", line 1020, in __init__
[rank0]:     _verify_param_shape_across_processes(self.process_group, parameters)
[rank0]:     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
[rank0]:   File "/opt/venv/lib/python3.13/site-packages/torch/distributed/utils.py", line 286, in _verify_param_shape_across_processes
[rank0]:     return dist._verify_params_across_processes(process_group, tensors, logger)
[rank0]:            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
[rank0]: RuntimeError: DDP expects same model across all ranks, but Rank 0 has 51 params, while rank 1 has inconsistent 0 params.
[rank1]: Traceback (most recent call last):
[rank1]:   File "/workspace/experiments/scripts/train.py", line 1647, in <module>
[rank1]:     main()
[rank1]:     ~~~~^^
[rank1]:   File "/workspace/experiments/scripts/train.py", line 692, in main
[rank1]:     model = DDP(
[rank1]:         model, device_ids=[ddp_local_rank] if device_type == "cuda" else None
[rank1]:     )
[rank1]:   File "/opt/venv/lib/python3.13/site-packages/torch/nn/parallel/distributed.py", line 1020, in __init__
[rank1]:     _verify_param_shape_across_processes(self.process_group, parameters)
[rank1]:     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
[rank1]:   File "/opt/venv/lib/python3.13/site-packages/torch/distributed/utils.py", line 286, in _verify_param_shape_across_processes
[rank1]:     return dist._verify_params_across_processes(process_group, tensors, logger)
[rank1]:            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
[rank1]: RuntimeError: DDP expects same model across all ranks, but Rank 1 has 51 params, while rank 0 has inconsistent 0 params.
[rank3]: Traceback (most recent call last):
[rank3]:   File "/workspace/experiments/scripts/train.py", line 1647, in <module>
[rank3]:     main()
[rank3]:     ~~~~^^
[rank3]:   File "/workspace/experiments/scripts/train.py", line 692, in main
[rank3]:     model = DDP(
[rank3]:         model, device_ids=[ddp_local_rank] if device_type == "cuda" else None
[rank3]:     )
[rank3]:   File "/opt/venv/lib/python3.13/site-packages/torch/nn/parallel/distributed.py", line 1020, in __init__
[rank3]:     _verify_param_shape_across_processes(self.process_group, parameters)
[rank3]:     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
[rank3]:   File "/opt/venv/lib/python3.13/site-packages/torch/distributed/utils.py", line 286, in _verify_param_shape_across_processes
[rank3]:     return dist._verify_params_across_processes(process_group, tensors, logger)
[rank3]:            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
[rank3]: RuntimeError: DDP expects same model across all ranks, but Rank 3 has 51 params, while rank 0 has inconsistent 0 params.
[rank3]:[E606 09:10:08.363069113 ProcessGroupNCCL.cpp:818] [Rank 3] Some NCCL operations have failed or timed out. Due to the asynchronous nature of CUDA kernels, subsequent GPU operations might run on corrupted/incomplete data.
[rank3]:[E606 09:10:08.363090292 ProcessGroupNCCL.cpp:832] [Rank 3] To avoid data inconsistency, we are taking the entire process down.
[rank1]:[E606 09:10:08.364365475 ProcessGroupNCCL.cpp:818] [Rank 1] Some NCCL operations have failed or timed out. Due to the asynchronous nature of CUDA kernels, subsequent GPU operations might run on corrupted/incomplete data.
[rank1]:[E606 09:10:08.364376295 ProcessGroupNCCL.cpp:832] [Rank 1] To avoid data inconsistency, we are taking the entire process down.
[rank1]:[E606 09:10:08.365393482 ProcessGroupNCCL.cpp:2197] [PG ID 0 PG GUID 0(default_pg) Rank 1] Process group watchdog thread terminated with exception: [Rank 1] Watchdog caught collective operation timeout: WorkNCCL(SeqNum=1, OpType=ALLGATHER, NumelIn=1, NumelOut=4, Timeout(ms)=600000) ran for 600004 milliseconds before timing out.
Exception raised from checkTimeout at /pytorch/torch/csrc/distributed/c10d/ProcessGroupNCCL.cpp:760 (most recent call first):
frame #0: c10::Error::Error(c10::SourceLocation, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >) + 0x9d (0x14e75056afad in /opt/venv/lib/python3.13/site-packages/torch/lib/libc10.so)
frame #1: c10d::ProcessGroupNCCL::WorkNCCL::checkTimeout(std::optional<std::chrono::duration<long, std::ratio<1l, 1000l> > >) + 0x2b1 (0x14e69c5e3611 in /opt/venv/lib/python3.13/site-packages/torch/lib/libtorch_cuda.so)
frame #2: c10d::ProcessGroupNCCL::Watchdog::runLoop() + 0x14f1 (0x14e69c5ef3b1 in /opt/venv/lib/python3.13/site-packages/torch/lib/libtorch_cuda.so)
frame #3: c10d::ProcessGroupNCCL::Watchdog::run() + 0x157 (0x14e69c5f0847 in /opt/venv/lib/python3.13/site-packages/torch/lib/libtorch_cuda.so)
frame #4: <unknown function> + 0xdc253 (0x14e75398c253 in /lib/x86_64-linux-gnu/libstdc++.so.6)
frame #5: <unknown function> + 0x94ac3 (0x14e756266ac3 in /lib/x86_64-linux-gnu/libc.so.6)
frame #6: <unknown function> + 0x126850 (0x14e7562f8850 in /lib/x86_64-linux-gnu/libc.so.6)

[rank3]:[E606 09:10:08.365403032 ProcessGroupNCCL.cpp:2197] [PG ID 0 PG GUID 0(default_pg) Rank 3] Process group watchdog thread terminated with exception: [Rank 3] Watchdog caught collective operation timeout: WorkNCCL(SeqNum=1, OpType=ALLGATHER, NumelIn=1, NumelOut=4, Timeout(ms)=600000) ran for 600004 milliseconds before timing out.
Exception raised from checkTimeout at /pytorch/torch/csrc/distributed/c10d/ProcessGroupNCCL.cpp:760 (most recent call first):
frame #0: c10::Error::Error(c10::SourceLocation, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >) + 0x9d (0x151bb3b6afad in /opt/venv/lib/python3.13/site-packages/torch/lib/libc10.so)
frame #1: c10d::ProcessGroupNCCL::WorkNCCL::checkTimeout(std::optional<std::chrono::duration<long, std::ratio<1l, 1000l> > >) + 0x2b1 (0x151affbe3611 in /opt/venv/lib/python3.13/site-packages/torch/lib/libtorch_cuda.so)
frame #2: c10d::ProcessGroupNCCL::Watchdog::runLoop() + 0x14f1 (0x151affbef3b1 in /opt/venv/lib/python3.13/site-packages/torch/lib/libtorch_cuda.so)
frame #3: c10d::ProcessGroupNCCL::Watchdog::run() + 0x157 (0x151affbf0847 in /opt/venv/lib/python3.13/site-packages/torch/lib/libtorch_cuda.so)
frame #4: <unknown function> + 0xdc253 (0x151bb6fdb253 in /lib/x86_64-linux-gnu/libstdc++.so.6)
frame #5: <unknown function> + 0x94ac3 (0x151bb98b5ac3 in /lib/x86_64-linux-gnu/libc.so.6)
frame #6: <unknown function> + 0x126850 (0x151bb9947850 in /lib/x86_64-linux-gnu/libc.so.6)

terminate called after throwing an instance of 'c10::DistBackendError'
terminate called after throwing an instance of 'c10::DistBackendError'
  what():  [PG ID 0 PG GUID 0(default_pg) Rank 1] Process group watchdog thread terminated with exception: [Rank 1] Watchdog caught collective operation timeout: WorkNCCL(SeqNum=1, OpType=ALLGATHER, NumelIn=1, NumelOut=4, Timeout(ms)=600000) ran for 600004 milliseconds before timing out.
Exception raised from checkTimeout at /pytorch/torch/csrc/distributed/c10d/ProcessGroupNCCL.cpp:760 (most recent call first):
frame #0: c10::Error::Error(c10::SourceLocation, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >) + 0x9d (0x14e75056afad in /opt/venv/lib/python3.13/site-packages/torch/lib/libc10.so)
frame #1: c10d::ProcessGroupNCCL::WorkNCCL::checkTimeout(std::optional<std::chrono::duration<long, std::ratio<1l, 1000l> > >) + 0x2b1 (0x14e69c5e3611 in /opt/venv/lib/python3.13/site-packages/torch/lib/libtorch_cuda.so)
frame #2: c10d::ProcessGroupNCCL::Watchdog::runLoop() + 0x14f1 (0x14e69c5ef3b1 in /opt/venv/lib/python3.13/site-packages/torch/lib/libtorch_cuda.so)
frame #3: c10d::ProcessGroupNCCL::Watchdog::run() + 0x157 (0x14e69c5f0847 in /opt/venv/lib/python3.13/site-packages/torch/lib/libtorch_cuda.so)
frame #4: <unknown function> + 0xdc253 (0x14e75398c253 in /lib/x86_64-linux-gnu/libstdc++.so.6)
frame #5: <unknown function> + 0x94ac3 (0x14e756266ac3 in /lib/x86_64-linux-gnu/libc.so.6)
frame #6: <unknown function> + 0x126850 (0x14e7562f8850 in /lib/x86_64-linux-gnu/libc.so.6)

Exception raised from run at /pytorch/torch/csrc/distributed/c10d/ProcessGroupNCCL.cpp:2203 (most recent call first):
frame #0: c10::Error::Error(c10::SourceLocation, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >) + 0x9d (0x14e75056afad in /opt/venv/lib/python3.13/site-packages/torch/lib/libc10.so)
frame #1: <unknown function> + 0x6915eb (0x14e69bd6b5eb in /opt/venv/lib/python3.13/site-packages/torch/lib/libtorch_cuda.so)
frame #2: <unknown function> + 0xdc253 (0x14e75398c253 in /lib/x86_64-linux-gnu/libstdc++.so.6)
frame #3: <unknown function> + 0x94ac3 (0x14e756266ac3 in /lib/x86_64-linux-gnu/libc.so.6)
frame #4: <unknown function> + 0x126850 (0x14e7562f8850 in /lib/x86_64-linux-gnu/libc.so.6)

  what():  [PG ID 0 PG GUID 0(default_pg) Rank 3] Process group watchdog thread terminated with exception: [Rank 3] Watchdog caught collective operation timeout: WorkNCCL(SeqNum=1, OpType=ALLGATHER, NumelIn=1, NumelOut=4, Timeout(ms)=600000) ran for 600004 milliseconds before timing out.
Exception raised from checkTimeout at /pytorch/torch/csrc/distributed/c10d/ProcessGroupNCCL.cpp:760 (most recent call first):
frame #0: c10::Error::Error(c10::SourceLocation, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >) + 0x9d (0x151bb3b6afad in /opt/venv/lib/python3.13/site-packages/torch/lib/libc10.so)
frame #1: c10d::ProcessGroupNCCL::WorkNCCL::checkTimeout(std::optional<std::chrono::duration<long, std::ratio<1l, 1000l> > >) + 0x2b1 (0x151affbe3611 in /opt/venv/lib/python3.13/site-packages/torch/lib/libtorch_cuda.so)
frame #2: c10d::ProcessGroupNCCL::Watchdog::runLoop() + 0x14f1 (0x151affbef3b1 in /opt/venv/lib/python3.13/site-packages/torch/lib/libtorch_cuda.so)
frame #3: c10d::ProcessGroupNCCL::Watchdog::run() + 0x157 (0x151affbf0847 in /opt/venv/lib/python3.13/site-packages/torch/lib/libtorch_cuda.so)
frame #4: <unknown function> + 0xdc253 (0x151bb6fdb253 in /lib/x86_64-linux-gnu/libstdc++.so.6)
frame #5: <unknown function> + 0x94ac3 (0x151bb98b5ac3 in /lib/x86_64-linux-gnu/libc.so.6)
frame #6: <unknown function> + 0x126850 (0x151bb9947850 in /lib/x86_64-linux-gnu/libc.so.6)

Exception raised from run at /pytorch/torch/csrc/distributed/c10d/ProcessGroupNCCL.cpp:2203 (most recent call first):
frame #0: c10::Error::Error(c10::SourceLocation, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >) + 0x9d (0x151bb3b6afad in /opt/venv/lib/python3.13/site-packages/torch/lib/libc10.so)
frame #1: <unknown function> + 0x6915eb (0x151aff36b5eb in /opt/venv/lib/python3.13/site-packages/torch/lib/libtorch_cuda.so)
frame #2: <unknown function> + 0xdc253 (0x151bb6fdb253 in /lib/x86_64-linux-gnu/libstdc++.so.6)
frame #3: <unknown function> + 0x94ac3 (0x151bb98b5ac3 in /lib/x86_64-linux-gnu/libc.so.6)
frame #4: <unknown function> + 0x126850 (0x151bb9947850 in /lib/x86_64-linux-gnu/libc.so.6)

[rank0]:[E606 09:10:09.678935251 ProcessGroupNCCL.cpp:818] [Rank 0] Some NCCL operations have failed or timed out. Due to the asynchronous nature of CUDA kernels, subsequent GPU operations might run on corrupted/incomplete data.
[rank2]:[E606 09:10:09.678934391 ProcessGroupNCCL.cpp:818] [Rank 2] Some NCCL operations have failed or timed out. Due to the asynchronous nature of CUDA kernels, subsequent GPU operations might run on corrupted/incomplete data.
[rank2]:[E606 09:10:09.678949651 ProcessGroupNCCL.cpp:832] [Rank 2] To avoid data inconsistency, we are taking the entire process down.
[rank0]:[E606 09:10:09.678952891 ProcessGroupNCCL.cpp:832] [Rank 0] To avoid data inconsistency, we are taking the entire process down.
[rank2]:[E606 09:10:09.679494943 ProcessGroupNCCL.cpp:2197] [PG ID 0 PG GUID 0(default_pg) Rank 2] Process group watchdog thread terminated with exception: [Rank 2] Watchdog caught collective operation timeout: WorkNCCL(SeqNum=1, OpType=ALLGATHER, NumelIn=1, NumelOut=4, Timeout(ms)=600000) ran for 600004 milliseconds before timing out.
Exception raised from checkTimeout at /pytorch/torch/csrc/distributed/c10d/ProcessGroupNCCL.cpp:760 (most recent call first):
frame #0: c10::Error::Error(c10::SourceLocation, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >) + 0x9d (0x1498aa16afad in /opt/venv/lib/python3.13/site-packages/torch/lib/libc10.so)
frame #1: c10d::ProcessGroupNCCL::WorkNCCL::checkTimeout(std::optional<std::chrono::duration<long, std::ratio<1l, 1000l> > >) + 0x2b1 (0x1497f61e3611 in /opt/venv/lib/python3.13/site-packages/torch/lib/libtorch_cuda.so)
frame #2: c10d::ProcessGroupNCCL::Watchdog::runLoop() + 0x14f1 (0x1497f61ef3b1 in /opt/venv/lib/python3.13/site-packages/torch/lib/libtorch_cuda.so)
frame #3: c10d::ProcessGroupNCCL::Watchdog::run() + 0x157 (0x1497f61f0847 in /opt/venv/lib/python3.13/site-packages/torch/lib/libtorch_cuda.so)
frame #4: <unknown function> + 0xdc253 (0x1498ad55f253 in /lib/x86_64-linux-gnu/libstdc++.so.6)
frame #5: <unknown function> + 0x94ac3 (0x1498afe39ac3 in /lib/x86_64-linux-gnu/libc.so.6)
frame #6: <unknown function> + 0x126850 (0x1498afecb850 in /lib/x86_64-linux-gnu/libc.so.6)

[rank0]:[E606 09:10:09.679495013 ProcessGroupNCCL.cpp:2197] [PG ID 0 PG GUID 0(default_pg) Rank 0] Process group watchdog thread terminated with exception: [Rank 0] Watchdog caught collective operation timeout: WorkNCCL(SeqNum=1, OpType=ALLGATHER, NumelIn=1, NumelOut=4, Timeout(ms)=600000) ran for 600004 milliseconds before timing out.
Exception raised from checkTimeout at /pytorch/torch/csrc/distributed/c10d/ProcessGroupNCCL.cpp:760 (most recent call first):
frame #0: c10::Error::Error(c10::SourceLocation, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >) + 0x9d (0x148aab36afad in /opt/venv/lib/python3.13/site-packages/torch/lib/libc10.so)
frame #1: c10d::ProcessGroupNCCL::WorkNCCL::checkTimeout(std::optional<std::chrono::duration<long, std::ratio<1l, 1000l> > >) + 0x2b1 (0x1489f73e3611 in /opt/venv/lib/python3.13/site-packages/torch/lib/libtorch_cuda.so)
frame #2: c10d::ProcessGroupNCCL::Watchdog::runLoop() + 0x14f1 (0x1489f73ef3b1 in /opt/venv/lib/python3.13/site-packages/torch/lib/libtorch_cuda.so)
frame #3: c10d::ProcessGroupNCCL::Watchdog::run() + 0x157 (0x1489f73f0847 in /opt/venv/lib/python3.13/site-packages/torch/lib/libtorch_cuda.so)
frame #4: <unknown function> + 0xdc253 (0x148aae7e0253 in /lib/x86_64-linux-gnu/libstdc++.so.6)
frame #5: <unknown function> + 0x94ac3 (0x148ab10baac3 in /lib/x86_64-linux-gnu/libc.so.6)
frame #6: <unknown function> + 0x126850 (0x148ab114c850 in /lib/x86_64-linux-gnu/libc.so.6)

terminate called after throwing an instance of 'terminate called after throwing an instance of 'c10::DistBackendError'
c10::DistBackendError'
  what():  [PG ID 0 PG GUID 0(default_pg) Rank 2] Process group watchdog thread terminated with exception: [Rank 2] Watchdog caught collective operation timeout: WorkNCCL(SeqNum=1, OpType=ALLGATHER, NumelIn=1, NumelOut=4, Timeout(ms)=600000) ran for 600004 milliseconds before timing out.
Exception raised from checkTimeout at /pytorch/torch/csrc/distributed/c10d/ProcessGroupNCCL.cpp:760 (most recent call first):
frame #0: c10::Error::Error(c10::SourceLocation, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >) + 0x9d (0x1498aa16afad in /opt/venv/lib/python3.13/site-packages/torch/lib/libc10.so)
frame #1: c10d::ProcessGroupNCCL::WorkNCCL::checkTimeout(std::optional<std::chrono::duration<long, std::ratio<1l, 1000l> > >) + 0x2b1 (0x1497f61e3611 in /opt/venv/lib/python3.13/site-packages/torch/lib/libtorch_cuda.so)
frame #2: c10d::ProcessGroupNCCL::Watchdog::runLoop() + 0x14f1 (0x1497f61ef3b1 in /opt/venv/lib/python3.13/site-packages/torch/lib/libtorch_cuda.so)
frame #3: c10d::ProcessGroupNCCL::Watchdog::run() + 0x157 (0x1497f61f0847 in /opt/venv/lib/python3.13/site-packages/torch/lib/libtorch_cuda.so)
frame #4: <unknown function> + 0xdc253 (0x1498ad55f253 in /lib/x86_64-linux-gnu/libstdc++.so.6)
frame #5: <unknown function> + 0x94ac3 (0x1498afe39ac3 in /lib/x86_64-linux-gnu/libc.so.6)
frame #6: <unknown function> + 0x126850 (0x1498afecb850 in /lib/x86_64-linux-gnu/libc.so.6)

Exception raised from run at /pytorch/torch/csrc/distributed/c10d/ProcessGroupNCCL.cpp:2203 (most recent call first):
frame #0: c10::Error::Error(c10::SourceLocation, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >) + 0x9d (0x1498aa16afad in /opt/venv/lib/python3.13/site-packages/torch/lib/libc10.so)
frame #1: <unknown function> + 0x6915eb (0x1497f596b5eb in /opt/venv/lib/python3.13/site-packages/torch/lib/libtorch_cuda.so)
frame #2: <unknown function> + 0xdc253 (0x1498ad55f253 in /lib/x86_64-linux-gnu/libstdc++.so.6)
frame #3: <unknown function> + 0x94ac3 (0x1498afe39ac3 in /lib/x86_64-linux-gnu/libc.so.6)
frame #4: <unknown function> + 0x126850 (0x1498afecb850 in /lib/x86_64-linux-gnu/libc.so.6)
  what():
[PG ID 0 PG GUID 0(default_pg) Rank 0] Process group watchdog thread terminated with exception: [Rank 0] Watchdog caught collective operation timeout: WorkNCCL(SeqNum=1, OpType=ALLGATHER, NumelIn=1, NumelOut=4, Timeout(ms)=600000) ran for 600004 milliseconds before timing out.
Exception raised from checkTimeout at /pytorch/torch/csrc/distributed/c10d/ProcessGroupNCCL.cpp:760 (most recent call first):
frame #0: c10::Error::Error(c10::SourceLocation, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >) + 0x9d (0x148aab36afad in /opt/venv/lib/python3.13/site-packages/torch/lib/libc10.so)
frame #1: c10d::ProcessGroupNCCL::WorkNCCL::checkTimeout(std::optional<std::chrono::duration<long, std::ratio<1l, 1000l> > >) + 0x2b1 (0x1489f73e3611 in /opt/venv/lib/python3.13/site-packages/torch/lib/libtorch_cuda.so)
frame #2: c10d::ProcessGroupNCCL::Watchdog::runLoop() + 0x14f1 (0x1489f73ef3b1 in /opt/venv/lib/python3.13/site-packages/torch/lib/libtorch_cuda.so)
frame #3: c10d::ProcessGroupNCCL::Watchdog::run() + 0x157 (0x1489f73f0847 in /opt/venv/lib/python3.13/site-packages/torch/lib/libtorch_cuda.so)
frame #4: <unknown function> + 0xdc253 (0x148aae7e0253 in /lib/x86_64-linux-gnu/libstdc++.so.6)
frame #5: <unknown function> + 0x94ac3 (0x148ab10baac3 in /lib/x86_64-linux-gnu/libc.so.6)
frame #6: <unknown function> + 0x126850 (0x148ab114c850 in /lib/x86_64-linux-gnu/libc.so.6)

Exception raised from run at /pytorch/torch/csrc/distributed/c10d/ProcessGroupNCCL.cpp:2203 (most recent call first):
frame #0: c10::Error::Error(c10::SourceLocation, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >) + 0x9d (0x148aab36afad in /opt/venv/lib/python3.13/site-packages/torch/lib/libc10.so)
frame #1: <unknown function> + 0x6915eb (0x1489f6b6b5eb in /opt/venv/lib/python3.13/site-packages/torch/lib/libtorch_cuda.so)
frame #2: <unknown function> + 0xdc253 (0x148aae7e0253 in /lib/x86_64-linux-gnu/libstdc++.so.6)
frame #3: <unknown function> + 0x94ac3 (0x148ab10baac3 in /lib/x86_64-linux-gnu/libc.so.6)
frame #4: <unknown function> + 0x126850 (0x148ab114c850 in /lib/x86_64-linux-gnu/libc.so.6)

W0606 09:10:14.761000 2888264 torch/distributed/elastic/multiprocessing/api.py:1014] Sending process 2888343 closing signal SIGTERM
W0606 09:10:14.761000 2888264 torch/distributed/elastic/multiprocessing/api.py:1014] Sending process 2888344 closing signal SIGTERM
W0606 09:10:14.761000 2888264 torch/distributed/elastic/multiprocessing/api.py:1014] Sending process 2888345 closing signal SIGTERM
E0606 09:10:14.926000 2888264 torch/distributed/elastic/multiprocessing/api.py:988] failed (exitcode: -6) local_rank: 3 (pid: 2888346) of binary: /opt/venv/bin/python
Traceback (most recent call last):
  File "/opt/venv/bin/torchrun", line 10, in <module>
    sys.exit(main())
             ~~~~^^
  File "/opt/venv/lib/python3.13/site-packages/torch/distributed/elastic/multiprocessing/errors/__init__.py", line 367, in wrapper
    return f(*args, **kwargs)
  File "/opt/venv/lib/python3.13/site-packages/torch/distributed/run.py", line 1016, in main
    run(args)
    ~~~^^^^^^
  File "/opt/venv/lib/python3.13/site-packages/torch/distributed/run.py", line 1007, in run
    elastic_launch(
    ~~~~~~~~~~~~~~~
        config=config,
        ~~~~~~~~~~~~~~
        entrypoint=cmd,
        ~~~~~~~~~~~~~~~
    )(*cmd_args)
    ~^^^^^^^^^^^
  File "/opt/venv/lib/python3.13/site-packages/torch/distributed/launcher/api.py", line 191, in __call__
    return launch_agent(self._config, self._entrypoint, list(args))
  File "/opt/venv/lib/python3.13/site-packages/torch/distributed/launcher/api.py", line 371, in launch_agent
    raise ChildFailedError(
    ...<2 lines>...
    )
[...]
```
