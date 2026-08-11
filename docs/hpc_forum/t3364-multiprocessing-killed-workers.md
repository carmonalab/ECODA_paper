# Multiprocessing killed workers

- Source: https://hpc-community.unige.ch/t/3364

- Created: 2024-03-12T10:20:00.054Z

- Posts: 3

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yazdan.Salimi (2024-03-12T10:20:00.118Z)

Helllo everyone.
I am using nnunet for segmentation which uses multiprocessing during preprocessing, after few iterations of preprocessing it stops with this error.
I tried increasing the CPUS pertask to 12 and memory to 64 GB, still the issue happens.
I appreciate your help. the error is here:
```
  0%|          | 0/395 [00:00<?, ?it/s]
  0%|          | 1/395 [00:25<2:48:47, 25.71s/it]
  1%|          | 2/395 [00:26<1:12:25, 11.06s/it]
  1%|          | 3/395 [00:28<44:38,  6.83s/it]  
  1%|          | 4/395 [00:41<1:00:16,  9.25s/it]
  1%|â–         | 5/395 [00:41<38:40,  5.95s/it]  
  2%|â–         | 6/395 [00:41<26:34,  4.10s/it]
  2%|â–         | 7/395 [00:45<26:32,  4.10s/it]
  2%|â–         | 8/395 [00:58<44:14,  6.86s/it]
  2%|â–         | 9/395 [01:14<1:01:50,  9.61s/it]
  3%|â–Ž         | 10/395 [01:19<53:22,  8.32s/it] 
  3%|â–Ž         | 11/395 [01:29<56:25,  8.82s/it]
  3%|â–Ž         | 12/395 [01:31<41:53,  6.56s/it]
  3%|â–Ž         | 12/395 [01:36<51:19,  8.04s/it]
Traceback (most recent call last):
  File "/home/users/s/salimi/anaconda3/envs/nn/lib/python3.10/runpy.py", line 196, in _run_module_as_main
    return _run_code(code, main_globals, None,
  File "/home/users/s/salimi/anaconda3/envs/nn/lib/python3.10/runpy.py", line 86, in _run_code
    exec(code, run_globals)
  File "/home/users/s/salimi/anaconda3/envs/nn/lib/python3.10/site-packages/monai/apps/nnunet/__main__.py", line 20, in <module>
    fire.Fire({"nnUNetV2Runner": nnUNetV2Runner})
  File "/home/users/s/salimi/anaconda3/envs/nn/lib/python3.10/site-packages/fire/core.py", line 141, in Fire
    component_trace = _Fire(component, args, parsed_flag_args, context, name)
  File "/home/users/s/salimi/anaconda3/envs/nn/lib/python3.10/site-packages/fire/core.py", line 475, in _Fire
    component, remaining_args = _CallAndUpdateTrace(
  File "/home/users/s/salimi/anaconda3/envs/nn/lib/python3.10/site-packages/fire/core.py", line 691, in _CallAndUpdateTrace
    component = fn(*varargs, **kwargs)
  File "/home/users/s/salimi/anaconda3/envs/nn/lib/python3.10/site-packages/monai/apps/nnunet/nnunetv2_runner.py", line 956, in run
    self.plan_and_process()
  File "/home/users/s/salimi/anaconda3/envs/nn/lib/python3.10/site-packages/monai/apps/nnunet/nnunetv2_runner.py", line 486, in plan_and_process
    self.preprocess(c, n_proc, overwrite_plans_name, verbose)
  File "/home/users/s/salimi/anaconda3/envs/nn/lib/python3.10/site-packages/monai/apps/nnunet/nnunetv2_runner.py", line 406, in preprocess
    preprocess(
  File "/home/users/s/salimi/anaconda3/envs/nn/lib/python3.10/site-packages/nnunetv2/experiment_planning/plan_and_preprocess_api.py", line 142, in preprocess
    preprocess_dataset(d, plans_identifier, configurations, num_processes, verbose)
  File "/home/users/s/salimi/anaconda3/envs/nn/lib/python3.10/site-packages/nnunetv2/experiment_planning/plan_and_preprocess_api.py", line 121, in preprocess_dataset
    preprocessor.run(dataset_id, c, plans_identifier, num_processes=n)
  File "/home/users/s/salimi/anaconda3/envs/nn/lib/python3.10/site-packages/nnunetv2/preprocessing/preprocessors/default_preprocessor.py", line 246, in run
    raise RuntimeError('Some background worker is 6 feet under. Yuck. \n'
RuntimeError: Some background worker is 6 feet under. Yuck.
OK jokes aside.
One of your background processes is missing. This could be because of an error (look for an error message) or because it was killed by your OS due to running out of RAM. If you don't see an error message, out of RAM is likely the problem. In that case reducing the number of workers might help
slurmstepd: error: Detected 1 oom_kill event in StepId=32436650.0. Some of the step tasks have been OOM Killed.
srun: error: gpu006: task 0: Out Of Memory
```


## Post 2 by @Gael.Rossignol (2024-03-13T08:27:32.057Z)

Yazdan.Salimi:
> `32436650`
Dear Yazdan,
After checking your job on Yggdrasil cluster (please don’t forget to provide cluster name), I see the job used more than allocated memory :white_check_mark:
```
(yggdrasil)-[root@admin1 intel]$ seff 32436650
Job ID: 32436650
Cluster: yggdrasil
User/Group: salimi/hpc_users
State: FAILED (exit code 1)
Nodes: 1
Cores per node: 12
CPU Utilized: 00:12:42
CPU Efficiency: 27.73% of 00:45:48 core-walltime
Job Wall-clock time: 00:03:49
Memory Utilized: 43.09 GB
Memory Efficiency: 89.76% of 48.00 GB
```
In you output we can see the same information : “Detected 1 oom_kill event” ; oom is Out Of Memory.
Try to allocate more memory to the job can solve the issue.
Best regards,


## Post 3 by @Yazdan.Salimi (2024-04-03T11:50:57.751Z)

Thank you so much
selecting more ram solves it
