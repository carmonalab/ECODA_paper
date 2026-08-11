# [baobab] gpu022 GPU0 has a hardware issue

- Source: https://hpc-community.unige.ch/t/4076

- Created: 2025-09-01T16:50:44.128Z

- Posts: 7

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Ramon.CalvoGonzalez (2025-09-01T16:50:44.190Z)

## Primary informations
Username: calvogon
Cluster: baobab
## Description
The gpu:0 of gpu022 seems to have a hardware issue.
## Steps to Reproduce
With `nvidia-smi` it can be seen that the GPU 0 has an ERR! (in the power usage column).
```
±----------------------------------------------------------------------------------------+| NVIDIA-SMI 575.51.03              Driver Version: 575.51.03      CUDA Version: 12.9     ||-----------------------------------------±-----------------------±---------------------+| GPU  Name                 Persistence-M | Bus-Id          Disp.A | Volatile Uncorr. ECC || Fan  Temp   Perf          Pwr:Usage/Cap |           Memory-Usage | GPU-Util  Compute M. ||                                         |                        |               MIG M. ||=========================================+========================+======================||   0  NVIDIA A100-PCIE-40GB          On  |   00000000:01:00.0 Off |                    0 || N/A   39C    P0           ERR!  /  250W |       0MiB /  40960MiB |      0%      Default ||                                         |                        |             Disabled |±----------------------------------------±-----------------------±---------------------+|   1  NVIDIA A100-PCIE-40GB          On  |   00000000:22:00.0 Off |                    0 || N/A   33C    P0             39W /  250W |       0MiB /  40960MiB |      0%      Default ||                                         |                        |             Disabled |±----------------------------------------±-----------------------±---------------------+|   2  NVIDIA A100-PCIE-40GB          On  |   00000000:41:00.0 Off |                    0 || N/A   34C    P0             34W /  250W |       0MiB /  40960MiB |      0%      Default ||                                         |                        |             Disabled |±----------------------------------------±-----------------------±---------------------+|   3  NVIDIA A100-PCIE-40GB          On  |   00000000:61:00.0 Off |                    0 || N/A   34C    P0             37W /  250W |       0MiB /  40960MiB |      0%      Default ||                                         |                        |             Disabled |±----------------------------------------±-----------------------±---------------------+|   4  NVIDIA A100-PCIE-40GB          On  |   00000000:81:00.0 Off |                    0 || N/A   33C    P0             35W /  250W |       0MiB /  40960MiB |      0%      Default ||                                         |                        |             Disabled |±----------------------------------------±-----------------------±---------------------+|   5  NVIDIA A100-PCIE-40GB          On  |   00000000:A1:00.0 Off |                    0 || N/A   33C    P0             35W /  250W |       0MiB /  40960MiB |      0%      Default ||                                         |                        |             Disabled |±----------------------------------------±-----------------------±---------------------+|   6  NVIDIA A100-PCIE-40GB          On  |   00000000:C1:00.0 Off |                    0 || N/A   32C    P0             32W /  250W |       0MiB /  40960MiB |      0%      Default ||                                         |                        |             Disabled |±----------------------------------------±-----------------------±---------------------+|   7  NVIDIA A100-PCIE-40GB          On  |   00000000:E1:00.0 Off |                    0 || N/A   32C    P0             34W /  250W |       0MiB /  40960MiB |      0%      Default ||                                         |                        |             Disabled |±----------------------------------------±-----------------------±---------------------+

±----------------------------------------------------------------------------------------+| Processes:                                                                              ||  GPU   GI   CI              PID   Type   Process name                        GPU Memory ||        ID   ID                                                               Usage      ||=========================================================================================||  No running processes found                                                             |±----------------------------------------------------------------------------------------+
```
When I try to access the GPU, the CUDA driver throws an exception saying that the CUDA driver was unable to create a stream.
```
Python 3.13.5 (main, Jul  1 2025, 18:37:36) [Clang 20.1.4 ] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import jax
>>> jax.device_count()
2025-09-01 18:38:59.204949: W external/xla/xla/service/platform_util.cc:220] unable to create StreamExecutor for CUDA:0: CUDA error: : CUDA_ERROR_UNKNOWN: unknown error
Traceback (most recent call last):
  File "/home/users/c/calvogon/dino-jax/.venv/lib/python3.13/site-packages/jax/_src/xla_bridge.py", line 820, in backends
    backend = _init_backend(platform)
  File "/home/users/c/calvogon/dino-jax/.venv/lib/python3.13/site-packages/jax/_src/xla_bridge.py", line 904, in _init_backend
    backend = registration.factory()
  File "/home/users/c/calvogon/dino-jax/.venv/lib/python3.13/site-packages/jax/_src/xla_bridge.py", line 576, in factory
    return xla_client.make_c_api_client(plugin_name, updated_options, None)
           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/users/c/calvogon/dino-jax/.venv/lib/python3.13/site-packages/jaxlib/xla_client.py", line 156, in make_c_api_client
    return _xla.get_c_api_client(
           ~~~~~~~~~~~~~~~~~~~~~^
        plugin_name,
        ^^^^^^^^^^^^
    ...<2 lines>...
        transfer_server_factory,
        ^^^^^^^^^^^^^^^^^^^^^^^^
    )
    ^
jaxlib._jax.XlaRuntimeError: INTERNAL: no supported devices found for platform CUDA

During handling of the above exception, another exception occurred:

Traceback (most recent call last):
  File "<python-input-1>", line 1, in <module>
    jax.device_count()
    ~~~~~~~~~~~~~~~~^^
  File "/home/users/c/calvogon/dino-jax/.venv/lib/python3.13/site-packages/jax/_src/xla_bridge.py", line 983, in device_count
    return int(get_backend(backend).device_count())
               ~~~~~~~~~~~^^^^^^^^^
  File "/home/users/c/calvogon/dino-jax/.venv/lib/python3.13/site-packages/jax/_src/util.py", line 298, in wrapper
    return cached(config.trace_context() if trace_context_in_key else _ignore(),
                  *args, **kwargs)
  File "/home/users/c/calvogon/dino-jax/.venv/lib/python3.13/site-packages/jax/_src/util.py", line 292, in cached
    return f(*args, **kwargs)
  File "/home/users/c/calvogon/dino-jax/.venv/lib/python3.13/site-packages/jax/_src/xla_bridge.py", line 952, in get_backend
    return _get_backend_uncached(platform)
  File "/home/users/c/calvogon/dino-jax/.venv/lib/python3.13/site-packages/jax/_src/xla_bridge.py", line 931, in _get_backend_uncached
    bs = backends()
  File "/home/users/c/calvogon/dino-jax/.venv/lib/python3.13/site-packages/jax/_src/xla_bridge.py", line 836, in backends
    raise RuntimeError(err_msg)
RuntimeError: Unable to initialize backend 'cuda': INTERNAL: no supported devices found for platform CUDA (you may need to uninstall the failing plugin package, or set JAX_PLATFORMS=cpu to skip this backend.)
```


## Post 2 by @Ramon.CalvoGonzalez (2025-09-01T16:51:52.236Z)

There’s a problem with the formatting of the first message, so I’m attaching the output of `nvidia-smi` here instead.
Screenshot 2025-09-01 at 18.51.21
Screenshot 2025-09-01 at 18.51.21735×199 8.51 KB
[Screenshot 2025-09-01 at 18.51.21735×199 8.51 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/b00d90b81e27220a67c949fc03a1eb9e8dff2fe2.png)


## Post 3 by @Gael.Rossignol (2025-09-03T09:00:34.966Z)

Dear Ramon,
Thanks to report this issue, I will disable this GPU to allow usage of other GPU and I will contact support to replace the GPU.
Have a nice day,


## Post 4 by @Gael.Rossignol (2025-09-03T11:55:59.872Z)

Hello,
Server is now available with 7 gpu’s.
Best regards,


## Post 5 by @Ramon.CalvoGonzalez (2025-09-03T12:54:30.440Z)

Thank you! Is there a time estimate of how long will take for the replacement to arrive?


## Post 6 by @Ramon.CalvoGonzalez (2025-10-09T09:31:08.445Z)

Dear HPC team,
Is there any update on the replacement for this GPU? Thanks you.
Best regards,


## Post 7 by @Gael.Rossignol (2025-10-17T12:22:58.672Z)

Dear Ramon,
Sorry for misunderstanding, I didn’t explain well.
One of 8 gpu’s is not working fine but server is out of waranty because too old. Our provider can’t provide us a new one by paying because this kind of gpu is out of stock.
The unique way is to use other 7 remaining GPU’s.
Best regards,
