# Gpu008 Faulty GPU 3

- Source: https://hpc-community.unige.ch/t/4364

- Created: 2026-08-10T11:28:41.757Z

- Tags: bamboo

- Posts: 4

- Category: 13

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Raphael.Rubino (2026-08-10T11:28:41.824Z)

Hello,
On bamboo, the node gpu008 has a faulty GPU (id 3).
Here is the error message encountered when trying to use it:
```
CUDA error: uncorrectable ECC error encountered
Search for `cudaErrorECCUncorrectable' in https://docs.nvidia.com/cuda/cuda-runtime-api/group__CUDART__TYPES.html for more information.
```
Thank you


## Post 2 by @Adrien.Albert (2026-08-11T09:01:56.670Z)

Hello @Raphael.Rubino[@Raphael.Rubino](https://hpc-community.unige.ch/u/raphael.rubino)
The node has been drained from production.
Could you share the steps to reproduce the issue so that we can identify the root cause more precisely?


## Post 3 by @Raphael.Rubino (2026-08-11T11:18:51.140Z)

Sure @Adrien.Albert[@Adrien.Albert](https://hpc-community.unige.ch/u/adrien.albert),
Please log on the node gpu008 with GPU #3 allocated (I have tested with single GPU instances only).
Then, in a python console (I have tested with python 3.10):
```
import torch
# declare a tensor
t=torch.zeros((2,2))
# print tensor
print(t)
# printed result
> tensor([[0., 0.], [0., 0.]])
# move the tensor to GPU
t=t.to("cuda")
```
The expected result should be:
```
# print the tensor once moved to GPU
print(t)
# printed result
> tensor([[0., 0.],
        [0., 0.]], device='cuda:0')
```
Instead, we get the error reported above in this thread.
Thank you.


## Post 4 by @Ghasem.Hajianfar (2026-08-11T12:31:23.898Z)

Hi,
I also got same error.
