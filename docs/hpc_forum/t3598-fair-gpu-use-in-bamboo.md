# "Fair" GPU use in Bamboo

- Source: https://hpc-community.unige.ch/t/3598

- Created: 2024-08-20T08:04:50.473Z

- Tags: bamboo

- Posts: 10

- Category: 8

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @daniel.forerosanchez (2024-08-20T08:04:50.517Z)

Hi, I have used bamboo for a couple of weeks and found the short queues quite a relief, and while I understand that it is a matter of time before more users migrate their research to bamboo and the queues fill up, I want to raise the issue of fair GPU usage, specially since it is quite easy for a single user (as it seems to be mostly the case at this very moment) to use all GPU nodes thus effectively blocking other people from doing any work.
My particular use case is: I am still developing and debugging my GPU code but I am unable to do a short (5min) test because other people are using all gpu nodes and have 40 or so jobs pending (as was the case yesterday). Shouldn’t there be a safeguard against this? Wouldn’t it be possible to at least ensure that one debug GPU will be available so the max wait time for a single test is 15 min?
I am aware there are fair use policies and limits on queued jobs (10k). It seems 10k pending jobs may be a good number for CPUs but not really for GPUs since it is far greater than the total number of GPUs available.
Thanks in advance for your help.


## Post 2 by @Yann.Sagon (2024-08-20T12:29:51.841Z)

Dear @daniel.forerosanchez[@daniel.forerosanchez](https://hpc-community.unige.ch/u/daniel.forerosanchez) thanks for your feedback.
In fact, gpu001 is available from two partitions: debug-gpu and shared-gpu, which is obviously a problem. As GPUs are quite expensive, it is not a good idea to use the whole GPU node for debugging purposes. What we’ll try to do in the next maintenance is to split the node in two, i.e. reserve 4 of the 8 GPUs for debugging purposes.
In the meantime, you can use debug-gpu on Yggdrasil which is only dedicated to debugging.
Best


## Post 3 by @daniel.forerosanchez (2024-09-26T12:03:35.772Z)

Hi, I am still unable to get some time on the debug gpu on Bamboo. I hoped it would be easier since the last maintenance but it doesn’t seem to be the case. Was the node split in the end? Which other solutions are available? Moving clusters this often is not really viable and makes it very hard to version control.
Thanks,


## Post 4 by @Raphael.Rubino (2024-09-26T13:46:57.016Z)

Hello,
Sorry, I am using a lot of resources on bamboo right now, I will kill my jobs.
About your job, are 32 CPUs necessary per GPU? lowering that number would increase your job priority, I assume.
Best regards


## Post 5 by @daniel.forerosanchez (2024-09-27T08:41:18.666Z)

Thanks, it seems to make little difference to ask for less cores in any case. Still unable to get into the debug node but it seems they are down or something now.


## Post 6 by @Yann.Sagon (2024-09-30T12:42:40.381Z)

Dear @daniel.forerosanchez[@daniel.forerosanchez](https://hpc-community.unige.ch/u/daniel.forerosanchez) sorry it is still in our long todo list :face_exhaling:


## Post 7 by @daniel.forerosanchez (2024-10-15T13:09:13.494Z)

Is it possible to do something else in the meantime? It is still quite impossible to use the GPUs most of the time for small jobs or development purposes.


## Post 8 by @daniel.forerosanchez (2025-01-29T10:02:11.368Z)

Hi, are there any advancements on this? There is currently a single user with 7+ hours obs on all gpus, making it impossible for the rest of us to even debug our code.


## Post 9 by @daniel.forerosanchez (2025-01-30T08:55:39.388Z)

Hi, @Raphael.Rubino[@Raphael.Rubino](https://hpc-community.unige.ch/u/raphael.rubino) would it be possible for you to not use all the GPUs of the cluster, all the time? Thanks!


## Post 10 by @Yann.Sagon (2025-01-30T13:13:49.470Z)

Dear @daniel.forerosanchez[@daniel.forerosanchez](https://hpc-community.unige.ch/u/daniel.forerosanchez)
We understand that this is a problem for you and other people. After considering the options, we have decided to make the following changes to improve the situation:
- remove the `debug-gpu` partition
- replace it with a new partition called `public-interactive-gpu` and limit each user to only one GPU par job and only one job at a time on this partition. This partition will have a maximum wall time of 4 hours.
- Move `gpu001.bamboo` from the `public-gpu` partition to this new partition.
Hopefully this will meet the needs of most users.
Best regards
