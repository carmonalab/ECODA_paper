# Poor GPU performance on Baobab gpu022

- Source: https://hpc-community.unige.ch/t/4006

- Created: 2025-07-16T13:01:34.727Z

- Tags: baobab, bamboo

- Posts: 6

- Category: 8

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Ramon.CalvoGonzalez (2025-07-16T13:01:34.825Z)

Dear HPC team,
I have been using gpu022 lately and found that it seemed to be running slower than it should.
To investigate this issue, I ran the same benchmark (training a small LLM) on a single GPU from baobab’s gpu022 and bamboo gpu003. You can find attached the charts that compares the time to run a training epoch (lower is better) as well as the MFU (Model FLOP utilization, higher is better). At the beginning they have the same performance but quickly baobab’s A100 slows down considerably.
The GPU temperatures have never been above 70 degrees, and thermal throttling was not active. I checked this with nvidia-smi.
Do you expect this behaviour? I really think baobab’s gpu022  should be running at the same speed as bamboo’s gpu003.
I also have an unrelated question: gpu003 from bamboo has a CPU with only 16 cores to handle x4 A100 80Gb. Are there any plans on upgrading the CPU? This CPU is not able to keep all the GPUs busy when training big machine learning models.
Thanks for your time.
Best,
Ramon.
comparison_plots
comparison_plots1200×500 36.9 KB
[comparison_plots1200×500 36.9 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/278d8cbb581481aa6e4cc3d58f220cc32e6da0ff.png)


## Post 2 by @Yann.Sagon (2025-07-16T13:07:42.790Z)

Dear @Ramon.CalvoGonzalez[@Ramon.CalvoGonzalez](https://hpc-community.unige.ch/u/ramon.calvogonzalez)
One major difference between the Baobab and Bamboo clusters is their home storage performance. Bamboo has full SSD storage, whereas Baobab has old HDD storage.
Could you try using the local[local](https://doc.eresearch.unige.ch/hpc/storage_on_hpc#scratch_directory_local_on_each_node) ‘/scratch’ directory of the compute node if your dataset fits in it, and let us know if it works better?
You can also use the fast storage[fast storage](https://doc.eresearch.unige.ch/hpc/storage_on_hpc#fast_directory) on Baobab, which has an SSD backend. Let us know if you would like to try this; we can create a space for you.
Even though the CPU on gpu003.bamboo has fewer cores, its frequency is higher than that on gpu22.baobab; this may also be a factor.
We didn’t plan to upgrade the CPU on gpu003, but we’re glad you let us know that this is an issue. What we can do is try enabling HT (two threads per core). We disabled it on the compute node as it wasn’t working well.
Best regards
Yann


## Post 3 by @Ramon.CalvoGonzalez (2025-07-16T13:17:13.869Z)

Dear @Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon),
Thanks for the quick reply.
The jobs I have tried were run using the /scratch storage of both nodes. The data is really small and fits in RAM, so the storage cannot be the problem.
The GPU utilization was constantly 100% on both nodes, as well as the power consumption, which was always close to the maximum reported by nvidia-smi. This rules out the weak CPU hypothesis.
I can see that the base frequency of gpu022 is half to that of gpu003. But I am unsure if this could be the problem since the frequency seems to be multiplied to around 1100-1200MHz.
I think enabling HT on gpu003 would be beneficial for tasks that are IO bound, but I doubt it would change much if there is heavy data pre-processing on the CPU.
Best,
Ramon.


## Post 4 by @Ramon.CalvoGonzalez (2025-07-22T17:38:55.795Z)

Dear @Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon),
I just wanted to follow-up on this issue, I’m wondering if there’s been any progress on investigating the problem or if you need any additional information from me.
This issue is particularly important to us since it’s affecting our group’s dedicated node and impacting our research workflow. I’d really appreciate any update you might have on the timeline for resolution.
Thanks for your time.
Best regards,
Ramon.


## Post 5 by @Yann.Sagon (2025-07-23T08:03:52.013Z)

Dear @Ramon.CalvoGonzalez[@Ramon.CalvoGonzalez](https://hpc-community.unige.ch/u/ramon.calvogonzalez) the difference in speed between Baobab and Bamboo GPU node is less than 2, which is still on the same order of magnitude, but I understand your concern.
Having looked into the differences further, I found the following:
- The A100 GPU on gpu022 has 40 GB of RAM, whereas the GPU on gpu03.bamboo has 80 GB. The latter has faster memory. According to Nvidia[Nvidia](https://www.nvidia.com/en-us/data-center/a100/), there may be a speed boost depending on the application running.
- On gpu022, there is a parameter named NPS (Numa Per Socket) that is set to 1, whereas on gpu003, I guess it is set to 4. I have set gpu022 to drain in order to change this parameter. This may improve the bandwidth between the memory, CPU and PCI. After making this change, you can rerun your application to see if it improves the speed.


## Post 6 by @Ramon.CalvoGonzalez (2025-07-23T19:00:01.841Z)

Dear @Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon),
Thank you for your reply.
> The A100 GPU on gpu022 has 40 GB of RAM, whereas the GPU on gpu03.bamboo has 80 GB. The latter has faster memory. According to Nvidia[Nvidia](https://www.nvidia.com/en-us/data-center/a100/), there may be a speed boost depending on the application running.
While it is true that the A100 80Gb has ~25% more memory bandwidth, which would translate in improved performance on memory-bound tasks, the task I used to test the performance is actually compute-bound. I am certain of this because I was monitoring the GPU power usage and CUDA core utilization, which both were constantly at 100%.
Under compute-bound tasks, both versions of the A100 should behave the same, since they are essentially the same chip. In the NVIDIA specs, the FLOPS are equivalent for both.
Moreover, for the first iterations, both GPUs have the same performance. Which would be impossible if the 80Gb was faster because it had strictly better hardware.
I believe there is some kind of throttling kicking in on the A100 40Gb PCIe cards. But I have not been able to identify the reason.
> On gpu022, there is a parameter named NPS (Numa Per Socket) that is set to 1, whereas on gpu003, I guess it is set to 4. I have set gpu022 to drain in order to change this parameter. This may improve the bandwidth between the memory, CPU and PCI. After making this change, you can rerun your application to see if it improves the speed.
Thanks. I am currently running a job that will take a long time on gpu022. But I will try to profile gpu022 once this parameter has been changed. Although I am sure that this will not change much, since as I said before, the GPUs are compute-bound. Moreover, the dataset I’m loading from RAM to the GPU is only 3Mb in total. So it is extremely improbable that the PCIe is the bottleneck.
I have run more experiments on other A100s from baobab, both 40Gb and 80Gb models (image attached below). The 80Gb models are able to sustain their performance under constant load, but all the 40Gb models start with the same performance of the 80Gb models, and then quickly drop to a lower perf. Interestingly, gpu031’s performance decay is way slower.
image
image1200×500 63.1 KB
[image1200×500 63.1 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/fe874bcd4e0029a0bcccd04f1ea41435cef2c447.png)
I have rented a private node from Vast.ai to compare its performance with our HPC. Vast.ai offers a containerized environment from private individuals, which usually translates in inconsistent performance. But the GPU I rented (blue line) is able to sustain a performance similar to the ones in Bamboo when other common resources are not being drained.
image
image1200×500 58.5 KB
[image1200×500 58.5 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/fd87b585bb8d3770402c32f08ed78bc41e7bce7b.png)
I hope you found this information useful. Please let me know if I can help in any way.
Best regards,
Ramon.
