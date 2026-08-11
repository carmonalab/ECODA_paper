# Limited number of Matlab jobs on clusters

- Source: https://hpc-community.unige.ch/t/3542

- Created: 2024-07-16T13:23:41.513Z

- Posts: 4

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Nicolas.Cuny (2024-07-16T13:23:41.546Z)

Username: cuny
Cluster: yggdrasil, baobab
Description:
It seems that the number of job using matlab that can run in parallel on the clusters is limited. When I try to launch more than about 100 matlab jobs some jobs remain pending with license reason untill previous job end. Is there a way of getting rid of this limitation?


## Post 2 by @Adrien.Albert (2024-07-16T14:21:16.784Z)

Dear @Nicolas.Cuny[@Nicolas.Cuny](https://hpc-community.unige.ch/u/nicolas.cuny)
This doc should interest you
https://doc.eresearch.unige.ch/hpc/applications_and_libraries#compile_your_matlab_code[https://doc.eresearch.unige.ch/hpc/applications_and_libraries#compile_your_matlab_code](https://doc.eresearch.unige.ch/hpc/applications_and_libraries#compile_your_matlab_code)
> The idea of compiling Matlab code is to save on licenses usage. Indeed, once compiled, a Matlab code can be run without using any license.


## Post 3 by @Nicolas.Cuny (2024-07-17T11:20:19.676Z)

Thank you very much for your answer. I compiled my Matlab code and tried to run the resulting executable with the following command:
./run.sh args
but it appears that before the arguments of my function (args) I need to give to the executable the path where Matlab is installed on the cluster but I don’t find this path, could you tell me what is the path?
Thank you very much


## Post 4 by @Adrien.Albert (2024-07-17T13:19:25.214Z)

Hi @Nicolas.Cuny[@Nicolas.Cuny](https://hpc-community.unige.ch/u/nicolas.cuny)
The path depends of the version, but with the command which you can easily get it:
```
(baobab)-[alberta@login1 ~]$ ml MATLAB/2022a
(baobab)-[alberta@login1 ~]$ which matlab
/opt/ebsofts/MATLAB/2022a/bin/matlab
```
