# Sandbox open-on-demand on bamboo

- Source: https://hpc-community.unige.ch/t/3872

- Created: 2025-03-18T15:59:21.976Z

- Tags: bamboo, openondemand

- Posts: 3

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Lucille.Delisle1 (2025-03-18T15:59:22.013Z)

Hi,
I tried to do a sandbox app on bamboo and I got this error message:
image
image1198×209 27.6 KB
[image1198×209 27.6 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/0b61e127706f594e4aba18ece1592ea57a25f851.png)
Do you have any idea?
Thanks
Lucille


## Post 2 by @Adrien.Albert (2025-03-18T23:10:05.806Z)

Hi Lucile,
Developing apps requires some knowledge of Open OnDemand. You can find the documentation here:
https://osc.github.io/ood-documentation/latest/how-tos/app-development/interactive.html[https://osc.github.io/ood-documentation/latest/how-tos/app-development/interactive.html](https://osc.github.io/ood-documentation/latest/how-tos/app-development/interactive.html)
My first suggestion is to check the `forms.yml` file. You will find:
```
cluster: baobab
```
You need to replace it with the name of your cluster (`bamboo` in your case).
I am also working on exporting our apps to a public GitLab repository to enable community contributions (pull requests).


## Post 3 by @Lucille.Delisle1 (2025-03-19T06:19:17.461Z)

Thanks Adrien,
I was not aware of this parameter. I fixed it and it works.
Thank you so much for the public GitLab! I look forward.
