# What about fast storage evolution on Baobab

- Source: https://hpc-community.unige.ch/t/3433

- Created: 2024-05-01T08:41:59.151Z

- Posts: 6

- Category: 8

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @maciej.falkiewicz (2024-05-01T08:41:59.151Z)

Dear @support[@support](https://hpc-community.unige.ch/groups/support) ,
`/srv/fast/users` run out of space. Are there any plans to add more volume?
Best regards,
Maciej


## Post 2 by @Yann.Sagon (2024-05-02T08:58:47.637Z)

Dear @maciej.falkiewicz[@maciej.falkiewicz](https://hpc-community.unige.ch/u/maciej.falkiewicz)
The fast storage on Baobab is still a proof of concept and there is no plan to increase significantly the capacity in a short term.
As our fast storage is getting full, we’ll enforce quota on fast storage on Baobab
- 500GB per user
- 1TB per group
As the purpose of this storage is to be ephemeral, we’ll as well erase the whole content at each Baobab maintenance. It is up to the user to copy the data back and forth.
Thanks for your understanding,
best
Yann


## Post 3 by @maciej.falkiewicz (2024-05-02T09:36:44.382Z)

Sounds reasonable.
Is it possible for a user to check how much fast storage his or her group is currently using?
Best,
Maciej


## Post 4 by @Yann.Sagon (2024-05-02T13:59:25.881Z)

We’ll enable a way to check the quota soon. In the meantime, you need to use a tool like `ncdu`.


## Post 5 by @Yann.Sagon (2024-07-10T08:34:50.543Z)

@maciej.falkiewicz[@maciej.falkiewicz](https://hpc-community.unige.ch/u/maciej.falkiewicz) did you had the opportunity to give a try with Bamboo? Here the storage is quite fast!


## Post 6 by @maciej.falkiewicz (2024-07-29T09:49:37.452Z)

@Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon) my experience so far has been positive. However, there is currently a small number of Bamboo users, I will wait for a final assessment until the real traffic appears.
