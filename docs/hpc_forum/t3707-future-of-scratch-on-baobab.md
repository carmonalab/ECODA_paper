# Future of scratch on Baobab

- Source: https://hpc-community.unige.ch/t/3707

- Created: 2024-10-31T10:46:52.150Z

- Posts: 8

- Category: 8

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yann.Sagon (2024-10-31T10:46:52.242Z)

Baobab unreachable[Baobab unreachable](https://hpc-community.unige.ch//hpc-community.unige.ch/t/baobab-unreachable/3702/7)
> It is a bit hard to work on the cluster these past days…
> Are you thinking about a long term solution for scratch being full ?
My answer is targeted to every user and not specifically @Ludovic.Dumoulin[@Ludovic.Dumoulin](https://hpc-community.unige.ch/u/ludovic.dumoulin) of course!
The issue is that we have already more than 3 x 1PB of storage for scratch. This has a high cost and we can’t/don’t want to increase its size.
We are already migrating hundreds of TBs to new storage bay to solve part of this short-term problem. We are also investing in a new storage that will be accessible from the three clusters, so users will not have to copy their data three times.
Until now, users have had the freedom to use scratch space as they see fit. And then the user has to remember to clean up their unneeded stuff, as scratch is only for short-term storage. For some unknown :wink: reason, many users forgot to remove their unused stuff and we are forced to remind them to do their duty. This is time-consuming for us and unfair to users who behave correctly.
So let us ask what changes would be more efficient/fair:
- We set an identical quota on the scratch storage for every user: not exactly the purpose of scratch space, a lot of administrative tasks to increase/decrease the quotas for users requiring more space for a good reason.
- We charge for the use of the scratch space: not very popular, at least people will probably remember to clean their stuff! More administrative tasks for us.
- We do the cleaning ourselves: any file older than 3 months is discarded. We expect a lot of complaints.
We welcome your feedbacks on this situation, and most importantly, ask yourself if you really need those old files cluttering up your scratch space! We have now 57TB free on Baobab, let’s try to have 500TB by the end of the week! You can do it! :index_pointing_at_the_viewer:


## Post 2 by @Baptiste.Leforestier (2024-10-31T12:25:20.687Z)

Hi all, my 2 cents on this.
I have previously worked on an HPC cluster in Grenoble, the solution there is automated and way more radical than the above suggestions. Scratch is scratch and it should never be long term storage anyway, so any file older than 30 days is automatically deleted, and users will receive an automated reminder email listing all the files that are > 20 days old (see public doc here[here](https://gricad-doc.univ-grenoble-alpes.fr/hpc/data_management/silenus/)) for them to deal within the 10 remaining days until deletion.
There were complaints at first because of the policy change, but then users organized their work differently and now it just works. In comparison to that, a 3 months delay, perhaps preceded by a reminder with the list of data to clean up for anything older than 1 or 2 months does not sound that bad at all and leaves ample time for users to make a move.
A lot of the complaints at the time were mostly because of the policy change happening abruptly and perhaps not enough communication with respect to alternative storage solutions.
I believe quotas and charging for scratch likely will also bring about complaints + the additional work load for the HPC team. Eventually if Baobab just works reliably the solution is worth it. Yggdrasil and Bamboo should probably not be left out of that, otherwise the easy and sneaky solution is to transfer the petabytes of files to their scratch and just displace the problem to another cluster.
Have a great day all,
Baptiste.


## Post 3 by @Ludovic.Dumoulin (2024-10-31T12:55:36.734Z)

Hello,
Thank you! Options 1 and 3 look great! While three months is quite a bit, I think it’s necessary for large job arrays that need to access dataframes saved on scratch or values from previous simulations (considering both queue and simulation time).
Any solution works for me as long as it helps make baobab reliable again.
Best regards and thank you again


## Post 4 by @Adrien.Albert (2024-11-04T15:54:13.263Z)

Hi All,
Here’s a section on robinhood that was rolled out some time ago. We had communicated about the 3-month policy, but only for information purposes. has the data been sorted with the report available? (Unfortunetly no)
For the time being, this tool has been deactivated (as we don’t apply a deletion policy). (but it’s on our todo list)
https://doc.eresearch.unige.ch/hpc/storage_on_hpc?s[]=robinhood#robinhood[https://doc.eresearch.unige.ch/hpc/storage_on_hpc?s[]=robinhood#robinhood](https://doc.eresearch.unige.ch/hpc/storage_on_hpc?s%5B%5D=robinhood#robinhood)
I agree with this point, scratch is scratch and shouldn’t be used as a long-term solution. But the fact is that there is always abuse, and on our side it takes a lot of time to contact and very often beg for data deletion or migration.
BUT it’s research data and we’re aware of its importance, so we only delete it in emergencies (or for specific reasons). So this point is very touchy.
I’d like to know how many users would be in favor of implementing a policy of deleting old scratch data in order to preserve performance and space. Even if such a policy would of course have its victims and be a brutal way of learning, but it would also encourage good practice.


## Post 5 by @Baptiste.Leforestier (2024-11-04T16:13:01.367Z)

Something to keep in mind, and I don’t know how much data that represents on scratch, but it’s likely a large chunk of it is just abandoned… Given how little effect the messages have on how much scratch is freed.
With a large portion of researchers using the cluster being PhD students or post-docs, i.e. with temporary contracts <4-5 years at most, it might just be a case of people forgetting to clean after they have moved on with their career elsewhere.
Do you think there is an easy way to know that? That means some of the data there may be quite, quite old and abandoned anyway.
Have a nice evening all,
Best,


## Post 6 by @Adrien.Albert (2024-11-04T16:37:25.873Z)

Hello Baptiste
You’re partly right, for accounts that are no longer in the active directory, we delete the account and the data (with PI notification). So we don’t have any “abandoned” data belonging to a former Unige member, but we definitely have data abandoned by a current Unige member.
The last Robinhood report (about 1 year ago) showed that we had at least 75% of data that was more than 3 months old.


## Post 7 by @maciej.falkiewicz (2024-11-07T14:29:47.445Z)

Yann.Sagon:
> We do the cleaning ourselves: any file older than 3 months is discarded. We expect a lot of complaints.
It depends on what you consider as “older than X”, will it be determined by change time (`ctime`) or access time (`atime`). If you respect `atime` then this sounds like the perfect solution.
In fact, there is no need to reinvent the wheel. Simply copy the policy of HPC Canada (Scratch purging policy - Alliance Doc[Scratch purging policy - Alliance Doc](https://docs.alliancecan.ca/wiki/Scratch_purging_policy)).


## Post 8 by @Adrien.Albert (2024-11-07T16:24:27.164Z)

Hi @maciej.falkiewicz[@maciej.falkiewicz](https://hpc-community.unige.ch/u/maciej.falkiewicz)
Thank you for your suggestion. However, each HPC site has its own unique requirements and policies.
For now, We will discuss, evaluate, and validate these different approach internally, with final decisions made in agreement with our hierarchie.
Best Regards,
