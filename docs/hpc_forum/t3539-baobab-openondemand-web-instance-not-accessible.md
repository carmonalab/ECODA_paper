# Baobab OpenOnDemand web instance not accessible

- Source: https://hpc-community.unige.ch/t/3539

- Created: 2024-07-14T16:50:06.216Z

- Posts: 5

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Dongpeng.Pan (2024-07-14T16:50:06.275Z)

## Primary informations
Username: $pand
Cluster: $baobab
## Description
When I log on https://ondemand.baobab.hpc.unige.ch/ [https://ondemand.baobab.hpc.unige.ch/](https://ondemand.baobab.hpc.unige.ch/), there is a 403 Error. Yet the bamboo is accessible, but my primary data is on Boabab.
## Steps to Reproduce
Log on to https://openondemand.baobab.hpc.unige.ch[https://openondemand.baobab.hpc.unige.ch](https://openondemand.baobab.hpc.unige.ch/)
## Expected Result
Opening the OOD dashboard and applying for an Rstudio Server and lunch.
## Actual Result
Error 403


## Post 2 by @Adrien.Albert (2024-07-15T11:13:14.131Z)

Hi @Dongpeng.Pan[@Dongpeng.Pan](https://hpc-community.unige.ch/u/dongpeng.pan)
I am sorry I can not reproduce the error, I have access with public network and Unige Network. Are you using another network ?


## Post 3 by @Adrien.Albert (2024-07-15T11:14:12.372Z)

It seems the problem is solved:
```
pand     2556610 2556561  0 09:03 ?        00:00:54 Passenger RubyApp: /var/www/ood/apps/sys/dashboard (production)
```
Could you confirm me ?


## Post 4 by @Erica.Lastufka (2024-07-16T07:25:04.109Z)

I am experiencing this same error.
Username: lastufka
Cluster: baobab


## Post 5 by @Adrien.Albert (2024-07-16T12:27:54.105Z)

Hi @Erica.Lastufka[@Erica.Lastufka](https://hpc-community.unige.ch/u/erica.lastufka)
The problem appeared sporadically, but I was able to catch and fix it. I apologize for any inconvenience caused.
