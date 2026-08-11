# Syntax error ug_getNodeCharacteristicsSummary.py

- Source: https://hpc-community.unige.ch/t/4257

- Created: 2026-03-18T15:11:39.267Z

- Posts: 4

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @pablo.strasser1 (2026-03-18T15:11:39.324Z)

Bonjour,
Il semble y avoir une erreur dans le code python suivant.
L’erreur est propagé sur tous les noeuds du cluster.
(yggdrasil)-[strassep@login1 ~]$ ug_getNodeCharacteristicsSummary.py
File “/usr/local/sbin/ug_getNodeCharacteristicsSummary.py”, line 57
print(f"salut {cluster} {summary.get(“net_billing_per_year”, 0.0)}")
^
SyntaxError: f-string: unmatched ‘(’
(yggdrasil)-[strassep@login1 ~]$
Bonne journée.
Pablo


## Post 2 by @Yann.Sagon (2026-03-18T16:34:21.938Z)

Hello @pablo.strasser1[@pablo.strasser1](https://hpc-community.unige.ch/u/pablo.strasser1)
Yes I saw that: it should be fixed by the end of the week.
Best regards
Yann


## Post 3 by @Yann.Sagon (2026-03-19T14:27:00.565Z)

Dear @pablo.strasser1[@pablo.strasser1](https://hpc-community.unige.ch/u/pablo.strasser1)
the issue is now fixed. As extra bonus, it isn’t needed anymore to pre load modules. For this to work, you need to launch the script without the “.py” extension.
Ex:
- `ug_getNodeCharacteristicsSummary`
- `ug_slurm_usage_per_user`


## Post 4 by @pablo.strasser1 (2026-03-19T15:28:12.196Z)

Confirmed it work well.
Thanks
