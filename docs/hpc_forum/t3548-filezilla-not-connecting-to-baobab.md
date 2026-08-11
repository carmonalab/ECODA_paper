# FileZilla not connecting to Baobab

- Source: https://hpc-community.unige.ch/t/3548

- Created: 2024-07-18T09:33:52.563Z

- Posts: 4

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Maddalena.Bugatti (2024-07-18T09:33:52.607Z)

Good afternoon,
Since some days I cannot connect to Baobab using FileZilla.
I use these inputs:
Host: sftp://login2.baobab2.hpc.unige.ch
Username: Bugatti
Password: ***
Port: 22
The output I get are the following:
|Status:      |Connecting to login2.baobab2.hpc.unige.ch…|
|Response: |fzSftp started, protocol_version=11|
|Command: |open bugatti@login2.baobab2.hpc.unige.ch 22|
|Error:         |ssh_init: nodename nor servname provided, or not known|
|Error:         |Could not connect to server|
I know there is another ticket about a FileZilla problem very similar to mine, but I could not find a way to solve my issue.
Thank you very much in advance,
Best regards,
Maddalena Bugatti


## Post 2 by @Raphael.Rubino (2024-07-18T09:51:44.145Z)

Hello,
FTP works with baobab and the following host:
sftp://login1.baobab.hpc.unige.ch
Best regards


## Post 3 by @Maddalena.Bugatti (2024-07-18T09:54:22.163Z)

Gosh, I did not notice this written error!
Thank you very much and sorry for the stupid question.
Maddalena


## Post 4 by @Adrien.Albert (2024-07-18T15:20:16.832Z)

Thank you @Raphael.Rubino[@Raphael.Rubino](https://hpc-community.unige.ch/u/raphael.rubino) you beat me to it :+1:t3:


## Post 5 by @Adrien.Albert (2024-07-19T12:47:38.356Z)
