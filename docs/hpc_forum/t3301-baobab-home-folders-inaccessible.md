# Baobab home folders inaccessible

- Source: https://hpc-community.unige.ch/t/3301

- Created: 2024-02-09T15:39:45.731Z

- Posts: 12

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Anton.Hanke (2024-02-09T15:39:45.795Z)

Hi,
Just now the Baobab user home directories became inaccessible for both me and other people in my group.
Upon login I get following output and am dropped in the login2 node root directory:
```
Could not chdir to home directory /home/users/h/hankea: Communication error on send
/usr/bin/xauth:  error in locking authority file /home/users/h/hankea/.Xauthority
X11 connection rejected because of wrong authentication.
-bash: /home/users/h/hankea/.bash_profile: Communication error on send
```
Could you please help?
Best wishes,
Anton


## Post 2 by @maciej.falkiewicz (2024-02-09T16:27:03.840Z)

I can confirm. Looks like it is now resolved but jobs were affected :frowning:


## Post 3 by @Raphael.Rubino (2024-02-10T08:11:24.990Z)

Baobab is down for me now.


## Post 4 by @Jonathan.Mutal (2024-02-10T09:49:32.490Z)

Baobab is down for me too.


## Post 5 by @Vera.Chau (2024-02-10T11:25:32.224Z)

Baobab is down for me as well.


## Post 6 by @Adrien.Albert (2024-02-11T11:32:06.211Z)

Hi all,
I am on it, I keep you inform !


## Post 7 by @Lorenzo.Bini (2024-02-11T11:51:07.633Z)

Baobab is down for me too


## Post 8 by @Adrien.Albert (2024-02-11T12:09:27.243Z)

Dear All,
I suspect the problem is a continuation of Friday’s issue and is due to a user running a process on the login node. Unfortunately, there are no logs that could help me determine the gulty one.
I will certainly remain extremely attentive to the process running on the login nodes on the next week and take the necessary steps to prevent other users from disrupting the system and to maintain a healthy working environment.
Best Regards


## Post 9 by @Cody.Cardenas (2024-02-11T15:29:10.800Z)

Hello Adrien,
I know jobs were affected, but I still have a job running (PID 7322051). Not sure if it is helpful for you at all.
Something else I noticed that may or may not be connected to this issue is that I cannot log in with VSCode, but can SSH in over the terminal. Its probably low priority given the current situation.
Thanks for tackling this issue over the weekend.


## Post 10 by @Anton.Hanke (2024-02-11T16:46:42.571Z)

Thank you for taking care of this over the weekend.
I am also able to login `login2.baobab` via ssh.
At login I get following output:
```
Error: Given mountpoint is invalid: /home

[BeeGFS Control Tool Version: 7.4.2
Refer to the default config file (/etc/beegfs/beegfs-client.conf)
or visit http://www.beegfs.com to find out about configuration options.]
```
Still I am dropped into `pwd = my/home/dir`, but all my data is gone, only my link to scratch exists :frowning:
Maybe this is a result of trying to fix things/ un/remounting of the home directories.
I just wanted to document here.
Thanks


## Post 11 by @Adrien.Albert (2024-02-11T23:19:24.571Z)

Hi all;
I was focus on the stucked login node… but there seems to be another problem with Home, I’m having a quick look but I’m not sure I can solve it yet, and I have to wait until Monday morning.


## Post 12 by @Adrien.Albert (2024-02-11T23:40:49.296Z)

Hi Guys,
The home directory is back, the login node seems to be OK with module Filesystem mounted (another discovered issue).
My first suspicion was wrong, the login node issue mystifies me, and has blinded me to the beegfs problem.
We apology for any inconvenience caused. :pray:t3:
