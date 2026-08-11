# Cannot connect to Yggdrasil after the maintenance

- Source: https://hpc-community.unige.ch/t/4019

- Created: 2025-07-31T09:51:15.985Z

- Tags: yggdrasil

- Posts: 13

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Lena.Parc (2025-07-31T09:51:16.029Z)

Hello,
I haven’t been able to connect to Yggdrasil since I received the e-mail confirming that maintenance had been completed.
I manage to log in, but then I get this message:
"-bash: fork: retry: Resource temporarily unavailable "
Thanks,
Léna


## Post 2 by @Max.Briel (2025-07-31T12:00:32.601Z)

Hey!
I’m having a similar issue with specifically using SSH to connect to a host within Visual Studio Code to yggdrasil. SSH directly via shell works fine.
I suspect this has to do with the new limit on concurrent processes on the login nodes.
VS code remote spawns many processes at the same time on startup (Too many server processes at once · Issue #10408 · microsoft/vscode-remote-release · GitHub[Too many server processes at once · Issue #10408 · microsoft/vscode-remote-release · GitHub](https://github.com/microsoft/vscode-remote-release/issues/10408) ; Note: the fix mentioned does not work when connecting to yggdrasil).
I understand the need for not running notebooks/code on the login node, but this breaks file exploration and connecting to an interactive node and JupyterLab session through VS code. I use this extensively and I think others too.
I hope this can be picked up quickly :slight_smile:
Thank you!
Cheers,
Max Briel


## Post 3 by @Lena.Parc (2025-07-31T12:14:53.298Z)

Hey,
I also experience this problem when connecting via VS Code, and occasionally with a simple terminal as well. I do not execute notebooks or code on the login node; I only use VS Code to develop and browse the files, since its exploration and editing features are very helpful.
I hope we can find a solution that works for everyone.
Thanks,
Léna


## Post 4 by @Adrien.Albert (2025-07-31T13:33:05.260Z)

Hi Lena,
Indeed, the issue comes from the new `MaxTasks` limit on the login node.
As I saw in your beautiful `.bashrc` :joy_cat: , you’re activating a Conda environment at each login or session. That might explain the issue. Note that the limit is applied per user, so the more sessions you have, the more processes are counted against the limit.
Here’s some documentation that might interest you on how to use Conda within a container (with cotainr):
doc.eresearch.unige.ch[doc.eresearch.unige.ch](https://doc.eresearch.unige.ch/hpc/applications_and_libraries#how_to_create_a_conda_environment_in_a_container)
### hpc:applications_and_libraries [eResearch Doc][hpc:applications_and_libraries [eResearch Doc]](https://doc.eresearch.unige.ch/hpc/applications_and_libraries#how_to_create_a_conda_environment_in_a_container)
For now, I’ve increased the limit, could you try again and let me know how it goes?
@Lena.Parc[@Lena.Parc](https://hpc-community.unige.ch/u/lena.parc) - PS: Everything comes at a price, get ready for the new season :woman_dancing:


## Post 5 by @Lena.Parc (2025-07-31T13:42:28.587Z)

Hi Adri,
It is working now ! I will have a look to use Conda within a container… oups
:woman_dancing:


## Post 6 by @Adrien.Albert (2025-07-31T13:54:34.503Z)

Hi @Max.Briel[@Max.Briel](https://hpc-community.unige.ch/u/max.briel)
You raise an important point. Unfortunately, not everyone understands the difference between a login node and a compute node. We’ve already had to deal with several users who were convinced that the login node was the cluster itself :sweat_smile:
We listen to all user feedback. Some request stricter resource limits, while others find such limitations problematic. It’s therefore challenging to find a configuration that satisfies everyone.
For now, we’re taking an iterative approach and will try to find a reasonable and acceptable balance. After the limit increase, is it working for you ?
I understand the need for VS Code. As another alternative, I suggest using OpenOnDemand, which allows you to launch an Interactive Web Session with VS Code directly on a compute node. OpenOnDemand also comes with a built-in file explorer.
You can try it here:
:backhand_index_pointing_right: https://openondemand.yggdrasil.hpc.unige.ch/[https://openondemand.yggdrasil.hpc.unige.ch/](https://openondemand.yggdrasil.hpc.unige.ch/)
Best :slight_smile:


## Post 7 by @Max.Briel (2025-07-31T14:29:49.083Z)

It’s currently working.
I fully get the stricter resource usage on the login node, and it makes sense :slight_smile:
Thank you!


## Post 8 by @Alexander.Froch (2025-07-31T15:31:34.325Z)

Hi @Adrien.Albert[@Adrien.Albert](https://hpc-community.unige.ch/u/Adrien.Albert)
Would it be possible to also increase the limit on baobab? I’m experiencing similar issues there ..
Edit: I was able to start VSCode on Yggdrasil and it connected, but there is the issue that VSCode cannot open a shell (Pty unresponsive). I guess this could also be related to the fact, that the max number is reached?


## Post 9 by @Max.Briel (2025-07-31T15:34:45.105Z)

It looks like it did not resolve it and issues also arise in normal terminals as mentioned by Lena.


## Post 10 by @Lena.Parc (2025-07-31T15:36:00.002Z)

Hi again,
It is not working right now…
I tried without success to change my conda env to use them within a container but I have an error when trying to build my .sif (cotainr build).
Cheers,
Léna


## Post 11 by @Adrien.Albert (2025-07-31T18:05:07.700Z)

Hello All,
## Yggdrasil
@Lena.Parc[@Lena.Parc](https://hpc-community.unige.ch/u/lena.parc) : no computing tasks on the login node, no compiling, building, or anything like that. :wink:
I wanted to restrict you all to just basic commands like `vim`, `emacs`, `nano`, `ls`, `rsync`, and `Slurm` on the login node, but my kind colleagues keep stopping me from being too harsh,  so consider yourself lucky! :grinning_face_with_smiling_eyes: (joke… or not I keep the mystery)
@Max.Briel[@Max.Briel](https://hpc-community.unige.ch/u/max.briel) You have vscode-server processes running on the login node and you source Conda in your `.bashrc`. The total number of processes you’re spawning exceeds the `TaskLimit`.
I’ve increased the limit again for now, but after some quick tests, the more SSH-remote VS Code sessions I open, the more tasks get created. Even when I terminate a VS Code session (exit command in the VS Code terminal), the processes on the login node stay alive in sleep mode. What a nightmare… Increasing the limit won’t solve the root cause, VS Code just keeps taking more and more. Such a child…
I think VS Code keeps those sleeping processes around to reconnect to the same session, maybe like `tmux`, but I’m not an expert. :thinking:
To fully clean my VS Code processes and avoid fork issues, I run:
- List VS Code processes:
```
pgrep -f vscode
```
- Kill all VS Code processes:
```
pkill -f vscode
```
## Baobab
The recent reboot of the Baobab login node triggered an update to the process limits configuration.
I manually removed the limitations on Baobab for now. I plan to fix this properly next Monday, but for now, I’m off for the weekend! :mending_heart:


## Post 12 by @Max.Briel (2025-08-05T07:38:19.482Z)

Thank you for the info! I hadn’t realized that VS Code starts a server on the login node when connecting to it. Thank you for pointing this out! Bad VS code. Looking at the processes on the login node, I think this might not be common knowledge :sweat_smile:
OpenOnDemand is definitely the easiest solution. Because I prefer my own workspace in VS Code, I preferred to connect manually to an interactive node through VS Code.
For others seeing this:
> I had to redo the steps here, because my previous ssh key was not working correctly when trying to connect to a node: hpc:access_the_hpc_clusters [eResearch Doc][hpc:access_the_hpc_clusters [eResearch Doc]](https://doc.eresearch.unige.ch/hpc/access_the_hpc_clusters#ssh_tunnel_and_socks_proxy)
> Additionally, I had to manually add the new key to `.ssh/authorized_keys`.
> @Adrien.Albert[@Adrien.Albert](https://hpc-community.unige.ch/u/Adrien.Albert) the `ssh-copy-id` command didn’t seem to add the key to the file? Though this might have been me using the command wrong.
> After setting up my local ssh config file, this was enough for me to `salloc` resources and directly connect to the allocated node afterwards from within VS code, proxyjumping though the login node :smiley:
@Adrien.Albert[@Adrien.Albert](https://hpc-community.unige.ch/u/Adrien.Albert) About conda: when following the conda instructions here hpc:applications_and_libraries [eResearch Doc][hpc:applications_and_libraries [eResearch Doc]](https://doc.eresearch.unige.ch/hpc/applications_and_libraries#conda) , the line `conda active <env_name>` errors due to conda wanting you to run `conda init` first. Running this adds the conda initialization to the `.bashrc`. How should I proceed here? (shall I make a new forum post for this?)


## Post 13 by @Florian.Brueck (2025-08-06T15:15:38.575Z)

Hello,
unfortunately I am facing the same issues as my colleagues above. However, I cannot even connect to Yggdrasil via the terminal anymore. Everything worked fine for me until yesterday. Any hint what I can do?
