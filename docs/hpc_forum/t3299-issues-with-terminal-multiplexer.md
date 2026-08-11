# Issues with terminal multiplexer

- Source: https://hpc-community.unige.ch/t/3299

- Created: 2024-02-09T12:48:09.411Z

- Posts: 1

- Category: 1

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Julian.Weninger (2024-02-09T12:48:09.473Z)

Hello,
I’m using the terminal multiplexer ‘tmux’ to handle multiple command lines and was experiencing issues recently.
Usually, tmux provides multiple command lines that are nicely arranged on the screen and can be returned to after logout and from different machines. I use every of the command lines according to the guidelines, submitting jobs using sbatch and following their progress.
In the last month, launching new tmux sessions would fail and I returned to single command lines. However, I was notified today that tmux was using the login’s CPU excessively. According to the log, “tmux new -s name” used 100% CPU over the last week.
Is there a reason to not use terminal multiplexer like ‘tmux’ on the login node? Can you explain why opening a new session would highjack the login node’s CPU?
Thank you in advance,
Julian
