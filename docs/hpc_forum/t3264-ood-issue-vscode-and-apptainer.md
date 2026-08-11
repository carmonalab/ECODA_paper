# [OOD] issue vscode and apptainer

- Source: https://hpc-community.unige.ch/t/3264

- Created: 2024-01-23T13:44:11.895Z

- Posts: 4

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Adrien.Albert (2024-01-23T13:44:11.976Z)

Hi @John.Raine[@John.Raine](https://hpc-community.unige.ch/u/john.raine)
Following your comment: Can't delete .vscode-server - #15 by John.Raine[Can't delete .vscode-server - #15 by John.Raine](https://hpc-community.unige.ch/t/cant-delete-vscode-server/3255/15)
Could you give me more information about the context ?


## Post 2 by @John.Raine (2024-02-05T10:20:46.680Z)

Hi Adrien,
one of the big benefits of vscode is development inside an environment for debugging purposes, which requires running within the docker container which is eventually used for running the code.
Normally, this can be done with devcontainers natively in vs code, but over ssh remote this can only work natively with a docker or podman backend.
Unfortunately as docker isn’t an option and podman doesn’t seem to work in its place on the cluster with current permissions.
A workaround that has been used by many is to have the node.js server running vscode on the remote node to run within a container environment, and to connect to this environment.
With OOD there is no support for devcontainers to run the IDE in a docker environment for python development, with the node.js server itself running within apptainer with this plugin.
There has been discussion about this on the OOD discourse here[here](https://discourse.openondemand.org/t/vscode-showcase/2256).
But until there is some form of containerised workspaces that can run within an environment I can’t see this replacing the current workarounds (though these have also been unreliable of late).
If there is no environment support, there isn’t much benefit over using the VS code desktop client locally.
Cheers,
Johnny


## Post 3 by @Adrien.Albert (2024-02-09T12:56:51.711Z)

Hi @John.Raine[@John.Raine](https://hpc-community.unige.ch/u/john.raine)
John.Raine:
> Unfortunately as docker isn’t an option and podman doesn’t seem to work in its place on the cluster with current permissions.
Thank you for this explaination. What do you need to make it work ?


## Post 4 by @John.Raine (2024-02-09T16:27:50.088Z)

Hi Adrien,
not sure, it seems to have issues when building the containers within vs code, I imagine in rootless mode.
But I’ve not been able to pin down what the problem is to try and solve it.
Maybe you might have more luck from the log than I have had, I’ve put it in pastebin here[here](https://pastebin.com/jMCs7Yhh).
There aren’t many if any guides on setting up podman on a remote as a non-root user to follow, unfortunately.
Johnny
