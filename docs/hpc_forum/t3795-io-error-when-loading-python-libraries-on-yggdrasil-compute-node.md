# IO Error when loading python libraries on Yggdrasil compute node

- Source: https://hpc-community.unige.ch/t/3795

- Created: 2025-01-23T11:08:28.232Z

- Tags: yggdrasil

- Posts: 3

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @daniel.forerosanchez (2025-01-23T11:08:28.308Z)

When submitting a job to a GPU compute node, I am getting IO Error when importing libraries. The specific error point changes each time but always IO error.
```
srun: job 37768678 queued and waiting for resources
srun: job 37768678 has been allocated resources
Traceback (most recent call last):
  File "/home/users/d/dforeros/projects/diff-alpt-paper/src/cosmo_unit_hmc.py", line 1, in <module>                                                         
    import arviz as az
  File "/home/users/d/dforeros/.conda/envs/numpyro/lib/python3.11/site-packages/arviz/__init__.py", line 34, in <module>                                    
    from .plots import *
  File "/home/users/d/dforeros/.conda/envs/numpyro/lib/python3.11/site-packages/arviz/plots/__init__.py", line 3, in <module>                               
    from .autocorrplot import plot_autocorr
  File "/home/users/d/dforeros/.conda/envs/numpyro/lib/python3.11/site-packages/arviz/plots/autocorrplot.py", line 8, in <module>                           
    from .plot_utils import default_grid, filter_plotters_list, get_plotting_function                                                                       
  File "/home/users/d/dforeros/.conda/envs/numpyro/lib/python3.11/site-packages/arviz/plots/plot_utils.py", line 11, in <module>                            
    from scipy.stats import mode, rankdata
  File "/home/users/d/dforeros/.conda/envs/numpyro/lib/python3.11/site-packages/scipy/stats/__init__.py", line 605, in <module>                             
    from ._stats_py import *
  File "/home/users/d/dforeros/.conda/envs/numpyro/lib/python3.11/site-packages/scipy/stats/_stats_py.py", line 45, in <module>                             
    from . import distributions
  File "/home/users/d/dforeros/.conda/envs/numpyro/lib/python3.11/site-packages/scipy/stats/distributions.py", line 8, in <module>                          
    from ._distn_infrastructure import (rv_discrete, rv_continuous, rv_frozen)  # noqa: F401                                                                
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^                                                                              
  File "/home/users/d/dforeros/.conda/envs/numpyro/lib/python3.11/site-packages/scipy/stats/_distn_infrastructure.py", line 26, in <module>                 
    from scipy import integrate
  File "<frozen importlib._bootstrap>", line 1229, in _handle_fromlist
  File "/home/users/d/dforeros/.conda/envs/numpyro/lib/python3.11/site-packages/scipy/__init__.py", line 134, in __getattr__                                
    return _importlib.import_module(f'scipy.{name}')
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/users/d/dforeros/.conda/envs/numpyro/lib/python3.11/importlib/__init__.py", line 126, in import_module                                        
    return _bootstrap._gcd_import(name[level:], package, level)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/users/d/dforeros/.conda/envs/numpyro/lib/python3.11/site-packages/scipy/integrate/__init__.py", line 94, in <module>                          
    from ._quadrature import *
  File "<frozen importlib._bootstrap>", line 1176, in _find_and_load
  File "<frozen importlib._bootstrap>", line 1138, in _find_and_load_unlocked
  File "<frozen importlib._bootstrap>", line 1078, in _find_spec
  File "<frozen importlib._bootstrap_external>", line 1507, in find_spec
  File "<frozen importlib._bootstrap_external>", line 1479, in _get_spec
  File "<frozen importlib._bootstrap_external>", line 1619, in find_spec
  File "<frozen importlib._bootstrap_external>", line 1662, in _fill_cache
OSError: [Errno 121] Remote I/O error: '/home/users/d/dforeros/.conda/envs/numpyro/lib/python3.11/site-packages/scipy/integrate' 
```
what did you try:
Loading libraries on the login node (no problem)
what didn’t work:
Running the same script I’ve been using the last couple of weeks on a gpu compute node like
```
(numpyro) (yggdrasil)-[dforeros@login1 src]$ srun -t 4:00:00 --mem-per-cpu=8G -c12 -p shared-gpu --gpus=1 python bias_halos_hmc_mycw.py -tracer lrg -ztarget 0.5 -ngrid 360 -lbox 2000 -save_mocks -seeds 1 2 3 4 -train
```
path to the relevant files (logs, sbatch script, etc):
Various, I have tried this
`/home/users/d/dforeros/projects/y1-holimocks/src`
and this
`/home/users/d/dforeros/projects/diff-alpt-paper`
Thanks in advance


## Post 2 by @daniel.forerosanchez (2025-01-23T13:36:08.301Z)

Hi again, st sure if it is related but now I can’t even log in. It hangs after the password has been introduced:
```
(base) daniel@daniel-nb:~$ ssh -vv dforeros@login1.yggdrasil.hpc.unige.ch
OpenSSH_9.6p1 Ubuntu-3ubuntu13.5, OpenSSL 3.0.13 30 Jan 2024
debug1: Reading configuration data /home/daniel/.ssh/config
debug1: /home/daniel/.ssh/config line 5: Applying options for login1.yggdrasil.hpc.unige.ch
debug1: Reading configuration data /etc/ssh/ssh_config
debug1: /etc/ssh/ssh_config line 19: include /etc/ssh/ssh_config.d/*.conf matched no files
debug1: /etc/ssh/ssh_config line 21: Applying options for *
debug2: resolving "login1.yggdrasil.hpc.unige.ch" port 22
debug1: Connecting to login1.yggdrasil.hpc.unige.ch [129.194.64.11] port 22.
debug1: Connection established.
debug1: identity file /home/daniel/.ssh/id_rsa type -1
debug1: identity file /home/daniel/.ssh/id_rsa-cert type -1
debug1: identity file /home/daniel/.ssh/id_ecdsa type -1
debug1: identity file /home/daniel/.ssh/id_ecdsa-cert type -1
debug1: identity file /home/daniel/.ssh/id_ecdsa_sk type -1
debug1: identity file /home/daniel/.ssh/id_ecdsa_sk-cert type -1
debug1: identity file /home/daniel/.ssh/id_ed25519 type 3
debug1: identity file /home/daniel/.ssh/id_ed25519-cert type -1
debug1: identity file /home/daniel/.ssh/id_ed25519_sk type -1
debug1: identity file /home/daniel/.ssh/id_ed25519_sk-cert type -1
debug1: identity file /home/daniel/.ssh/id_xmss type -1
debug1: identity file /home/daniel/.ssh/id_xmss-cert type -1
debug1: identity file /home/daniel/.ssh/id_dsa type -1
debug1: identity file /home/daniel/.ssh/id_dsa-cert type -1
debug1: Local version string SSH-2.0-OpenSSH_9.6p1 Ubuntu-3ubuntu13.5
debug1: Remote protocol version 2.0, remote software version OpenSSH_8.0
debug1: compat_banner: match: OpenSSH_8.0 pat OpenSSH* compat 0x04000000
debug2: fd 3 setting O_NONBLOCK
debug1: Authenticating to login1.yggdrasil.hpc.unige.ch:22 as 'dforeros'
debug1: load_hostkeys: fopen /home/daniel/.ssh/known_hosts2: No such file or directory
debug1: load_hostkeys: fopen /etc/ssh/ssh_known_hosts: No such file or directory
debug1: load_hostkeys: fopen /etc/ssh/ssh_known_hosts2: No such file or directory
debug1: SSH2_MSG_KEXINIT sent
debug1: SSH2_MSG_KEXINIT received
debug2: local client KEXINIT proposal
debug2: KEX algorithms: sntrup761x25519-sha512@openssh.com,curve25519-sha256,curve25519-sha256@libssh.org,ecdh-sha2-nistp256,ecdh-sha2-nistp384,ecdh-sha2-nistp521,diffie-hellman-group-exchange-sha256,diffie-hellman-group16-sha512,diffie-hellman-group18-sha512,diffie-hellman-group14-sha256,ext-info-c,kex-strict-c-v00@openssh.com
debug2: host key algorithms: rsa-sha2-512-cert-v01@openssh.com,rsa-sha2-256-cert-v01@openssh.com,rsa-sha2-512,rsa-sha2-256,ssh-ed25519-cert-v01@openssh.com,ecdsa-sha2-nistp256-cert-v01@openssh.com,ecdsa-sha2-nistp384-cert-v01@openssh.com,ecdsa-sha2-nistp521-cert-v01@openssh.com,sk-ssh-ed25519-cert-v01@openssh.com,sk-ecdsa-sha2-nistp256-cert-v01@openssh.com,ssh-ed25519,ecdsa-sha2-nistp256,ecdsa-sha2-nistp384,ecdsa-sha2-nistp521,sk-ssh-ed25519@openssh.com,sk-ecdsa-sha2-nistp256@openssh.com
debug2: ciphers ctos: chacha20-poly1305@openssh.com,aes128-ctr,aes192-ctr,aes256-ctr,aes128-gcm@openssh.com,aes256-gcm@openssh.com
debug2: ciphers stoc: chacha20-poly1305@openssh.com,aes128-ctr,aes192-ctr,aes256-ctr,aes128-gcm@openssh.com,aes256-gcm@openssh.com
debug2: MACs ctos: umac-64-etm@openssh.com,umac-128-etm@openssh.com,hmac-sha2-256-etm@openssh.com,hmac-sha2-512-etm@openssh.com,hmac-sha1-etm@openssh.com,umac-64@openssh.com,umac-128@openssh.com,hmac-sha2-256,hmac-sha2-512,hmac-sha1
debug2: MACs stoc: umac-64-etm@openssh.com,umac-128-etm@openssh.com,hmac-sha2-256-etm@openssh.com,hmac-sha2-512-etm@openssh.com,hmac-sha1-etm@openssh.com,umac-64@openssh.com,umac-128@openssh.com,hmac-sha2-256,hmac-sha2-512,hmac-sha1
debug2: compression ctos: none,zlib@openssh.com,zlib
debug2: compression stoc: none,zlib@openssh.com,zlib
debug2: languages ctos: 
debug2: languages stoc: 
debug2: first_kex_follows 0 
debug2: reserved 0 
debug2: peer server KEXINIT proposal
debug2: KEX algorithms: curve25519-sha256,curve25519-sha256@libssh.org,ecdh-sha2-nistp256,ecdh-sha2-nistp384,ecdh-sha2-nistp521,diffie-hellman-group-exchange-sha256,diffie-hellman-group14-sha256,diffie-hellman-group16-sha512,diffie-hellman-group18-sha512,diffie-hellman-group-exchange-sha1,diffie-hellman-group14-sha1,kex-strict-s-v00@openssh.com
debug2: host key algorithms: rsa-sha2-512,rsa-sha2-256,ssh-rsa
debug2: ciphers ctos: aes256-gcm@openssh.com,chacha20-poly1305@openssh.com,aes256-ctr,aes256-cbc,aes128-gcm@openssh.com,aes128-ctr,aes128-cbc
debug2: ciphers stoc: aes256-gcm@openssh.com,chacha20-poly1305@openssh.com,aes256-ctr,aes256-cbc,aes128-gcm@openssh.com,aes128-ctr,aes128-cbc
debug2: MACs ctos: hmac-sha2-256-etm@openssh.com,hmac-sha1-etm@openssh.com,umac-128-etm@openssh.com,hmac-sha2-512-etm@openssh.com,hmac-sha2-256,hmac-sha1,umac-128@openssh.com,hmac-sha2-512
debug2: MACs stoc: hmac-sha2-256-etm@openssh.com,hmac-sha1-etm@openssh.com,umac-128-etm@openssh.com,hmac-sha2-512-etm@openssh.com,hmac-sha2-256,hmac-sha1,umac-128@openssh.com,hmac-sha2-512
debug2: compression ctos: none,zlib@openssh.com
debug2: compression stoc: none,zlib@openssh.com
debug2: languages ctos: 
debug2: languages stoc: 
debug2: first_kex_follows 0 
debug2: reserved 0 
debug1: kex: algorithm: curve25519-sha256
debug1: kex: host key algorithm: rsa-sha2-512
debug1: kex: server->client cipher: chacha20-poly1305@openssh.com MAC: <implicit> compression: none
debug1: kex: client->server cipher: chacha20-poly1305@openssh.com MAC: <implicit> compression: none
debug1: expecting SSH2_MSG_KEX_ECDH_REPLY
debug1: SSH2_MSG_KEX_ECDH_REPLY received
debug1: Server host key: ssh-rsa SHA256:tKqp4nljL+EGVKl8T0VF2nS36DkHVFMpLxQOPg/gKvg
debug1: load_hostkeys: fopen /home/daniel/.ssh/known_hosts2: No such file or directory
debug1: load_hostkeys: fopen /etc/ssh/ssh_known_hosts: No such file or directory
debug1: load_hostkeys: fopen /etc/ssh/ssh_known_hosts2: No such file or directory
debug1: Host 'login1.yggdrasil.hpc.unige.ch' is known and matches the RSA host key.
debug1: Found key in /home/daniel/.ssh/known_hosts:3
debug1: ssh_packet_send2_wrapped: resetting send seqnr 3
debug2: ssh_set_newkeys: mode 1
debug1: rekey out after 134217728 blocks
debug1: SSH2_MSG_NEWKEYS sent
debug1: expecting SSH2_MSG_NEWKEYS
debug1: ssh_packet_read_poll2: resetting read seqnr 3
debug1: SSH2_MSG_NEWKEYS received
debug2: ssh_set_newkeys: mode 0
debug1: rekey in after 134217728 blocks
debug1: SSH2_MSG_EXT_INFO received
debug1: kex_ext_info_client_parse: server-sig-algs=<ssh-ed25519,ssh-rsa,rsa-sha2-256,rsa-sha2-512,ssh-dss,ecdsa-sha2-nistp256,ecdsa-sha2-nistp384,ecdsa-sha2-nistp521>
debug2: service_accept: ssh-userauth
debug1: SSH2_MSG_SERVICE_ACCEPT received
debug1: Authentications that can continue: publickey,password,keyboard-interactive,hostbased
debug1: Next authentication method: publickey
debug1: get_agent_identities: bound agent to hostkey
debug1: get_agent_identities: agent returned 2 keys
debug1: Will attempt key: /home/daniel/.ssh/id_ed25519 ED25519 SHA256:7Yvw4JGJzwQztu0jg/ZFQugtJ0p/SFwhzNV5vJyE45g agent
debug1: Will attempt key:  RSA SHA256:7ihm+UqhHIiSUQsbb4qO6WyIWtoslW3XBUEfMzScR6U agent
debug1: Will attempt key: /home/daniel/.ssh/id_rsa 
debug1: Will attempt key: /home/daniel/.ssh/id_ecdsa 
debug1: Will attempt key: /home/daniel/.ssh/id_ecdsa_sk 
debug1: Will attempt key: /home/daniel/.ssh/id_ed25519_sk 
debug1: Will attempt key: /home/daniel/.ssh/id_xmss 
debug1: Will attempt key: /home/daniel/.ssh/id_dsa 
debug2: pubkey_prepare: done
debug1: Offering public key: /home/daniel/.ssh/id_ed25519 ED25519 SHA256:7Yvw4JGJzwQztu0jg/ZFQugtJ0p/SFwhzNV5vJyE45g agent
debug2: we sent a publickey packet, wait for reply
debug1: Authentications that can continue: publickey,password,keyboard-interactive,hostbased
debug1: Offering public key:  RSA SHA256:7ihm+UqhHIiSUQsbb4qO6WyIWtoslW3XBUEfMzScR6U agent
debug2: we sent a publickey packet, wait for reply
debug1: Authentications that can continue: publickey,password,keyboard-interactive,hostbased
debug1: Trying private key: /home/daniel/.ssh/id_rsa
debug1: Trying private key: /home/daniel/.ssh/id_ecdsa
debug1: Trying private key: /home/daniel/.ssh/id_ecdsa_sk
debug1: Trying private key: /home/daniel/.ssh/id_ed25519_sk
debug1: Trying private key: /home/daniel/.ssh/id_xmss
debug1: Trying private key: /home/daniel/.ssh/id_dsa
debug2: we did not send a packet, disable method
debug1: Next authentication method: keyboard-interactive
debug2: userauth_kbdint
debug2: we sent a keyboard-interactive packet, wait for reply
debug2: input_userauth_info_req: entering
debug2: input_userauth_info_req: num_prompts 1
(dforeros@login1.yggdrasil.hpc.unige.ch) Password: 
debug2: input_userauth_info_req: entering
debug2: input_userauth_info_req: num_prompts 0
Authenticated to login1.yggdrasil.hpc.unige.ch ([129.194.64.11]:22) using "keyboard-interactive".
debug1: channel 0: new session [client-session] (inactive timeout: 0)
debug2: channel 0: send open
debug1: Requesting no-more-sessions@openssh.com
debug1: Entering interactive session.
debug1: pledge: filesystem
debug1: client_input_global_request: rtype hostkeys-00@openssh.com want_reply 0
debug1: client_input_hostkeys: searching /home/daniel/.ssh/known_hosts for login1.yggdrasil.hpc.unige.ch / (none)
debug1: client_input_hostkeys: searching /home/daniel/.ssh/known_hosts2 for login1.yggdrasil.hpc.unige.ch / (none)
debug1: client_input_hostkeys: hostkeys file /home/daniel/.ssh/known_hosts2 does not exist
debug1: client_input_hostkeys: no new or deprecated keys from server
debug2: channel_input_open_confirmation: channel 0: callback start
debug2: fd 3 setting TCP_NODELAY
debug2: client_session2_setup: id 0
debug2: channel 0: request pty-req confirm 1
debug1: Sending environment.
debug1: channel 0: setting env LANG = "en_US.UTF-8"
debug2: channel 0: request env confirm 0
debug2: channel 0: request shell confirm 1
debug1: pledge: fork
debug2: channel_input_open_confirmation: channel 0: callback done
debug2: channel 0: open confirm rwindow 0 rmax 32768
debug2: channel_input_status_confirm: type 99 id 0
debug2: PTY allocation request accepted on channel 0
debug2: channel 0: rcvd adjust 2097152
debug2: channel_input_status_confirm: type 99 id 0
debug2: shell request accepted on channel 0
debug2: client_check_window_change: changed
debug2: channel 0: request window-change confirm 0
debug2: client_check_window_change: changed
debug2: channel 0: request window-change confirm 0
```


## Post 3 by @Adrien.Albert (2025-01-24T09:24:56.334Z)

I/O remote error => duplicate issue Login Issues Yggdrasil - #10 by Adrien.Albert[Login Issues Yggdrasil - #10 by Adrien.Albert](https://hpc-community.unige.ch/t/login-issues-yggdrasil/3796/10)
Close thread


## Post 4 by @Adrien.Albert (2025-01-24T09:25:04.044Z)
