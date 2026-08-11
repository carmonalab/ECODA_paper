# Can't connect to Yggdrasil

- Source: https://hpc-community.unige.ch/t/3602

- Created: 2024-08-22T14:23:00.192Z

- Tags: yggdrasil

- Posts: 5

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @daniel.forerosanchez (2024-08-22T14:23:00.232Z)

Hi, I can’t connect to Yggdrasil, and there seems to be no message regarding any maintenance.
ssh log
```
(base) daniel@daniel-nb:~$ ssh -vv dforeros@login1.yggdrasil.hpc.unige.ch
OpenSSH_8.9p1 Ubuntu-3ubuntu0.10, OpenSSL 3.0.2 15 Mar 2022
debug1: Reading configuration data /home/daniel/.ssh/config
debug1: /home/daniel/.ssh/config line 45: Applying options for login1.yggdrasil.hpc.unige.ch
debug1: Reading configuration data /etc/ssh/ssh_config
debug1: /etc/ssh/ssh_config line 19: include /etc/ssh/ssh_config.d/*.conf matched no files
debug1: /etc/ssh/ssh_config line 21: Applying options for *
debug2: resolving "login1.yggdrasil.hpc.unige.ch" port 22
debug1: Connecting to login1.yggdrasil.hpc.unige.ch [129.194.64.11] port 22.
debug1: Connection established.
debug1: identity file /home/daniel/.ssh/id_rsa type 0
debug1: identity file /home/daniel/.ssh/id_rsa-cert type -1
debug1: identity file /home/daniel/.ssh/id_ecdsa type -1
debug1: identity file /home/daniel/.ssh/id_ecdsa-cert type -1
debug1: identity file /home/daniel/.ssh/id_ecdsa_sk type -1
debug1: identity file /home/daniel/.ssh/id_ecdsa_sk-cert type -1
debug1: identity file /home/daniel/.ssh/id_ed25519 type -1
debug1: identity file /home/daniel/.ssh/id_ed25519-cert type -1
debug1: identity file /home/daniel/.ssh/id_ed25519_sk type -1
debug1: identity file /home/daniel/.ssh/id_ed25519_sk-cert type -1
debug1: identity file /home/daniel/.ssh/id_xmss type -1
debug1: identity file /home/daniel/.ssh/id_xmss-cert type -1
debug1: identity file /home/daniel/.ssh/id_dsa type -1
debug1: identity file /home/daniel/.ssh/id_dsa-cert type -1
debug1: Local version string SSH-2.0-OpenSSH_8.9p1 Ubuntu-3ubuntu0.10
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
debug2: KEX algorithms: curve25519-sha256,curve25519-sha256@libssh.org,ecdh-sha2-nistp256,ecdh-sha2-nistp384,ecdh-sha2-nistp521,sntrup761x25519-sha512@openssh.com,diffie-hellman-group-exchange-sha256,diffie-hellman-group16-sha512,diffie-hellman-group18-sha512,diffie-hellman-group14-sha256,ext-info-c,kex-strict-c-v00@openssh.com
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
debug1: Found key in /home/daniel/.ssh/known_hosts:29
debug1: ssh_packet_send2_wrapped: resetting send seqnr 3
debug2: ssh_set_newkeys: mode 1
debug1: rekey out after 134217728 blocks
debug1: SSH2_MSG_NEWKEYS sent
debug1: expecting SSH2_MSG_NEWKEYS
debug1: ssh_packet_read_poll2: resetting read seqnr 3
debug1: SSH2_MSG_NEWKEYS received
debug2: ssh_set_newkeys: mode 0
debug1: rekey in after 134217728 blocks
debug1: get_agent_identities: bound agent to hostkey
debug1: get_agent_identities: agent returned 7 keys
debug1: Will attempt key: /home/daniel/.ssh/id_rsa RSA SHA256:A6pKrXedkwxgG0XwUHtzmeEsl4Qa3WJ3+udkTrFsXOo agent
debug1: Will attempt key: /home/daniel/.ssh/nersc RSA SHA256:f8qX3Fjz8wLTMaOMF9ZvLgMnQrkWFJaNqYreuJcStWU agent
debug1: Will attempt key: /home/daniel/.ssh/nersc RSA-CERT SHA256:f8qX3Fjz8wLTMaOMF9ZvLgMnQrkWFJaNqYreuJcStWU agent
debug1: Will attempt key: /home/daniel/.ssh/nersc RSA SHA256:ufWUf/5qmlDD+yeBN1YJuZZhLiHwZYeKB3L5IXG39tE agent
debug1: Will attempt key: /home/daniel/.ssh/nersc RSA-CERT SHA256:ufWUf/5qmlDD+yeBN1YJuZZhLiHwZYeKB3L5IXG39tE agent
debug1: Will attempt key: /home/daniel/.ssh/nersc RSA SHA256:9TiSKsh37060I3mt8RBlkpGMj7AL3byIl197AHUAQu0 agent
debug1: Will attempt key: /home/daniel/.ssh/nersc RSA-CERT SHA256:9TiSKsh37060I3mt8RBlkpGMj7AL3byIl197AHUAQu0 agent
debug1: Will attempt key: /home/daniel/.ssh/id_ecdsa 
debug1: Will attempt key: /home/daniel/.ssh/id_ecdsa_sk 
debug1: Will attempt key: /home/daniel/.ssh/id_ed25519 
debug1: Will attempt key: /home/daniel/.ssh/id_ed25519_sk 
debug1: Will attempt key: /home/daniel/.ssh/id_xmss 
debug1: Will attempt key: /home/daniel/.ssh/id_dsa 
debug2: pubkey_prepare: done
debug1: SSH2_MSG_EXT_INFO received
debug1: kex_input_ext_info: server-sig-algs=<ssh-ed25519,ssh-rsa,rsa-sha2-256,rsa-sha2-512,ssh-dss,ecdsa-sha2-nistp256,ecdsa-sha2-nistp384,ecdsa-sha2-nistp521>
debug2: service_accept: ssh-userauth
debug1: SSH2_MSG_SERVICE_ACCEPT received
debug1: Authentications that can continue: publickey,password,keyboard-interactive,hostbased
debug1: Next authentication method: publickey
debug1: Offering public key: /home/daniel/.ssh/id_rsa RSA SHA256:A6pKrXedkwxgG0XwUHtzmeEsl4Qa3WJ3+udkTrFsXOo agent
debug2: we sent a publickey packet, wait for reply
debug1: Server accepts key: /home/daniel/.ssh/id_rsa RSA SHA256:A6pKrXedkwxgG0XwUHtzmeEsl4Qa3WJ3+udkTrFsXOo agent
Authenticated to login1.yggdrasil.hpc.unige.ch ([129.194.64.11]:22) using "publickey".
debug1: channel 0: new [client-session]
debug2: channel 0: send open
debug1: Requesting no-more-sessions@openssh.com
debug1: Entering interactive session.
debug1: pledge: filesystem
debug1: client_input_global_request: rtype hostkeys-00@openssh.com want_reply 0
debug1: client_input_hostkeys: searching /home/daniel/.ssh/known_hosts for login1.yggdrasil.hpc.unige.ch / (none)
debug1: client_input_hostkeys: searching /home/daniel/.ssh/known_hosts2 for login1.yggdrasil.hpc.unige.ch / (none)
debug1: client_input_hostkeys: hostkeys file /home/daniel/.ssh/known_hosts2 does not exist
debug1: client_input_hostkeys: no new or deprecated keys from server
debug1: Remote: /usr/bin/sss_ssh_authorizedkeys:1: key options: agent-forwarding port-forwarding pty user-rc x11-forwarding
debug1: Remote: /usr/bin/sss_ssh_authorizedkeys:1: key options: agent-forwarding port-forwarding pty user-rc x11-forwarding
debug2: channel_input_open_confirmation: channel 0: callback start
debug2: fd 3 setting TCP_NODELAY
debug2: client_session2_setup: id 0
debug2: channel 0: request pty-req confirm 1
debug1: Sending environment.
debug1: channel 0: setting env LC_ADDRESS = "es_CO.UTF-8"
debug2: channel 0: request env confirm 0
debug1: channel 0: setting env LC_NAME = "es_CO.UTF-8"
debug2: channel 0: request env confirm 0
debug1: channel 0: setting env LC_MONETARY = "es_CO.UTF-8"
debug2: channel 0: request env confirm 0
debug1: channel 0: setting env LC_PAPER = "es_CO.UTF-8"
debug2: channel 0: request env confirm 0
debug1: channel 0: setting env LANG = "en_US.UTF-8"
debug2: channel 0: request env confirm 0
debug1: channel 0: setting env LC_IDENTIFICATION = "es_CO.UTF-8"
debug2: channel 0: request env confirm 0
debug1: channel 0: setting env LC_TELEPHONE = "es_CO.UTF-8"
debug2: channel 0: request env confirm 0
debug1: channel 0: setting env LC_MEASUREMENT = "es_CO.UTF-8"
debug2: channel 0: request env confirm 0
debug1: channel 0: setting env LC_TIME = "es_CO.UTF-8"
debug2: channel 0: request env confirm 0
debug1: channel 0: setting env LC_NUMERIC = "es_CO.UTF-8"
debug2: channel 0: request env confirm 0
debug2: channel 0: request shell confirm 1
debug2: channel_input_open_confirmation: channel 0: callback done
debug2: channel 0: open confirm rwindow 0 rmax 32768
debug2: channel_input_status_confirm: type 99 id 0
debug2: PTY allocation request accepted on channel 0
debug2: channel 0: rcvd adjust 2097152
debug2: channel_input_status_confirm: type 99 id 0
debug2: shell request accepted on channel 0
```


## Post 2 by @William.Ceva (2024-08-22T14:34:39.040Z)

Same, apologies for the duplicate thread: Connection to yggdrasil is never established[Connection to yggdrasil is never established](https://hpc-community.unige.ch/t/connection-to-yggdrasil-is-never-established/3603)


## Post 3 by @Adrien.Albert (2024-08-22T15:19:32.103Z)

Hi @daniel.forerosanchez[@daniel.forerosanchez](https://hpc-community.unige.ch/u/daniel.forerosanchez)
It seems that a user’s gvfs mount point blocked the login1 file system listing.
The login node has been restarted, so you should be able to log in again.


## Post 4 by @Adrien.Albert (2024-08-22T15:23:22.681Z)

Edit: The home filesystem is failed after reboot => investigating…


## Post 5 by @Adrien.Albert (2024-08-22T15:35:51.029Z)

Problem solved !!
Let me know if anything’s wrong!
