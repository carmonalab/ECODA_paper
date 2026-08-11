# Connection baobab "failed"

- Source: https://hpc-community.unige.ch/t/4066

- Created: 2025-08-26T12:33:02.334Z

- Tags: baobab

- Posts: 3

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Thibaut.Chataing (2025-08-26T12:33:02.404Z)

Hi,
I can’t connect to baobab by ssh. When I try to log the connection is immediately closed.
I connect to baobab in ssh though a key. So, I use the following command :
```
$ ssh chataint@login1.baobab.hpc.unige.ch
```
It “works”. After some times I receive this message (with the closed line)
> Last login: Tue Aug 26 11:18:31 2025 from app6.baobab
---
> |  _ \            | |         | |
> | |) | __ _  ___ | |__   __ | |_
> |  _ < /  `|/ _ \| '_ \ / _` | ’
> | |) | (| | () | |) | (| | |) |
> |/ _,|__/|./ _,|_./
> _             _      __
> | |           ()    / |
> | | __   __ _ _ _ __ | |
> | |/ _ \ / ` | | ’ | |
> | | () | (| | | | | | |
> ||___/ _, ||| |||
> / |
> |_/
> Documentation: hpc:start [eResearch Doc][hpc:start [eResearch Doc]](https://doc.eresearch.unige.ch/hpc/start)
> Forum: https://hpc-community.unige.ch/[https://hpc-community.unige.ch/](https://hpc-community.unige.ch/)
> OpenOndemand: https://openondemand.baobab.hpc.unige.ch/[https://openondemand.baobab.hpc.unige.ch/](https://openondemand.baobab.hpc.unige.ch/)
> support: hpc:start [eResearch Doc][hpc:start [eResearch Doc]](https://doc.eresearch.unige.ch/hpc/start#support%5C_-%5C_get_help)
> Connection to login1.baobab.hpc.unige.ch closed.
On the other hand, I manage to connect though OnDemand.
My quota is ok :
```
  user/group                 ||           size          ||    chunk files
```
> storage     |   name        |  id  ||    used    |    hard    ||  used   |  hard
> ----------------------------|------||------------|------------||---------|---------
> home        |       chataint|418159||  518.97 GiB| 1024.00 GiB||   365828|unlimited
> scratch     |       chataint|418159||      0 Byte|   unlimited||        0| 10000000
Here are the detail of the ssh connection with the debug argument -vvv :
```
PS C:\Users\chataint\Documents> ssh chataint@login1.baobab.hpc.unige.ch -vvv
OpenSSH_for_Windows_9.5p1, LibreSSL 3.8.2
debug1: Reading configuration data C:\\Users\\chataint/.ssh/config
debug1: C:\\Users\\chataint/.ssh/config line 1: Applying options for login1.baobab.hpc.unige.ch
debug3: Failed to open file:C:/ProgramData/ssh/ssh_config error:2
debug3: expanded UserKnownHostsFile '~/.ssh/known_hosts' -> 'C:\\Users\\chataint/.ssh/known_hosts'
debug3: expanded UserKnownHostsFile '~/.ssh/known_hosts2' -> 'C:\\Users\\chataint/.ssh/known_hosts2'
debug2: resolving "login1.baobab.hpc.unige.ch" port 22
debug3: resolve_host: lookup login1.baobab.hpc.unige.ch:22
debug3: ssh_connect_direct: entering
debug1: Connecting to login1.baobab.hpc.unige.ch [129.194.9.190] port 22.
debug1: Connection established.
debug1: identity file C:\\Users\\chataint/.ssh/id_rsa type 0
debug3: Failed to open file:C:/Users/chataint/.ssh/id_rsa-cert error:2
debug3: Failed to open file:C:/Users/chataint/.ssh/id_rsa-cert.pub error:2
debug3: failed to open file:C:/Users/chataint/.ssh/id_rsa-cert error:2
debug1: identity file C:\\Users\\chataint/.ssh/id_rsa-cert type -1
debug1: identity file C:\\Users\\chataint/.ssh/id_ecdsa type 2
debug3: Failed to open file:C:/Users/chataint/.ssh/id_ecdsa-cert error:2
debug3: Failed to open file:C:/Users/chataint/.ssh/id_ecdsa-cert.pub error:2
debug3: failed to open file:C:/Users/chataint/.ssh/id_ecdsa-cert error:2
debug1: identity file C:\\Users\\chataint/.ssh/id_ecdsa-cert type -1
debug3: Failed to open file:C:/Users/chataint/.ssh/id_ecdsa_sk error:2
debug3: Failed to open file:C:/Users/chataint/.ssh/id_ecdsa_sk.pub error:2
debug3: failed to open file:C:/Users/chataint/.ssh/id_ecdsa_sk error:2
debug1: identity file C:\\Users\\chataint/.ssh/id_ecdsa_sk type -1
debug3: Failed to open file:C:/Users/chataint/.ssh/id_ecdsa_sk-cert error:2
debug3: Failed to open file:C:/Users/chataint/.ssh/id_ecdsa_sk-cert.pub error:2
debug3: failed to open file:C:/Users/chataint/.ssh/id_ecdsa_sk-cert error:2
debug1: identity file C:\\Users\\chataint/.ssh/id_ecdsa_sk-cert type -1
debug1: identity file C:\\Users\\chataint/.ssh/id_ed25519 type 3
debug3: Failed to open file:C:/Users/chataint/.ssh/id_ed25519-cert error:2
debug3: Failed to open file:C:/Users/chataint/.ssh/id_ed25519-cert.pub error:2
debug3: failed to open file:C:/Users/chataint/.ssh/id_ed25519-cert error:2
debug1: identity file C:\\Users\\chataint/.ssh/id_ed25519-cert type -1
debug3: Failed to open file:C:/Users/chataint/.ssh/id_ed25519_sk error:2
debug3: Failed to open file:C:/Users/chataint/.ssh/id_ed25519_sk.pub error:2
debug3: failed to open file:C:/Users/chataint/.ssh/id_ed25519_sk error:2
debug1: identity file C:\\Users\\chataint/.ssh/id_ed25519_sk type -1
debug3: Failed to open file:C:/Users/chataint/.ssh/id_ed25519_sk-cert error:2
debug3: Failed to open file:C:/Users/chataint/.ssh/id_ed25519_sk-cert.pub error:2
debug3: failed to open file:C:/Users/chataint/.ssh/id_ed25519_sk-cert error:2
debug1: identity file C:\\Users\\chataint/.ssh/id_ed25519_sk-cert type -1
debug3: Failed to open file:C:/Users/chataint/.ssh/id_xmss error:2
debug3: Failed to open file:C:/Users/chataint/.ssh/id_xmss.pub error:2
debug3: failed to open file:C:/Users/chataint/.ssh/id_xmss error:2
debug1: identity file C:\\Users\\chataint/.ssh/id_xmss type -1
debug3: Failed to open file:C:/Users/chataint/.ssh/id_xmss-cert error:2
debug3: Failed to open file:C:/Users/chataint/.ssh/id_xmss-cert.pub error:2
debug3: failed to open file:C:/Users/chataint/.ssh/id_xmss-cert error:2
debug1: identity file C:\\Users\\chataint/.ssh/id_xmss-cert type -1
debug3: Failed to open file:C:/Users/chataint/.ssh/id_dsa error:2
debug3: Failed to open file:C:/Users/chataint/.ssh/id_dsa.pub error:2
debug3: failed to open file:C:/Users/chataint/.ssh/id_dsa error:2
debug1: identity file C:\\Users\\chataint/.ssh/id_dsa type -1
debug3: Failed to open file:C:/Users/chataint/.ssh/id_dsa-cert error:2
debug3: Failed to open file:C:/Users/chataint/.ssh/id_dsa-cert.pub error:2
debug3: failed to open file:C:/Users/chataint/.ssh/id_dsa-cert error:2
debug1: identity file C:\\Users\\chataint/.ssh/id_dsa-cert type -1
debug1: Local version string SSH-2.0-OpenSSH_for_Windows_9.5
debug1: Remote protocol version 2.0, remote software version OpenSSH_8.7
debug1: compat_banner: match: OpenSSH_8.7 pat OpenSSH* compat 0x04000000
debug2: fd 3 setting O_NONBLOCK
debug1: Authenticating to login1.baobab.hpc.unige.ch:22 as 'chataint'
debug3: record_hostkey: found key type RSA in file C:\\Users\\chataint/.ssh/known_hosts:7
debug3: load_hostkeys_file: loaded 1 keys from login1.baobab.hpc.unige.ch
debug3: Failed to open file:C:/Users/chataint/.ssh/known_hosts2 error:2
debug1: load_hostkeys: fopen C:\\Users\\chataint/.ssh/known_hosts2: No such file or directory
debug3: Failed to open file:C:/ProgramData/ssh/ssh_known_hosts error:2
debug1: load_hostkeys: fopen __PROGRAMDATA__\\ssh/ssh_known_hosts: No such file or directory
debug3: Failed to open file:C:/ProgramData/ssh/ssh_known_hosts2 error:2
debug1: load_hostkeys: fopen __PROGRAMDATA__\\ssh/ssh_known_hosts2: No such file or directory
debug3: order_hostkeyalgs: prefer hostkeyalgs: rsa-sha2-512-cert-v01@openssh.com,rsa-sha2-256-cert-v01@openssh.com,rsa-sha2-512,rsa-sha2-256
debug3: send packet: type 20
debug1: SSH2_MSG_KEXINIT sent
debug3: receive packet: type 20
debug1: SSH2_MSG_KEXINIT received
debug2: local client KEXINIT proposal
debug2: KEX algorithms: curve25519-sha256,curve25519-sha256@libssh.org,ecdh-sha2-nistp256,ecdh-sha2-nistp384,ecdh-sha2-nistp521,diffie-hellman-group-exchange-sha256,diffie-hellman-group16-sha512,diffie-hellman-group18-sha512,diffie-hellman-group14-sha256,ext-info-c,kex-strict-c-v00@openssh.com
debug2: host key algorithms: rsa-sha2-512-cert-v01@openssh.com,rsa-sha2-256-cert-v01@openssh.com,rsa-sha2-512,rsa-sha2-256,ssh-ed25519-cert-v01@openssh.com,ecdsa-sha2-nistp256-cert-v01@openssh.com,ecdsa-sha2-nistp384-cert-v01@openssh.com,ecdsa-sha2-nistp521-cert-v01@openssh.com,sk-ssh-ed25519-cert-v01@openssh.com,sk-ecdsa-sha2-nistp256-cert-v01@openssh.com,ssh-ed25519,ecdsa-sha2-nistp256,ecdsa-sha2-nistp384,ecdsa-sha2-nistp521,sk-ssh-ed25519@openssh.com,sk-ecdsa-sha2-nistp256@openssh.com
debug2: ciphers ctos: chacha20-poly1305@openssh.com,aes128-ctr,aes192-ctr,aes256-ctr,aes128-gcm@openssh.com,aes256-gcm@openssh.com
debug2: ciphers stoc: chacha20-poly1305@openssh.com,aes128-ctr,aes192-ctr,aes256-ctr,aes128-gcm@openssh.com,aes256-gcm@openssh.com
debug2: MACs ctos: umac-64-etm@openssh.com,umac-128-etm@openssh.com,hmac-sha2-256-etm@openssh.com,hmac-sha2-512-etm@openssh.com,umac-64@openssh.com,umac-128@openssh.com,hmac-sha2-256,hmac-sha2-512
debug2: MACs stoc: umac-64-etm@openssh.com,umac-128-etm@openssh.com,hmac-sha2-256-etm@openssh.com,hmac-sha2-512-etm@openssh.com,umac-64@openssh.com,umac-128@openssh.com,hmac-sha2-256,hmac-sha2-512
debug2: compression ctos: none,zlib@openssh.com,zlib
debug2: compression stoc: none,zlib@openssh.com,zlib
debug2: languages ctos:
debug2: languages stoc:
debug2: first_kex_follows 0
debug2: reserved 0
debug2: peer server KEXINIT proposal
debug2: KEX algorithms: curve25519-sha256,curve25519-sha256@libssh.org,ecdh-sha2-nistp256,ecdh-sha2-nistp384,ecdh-sha2-nistp521,diffie-hellman-group-exchange-sha256,diffie-hellman-group16-sha512,diffie-hellman-group18-sha512,diffie-hellman-group14-sha256,kex-strict-s-v00@openssh.com
debug2: host key algorithms: rsa-sha2-512,rsa-sha2-256,ssh-rsa
debug2: ciphers ctos: chacha20-poly1305@openssh.com,aes128-ctr,aes192-ctr,aes256-ctr,aes128-gcm@openssh.com,aes256-gcm@openssh.com
debug2: ciphers stoc: chacha20-poly1305@openssh.com,aes128-ctr,aes192-ctr,aes256-ctr,aes128-gcm@openssh.com,aes256-gcm@openssh.com
debug2: MACs ctos: umac-64-etm@openssh.com,umac-128-etm@openssh.com,hmac-sha2-256-etm@openssh.com,hmac-sha2-512-etm@openssh.com,hmac-sha1-etm@openssh.com,umac-64@openssh.com,umac-128@openssh.com,hmac-sha2-256,hmac-sha2-512,hmac-sha1
debug2: MACs stoc: umac-64-etm@openssh.com,umac-128-etm@openssh.com,hmac-sha2-256-etm@openssh.com,hmac-sha2-512-etm@openssh.com,hmac-sha1-etm@openssh.com,umac-64@openssh.com,umac-128@openssh.com,hmac-sha2-256,hmac-sha2-512,hmac-sha1
debug2: compression ctos: none,zlib@openssh.com
debug2: compression stoc: none,zlib@openssh.com
debug2: languages ctos:
debug2: languages stoc:
debug2: first_kex_follows 0
debug2: reserved 0
debug3: kex_choose_conf: will use strict KEX ordering
debug1: kex: algorithm: curve25519-sha256
debug1: kex: host key algorithm: rsa-sha2-512
debug1: kex: server->client cipher: chacha20-poly1305@openssh.com MAC: <implicit> compression: none
debug1: kex: client->server cipher: chacha20-poly1305@openssh.com MAC: <implicit> compression: none
debug3: send packet: type 30
debug1: expecting SSH2_MSG_KEX_ECDH_REPLY
debug3: receive packet: type 31
debug1: SSH2_MSG_KEX_ECDH_REPLY received
debug1: Server host key: ssh-rsa SHA256:tKqp4nljL+EGVKl8T0VF2nS36DkHVFMpLxQOPg/gKvg
debug3: record_hostkey: found key type RSA in file C:\\Users\\chataint/.ssh/known_hosts:7
debug3: load_hostkeys_file: loaded 1 keys from login1.baobab.hpc.unige.ch
debug3: Failed to open file:C:/Users/chataint/.ssh/known_hosts2 error:2
debug1: load_hostkeys: fopen C:\\Users\\chataint/.ssh/known_hosts2: No such file or directory
debug3: Failed to open file:C:/ProgramData/ssh/ssh_known_hosts error:2
debug1: load_hostkeys: fopen __PROGRAMDATA__\\ssh/ssh_known_hosts: No such file or directory
debug3: Failed to open file:C:/ProgramData/ssh/ssh_known_hosts2 error:2
debug1: load_hostkeys: fopen __PROGRAMDATA__\\ssh/ssh_known_hosts2: No such file or directory
debug1: Host 'login1.baobab.hpc.unige.ch' is known and matches the RSA host key.
debug1: Found key in C:\\Users\\chataint/.ssh/known_hosts:7
debug3: send packet: type 21
debug1: ssh_packet_send2_wrapped: resetting send seqnr 3
debug2: ssh_set_newkeys: mode 1
debug1: rekey out after 134217728 blocks
debug1: SSH2_MSG_NEWKEYS sent
debug1: expecting SSH2_MSG_NEWKEYS
debug3: receive packet: type 21
debug1: ssh_packet_read_poll2: resetting read seqnr 3
debug1: SSH2_MSG_NEWKEYS received
debug2: ssh_set_newkeys: mode 0
debug1: rekey in after 134217728 blocks
debug3: ssh_get_authentication_socket_path: path '\\\\.\\pipe\\openssh-ssh-agent'
debug3: unable to connect to pipe \\\\.\\pipe\\openssh-ssh-agent, error: 2
debug1: get_agent_identities: ssh_get_authentication_socket: No such file or directory
debug1: Will attempt key: C:\\Users\\chataint/.ssh/id_rsa RSA SHA256:WcKKmxoxx2NZ3N+z7J5SIrJ4qFbjCdNSuAkuIQrwykk
debug1: Will attempt key: C:\\Users\\chataint/.ssh/id_ecdsa ECDSA SHA256:ziABSoHMXiWyB2Anxb5c0ALAJ5j+mscFYWWZ3G4ajQU
debug1: Will attempt key: C:\\Users\\chataint/.ssh/id_ecdsa_sk
debug1: Will attempt key: C:\\Users\\chataint/.ssh/id_ed25519 ED25519 SHA256:sb2IelnPfxIXGXeUF80/Z7f+UOGijjoUOEmuIa2lMxU
debug1: Will attempt key: C:\\Users\\chataint/.ssh/id_ed25519_sk
debug1: Will attempt key: C:\\Users\\chataint/.ssh/id_xmss
debug1: Will attempt key: C:\\Users\\chataint/.ssh/id_dsa
debug2: pubkey_prepare: done
debug3: send packet: type 5
debug3: receive packet: type 7
debug1: SSH2_MSG_EXT_INFO received
debug1: kex_input_ext_info: server-sig-algs=<ssh-ed25519,sk-ssh-ed25519@openssh.com,ssh-rsa,rsa-sha2-256,rsa-sha2-512,ssh-dss,ecdsa-sha2-nistp256,ecdsa-sha2-nistp384,ecdsa-sha2-nistp521,sk-ecdsa-sha2-nistp256@openssh.com,webauthn-sk-ecdsa-sha2-nistp256@openssh.com>
debug3: receive packet: type 6
debug2: service_accept: ssh-userauth
debug1: SSH2_MSG_SERVICE_ACCEPT received
debug3: send packet: type 50
debug3: receive packet: type 51
debug1: Authentications that can continue: publickey,password,keyboard-interactive,hostbased
debug3: start over, passed a different list publickey,password,keyboard-interactive,hostbased
debug3: preferred publickey,keyboard-interactive,password
debug3: authmethod_lookup publickey
debug3: remaining preferred: keyboard-interactive,password
debug3: authmethod_is_enabled publickey
debug1: Next authentication method: publickey
debug1: Offering public key: C:\\Users\\chataint/.ssh/id_rsa RSA SHA256:WcKKmxoxx2NZ3N+z7J5SIrJ4qFbjCdNSuAkuIQrwykk
debug3: send packet: type 50
debug2: we sent a publickey packet, wait for reply
debug3: receive packet: type 51
debug1: Authentications that can continue: publickey,password,keyboard-interactive,hostbased
debug1: Offering public key: C:\\Users\\chataint/.ssh/id_ecdsa ECDSA SHA256:ziABSoHMXiWyB2Anxb5c0ALAJ5j+mscFYWWZ3G4ajQU
debug3: send packet: type 50
debug2: we sent a publickey packet, wait for reply
debug3: receive packet: type 51
debug1: Authentications that can continue: publickey,password,keyboard-interactive,hostbased
debug1: Trying private key: C:\\Users\\chataint/.ssh/id_ecdsa_sk
debug3: no such identity: C:\\Users\\chataint/.ssh/id_ecdsa_sk: No such file or directory
debug1: Offering public key: C:\\Users\\chataint/.ssh/id_ed25519 ED25519 SHA256:sb2IelnPfxIXGXeUF80/Z7f+UOGijjoUOEmuIa2lMxU
debug3: send packet: type 50
debug2: we sent a publickey packet, wait for reply
debug3: receive packet: type 60
debug1: Server accepts key: C:\\Users\\chataint/.ssh/id_ed25519 ED25519 SHA256:sb2IelnPfxIXGXeUF80/Z7f+UOGijjoUOEmuIa2lMxU
debug3: sign_and_send_pubkey: using publickey with ED25519 SHA256:sb2IelnPfxIXGXeUF80/Z7f+UOGijjoUOEmuIa2lMxU
debug3: sign_and_send_pubkey: signing using ssh-ed25519 SHA256:sb2IelnPfxIXGXeUF80/Z7f+UOGijjoUOEmuIa2lMxU
debug3: send packet: type 50
debug3: receive packet: type 52
Authenticated to login1.baobab.hpc.unige.ch ([129.194.9.190]:22) using "publickey".
debug1: channel 0: new session [client-session] (inactive timeout: 0)
debug3: ssh_session2_open: channel_new: 0
debug2: channel 0: send open
debug3: send packet: type 90
debug1: Requesting no-more-sessions@openssh.com
debug3: send packet: type 80
debug1: Entering interactive session.
debug1: pledge: filesystem
debug3: client_repledge: enter
debug1: ENABLE_VIRTUAL_TERMINAL_INPUT is supported. Reading the VTSequence from console
debug3: This windows OS supports conpty
debug1: ENABLE_VIRTUAL_TERMINAL_PROCESSING is supported. Console supports the ansi parsing
debug3: Successfully set console output code page from:65001 to 65001
debug3: Successfully set console input code page from:437 to 65001
debug3: receive packet: type 80
debug1: client_input_global_request: rtype hostkeys-00@openssh.com want_reply 0
debug3: client_input_hostkeys: received RSA key SHA256:tKqp4nljL+EGVKl8T0VF2nS36DkHVFMpLxQOPg/gKvg
debug1: client_input_hostkeys: searching C:\\Users\\chataint/.ssh/known_hosts for login1.baobab.hpc.unige.ch / (none)
debug3: hostkeys_foreach: reading file "C:\\Users\\chataint/.ssh/known_hosts"
debug3: hostkeys_find: found ssh-rsa key at C:\\Users\\chataint/.ssh/known_hosts:7
debug3: hostkeys_find: found ssh-rsa key under different name/addr at C:\\Users\\chataint/.ssh/known_hosts:8
debug3: hostkeys_find: found ssh-rsa key under different name/addr at C:\\Users\\chataint/.ssh/known_hosts:9
debug3: hostkeys_find: found ssh-rsa key under different name/addr at C:\\Users\\chataint/.ssh/known_hosts:10
debug3: hostkeys_find: found ssh-rsa key under different name/addr at C:\\Users\\chataint/.ssh/known_hosts:11
debug3: hostkeys_find: found ssh-rsa key under different name/addr at C:\\Users\\chataint/.ssh/known_hosts:12
debug3: hostkeys_find: found ssh-rsa key under different name/addr at C:\\Users\\chataint/.ssh/known_hosts:13
debug3: hostkeys_find: found ssh-rsa key under different name/addr at C:\\Users\\chataint/.ssh/known_hosts:14
debug3: hostkeys_find: found ssh-rsa key under different name/addr at C:\\Users\\chataint/.ssh/known_hosts:15
debug3: hostkeys_find: found ssh-rsa key under different name/addr at C:\\Users\\chataint/.ssh/known_hosts:16
debug3: hostkeys_find: found ssh-rsa key under different name/addr at C:\\Users\\chataint/.ssh/known_hosts:17
debug3: hostkeys_find: found ssh-rsa key under different name/addr at C:\\Users\\chataint/.ssh/known_hosts:18
debug3: hostkeys_find: found ssh-rsa key under different name/addr at C:\\Users\\chataint/.ssh/known_hosts:19
debug1: client_input_hostkeys: searching C:\\Users\\chataint/.ssh/known_hosts2 for login1.baobab.hpc.unige.ch / (none)
debug3: Failed to open file:C:/Users/chataint/.ssh/known_hosts2 error:2
debug1: client_input_hostkeys: hostkeys file C:\\Users\\chataint/.ssh/known_hosts2 does not exist
debug3: client_input_hostkeys: 1 server keys: 0 new, 1 retained, 0 incomplete match. 0 to remove
debug1: client_input_hostkeys: no new or deprecated keys from server
debug3: client_repledge: enter
debug3: receive packet: type 4
debug1: Remote: /usr/bin/sss_ssh_authorizedkeys:1: key options: agent-forwarding port-forwarding pty user-rc x11-forwarding
debug3: receive packet: type 4
debug1: Remote: /usr/bin/sss_ssh_authorizedkeys:1: key options: agent-forwarding port-forwarding pty user-rc x11-forwarding
debug3: receive packet: type 91
debug2: channel_input_open_confirmation: channel 0: callback start
debug2: fd 3 setting TCP_NODELAY
debug2: client_session2_setup: id 0
debug2: channel 0: request pty-req confirm 1
debug3: send packet: type 98
debug2: channel 0: request shell confirm 1
debug3: send packet: type 98
debug3: client_repledge: enter
debug1: pledge: fork
debug2: channel_input_open_confirmation: channel 0: callback done
debug2: channel 0: open confirm rwindow 0 rmax 32768
debug3: receive packet: type 99
debug2: channel_input_status_confirm: type 99 id 0
debug2: PTY allocation request accepted on channel 0
debug2: channel 0: rcvd adjust 2097152
debug3: receive packet: type 99
debug2: channel_input_status_confirm: type 99 id 0
debug2: shell request accepted on channel 0
Last login: Tue Aug 26 14:15:45 2025 from 10.20.211.149
 ____              _           _
|  _ \            | |         | |
| |_) | __ _  ___ | |__   __ _| |__
|  _ < / _` |/ _ \| '_ \ / _` | '_ \
| |_) | (_| | (_) | |_) | (_| | |_) |
|____/ \__,_|\___/|_.__/ \__,_|_.__/
             _             _      __
            | |           (_)    /_ |
            | | ___   __ _ _ _ __ | |
            | |/ _ \ / _` | | '_ \| |
            | | (_) | (_| | | | | | |
            |_|\___/ \__, |_|_| |_|_|
                      __/ |
                      |___/

 Documentation: https://doc.eresearch.unige.ch/hpc/start
 Forum: https://hpc-community.unige.ch/
 OpenOndemand: https://openondemand.baobab.hpc.unige.ch/
 support: https://doc.eresearch.unige.ch/hpc/start#support_-_get_help

debug3: receive packet: type 96
debug2: channel 0: rcvd eof
debug2: channel 0: output open -> drain
debug2: channel 0: obuf empty
debug2: chan_shutdown_write: channel 0: (i0 o1 sock -1 wfd 5 efd 6 [write])
debug2: channel 0: output drain -> closed
debug3: receive packet: type 98
debug1: client_input_channel_req: channel 0 rtype exit-status reply 0
debug3: receive packet: type 98
debug1: client_input_channel_req: channel 0 rtype eow@openssh.com reply 0
debug2: channel 0: rcvd eow
debug2: chan_shutdown_read: channel 0: (i0 o3 sock -1 wfd 4 efd 6 [write])
debug2: channel 0: input open -> closed
debug3: receive packet: type 97
debug2: channel 0: rcvd close
debug3: channel 0: will not send data after close
debug2: channel 0: almost dead
debug2: channel 0: gc: notify user
debug3: Successfully set console output code page from 65001 to 65001
debug3: Successfully set console input code page from 65001 to 437
debug2: channel 0: gc: user detached
debug2: channel 0: send close
debug3: send packet: type 97
debug2: channel 0: is dead
debug2: channel 0: garbage collecting
debug1: channel 0: free: client-session, nchannels 1
debug3: channel 0: status: The following connections are open:
  #0 client-session (t4 [session] r0 i3/0 o3/0 e[write]/0 fd -1/-1/6 sock -1 cc -1 io 0x00/0x00)

debug3: send packet: type 1
Connection to login1.baobab.hpc.unige.ch closed.
Transferred: sent 2832, received 3992 bytes, in 121.1 seconds
Bytes per second: sent 23.4, received 33.0
debug1: Exit status 254
```
Can you help me, pls ?


## Post 2 by @Adrien.Albert (2025-08-27T08:32:02.774Z)

Hi @Thibaut.Chataing[@Thibaut.Chataing](https://hpc-community.unige.ch/u/thibaut.chataing)
There is many errors, what do you use to connect on clusters ?
I have cleaned all remaining processes on login1.baobab and moved your `bashrc` `bash_profile` could you try again ?


## Post 3 by @Thibaut.Chataing (2025-08-27T09:32:50.654Z)

Hi,
It seems the cleaning work. It is now working. Thank you
To connect, i log with the terminal through the command ssh chataint@login1.baobab.hpc.unige.ch in windows.
