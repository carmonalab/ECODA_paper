# Mandatory SSH Key Authentication for Baobab HPC Access

- Source: https://hpc-community.unige.ch/t/4318

- Created: 2026-06-19T06:16:55.914Z

- Posts: 16

- Category: 6

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yann.Sagon (2026-06-19T06:16:56.050Z)

Dear users,
This post is the same as the email we sent the 15th of June to every users.
For security reasons, we will require all users to authenticate using SSH keys instead of passwords when connecting to our Baobab HPC systems, this is valid for the three clusters Baobab, Bamboo, Yggdrasil.
:date: Effective date: June 22, 2026
You can already start using SSH keys today.
---
## :white_check_mark: Why this change?
Using an SSH key provides an additional layer of security.
In most phishing incidents, compromised users unintentionally reveal their password, but not their SSH private key. This significantly reduces the risk of unauthorized access.
---
## :key: What you need to do
- Start using SSH key authentication as soon as possible
- Make sure your public key is registered:
  - UNIGE users: https://my-account.unige.ch/[https://my-account.unige.ch/](https://my-account.unige.ch/)
  - Outsider users: https://applicant.unige.ch/main/outsider-info/update[https://applicant.unige.ch/main/outsider-info/update](https://applicant.unige.ch/main/outsider-info/update)
---
## :warning: Security requirements
- :locked: Protect your private key carefully
  - It must be secured with a strong passphrase
- :prohibited: Do NOT share your SSH key
  - SSH keys are strictly personal
  - If you need shared access, please contact us for alternative solutions
---
## :light_bulb: Recommended: Use an SSH agent
We strongly recommend using an SSH agent:
- It allows you to unlock your private key once per session
- It works seamlessly with many tools such as:
  - `ssh`
  - FileZilla
  - VS Code
  - and others
---
## :desktop_computer: Tools by platform
### Windows
Download: Download PuTTY: latest release (0.84)[Download PuTTY: latest release (0.84)](https://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html)
- PuTTY – SSH client
- Pageant – SSH agent
- PuTTYgen – SSH key generator
### macOS / Linux
Use native tools:
- `ssh`
- `ssh-agent`
- `ssh-keygen`
---
## :counterclockwise_arrows_button: Alternative access (if needed)
Users who cannot or do not wish to use SSH keys can still access the platform via:
:right_arrow: Open OnDemand
(using switch edu-ID + 2FA)
---
Thank you for your cooperation in improving the security of our HPC infrastructure.
Best regards,
HPC Team


## Post 2 by @Yann.Sagon (2026-06-19T06:19:49.340Z)

You can check which ssh key is currently enrolled with your username, directly from one of the clusters with the following command:
```
(bamboo)-[root@login1 ~]$ /usr/bin/sss_ssh_authorizedkeys sagon
ssh-rsa AAAAB3NzaC1yc2EAAAADAQABAAABAQCjTrDWG5vD2eOYuiJOnZ6n4wss3pYvbtJCP0S04PPHl+8oOvBSl9/5jV8Um8IsAtMniMAQ5mgxfIJSlsGY7J5G6IqH1brWDvPj7pYXsxGih3kidOpXI1sPzock9ocwlk5sgcVJkuSKsn/6woacHQlVajGNE7QG3/ZbuJ+EuL9gO0WUmEhu72bxVdA4/AaV0o1nA8F2xfy09t63OsqeUV1he25EqGkeI1v9J6dMg1ckXpHM8sG6KryZCrlysDu2DJ9zaQar2p3KTKcFnpn8oQr98D2pYA67ywrEZOWPntlYwkvrGztVs1JdwD4nbf5P080BBscGOK4ia+BClk715aLF rsa-key-20240228 Yann.Sagon@unige.ch 2024 registration
```
Best regards
Yann


## Post 3 by @Davide.Piras (2026-06-24T12:52:09.052Z)

Hello, I think I did everything correctly but I’m unable to login. Is there maybe some time one has to wait between copying the public key on the platform and logging in?
Thank you!


## Post 4 by @Davide.Piras (2026-06-24T12:59:11.261Z)

It worked now without me doing anything, so I assume it was just a matter of time.


## Post 5 by @Erica.Lastufka (2026-06-24T16:57:17.949Z)

Hello, I’m unable to register my public key to my UNIGE account. Whenever I copy the contents of ~/.ssh/id_rsa.pub to the box, it says wrong value or format. What exactly is it expecting?


## Post 6 by @Yann.Sagon (2026-06-25T06:37:50.808Z)

Important note:
there is a synchronization delay of about 5 minutes once you add/modify your ssh key to your account. Please be patient before trying to connect or your account risk to be banned due to too many authentification failure, but don’t worry the ban is temporary only and your account is unblock automatically. hpc:faq [eResearch Doc][hpc:faq [eResearch Doc]](https://doc.eresearch.unige.ch/hpc/faq#when_i_tried_to_connect_to_the)


## Post 7 by @Yann.Sagon (2026-06-25T07:54:46.154Z)

@Erica.Lastufka[@Erica.Lastufka](https://hpc-community.unige.ch/u/erica.lastufka) what is the format of your ssh pub key?


## Post 8 by @Erica.Lastufka (2026-06-25T08:19:41.189Z)

I changed the format to ed25519 and it was accepted by my-account.unige


## Post 9 by @Gael.Rossignol (2026-06-26T12:29:23.946Z)

Update: We have identified that the SSSD system utilizes a large cache. To optimize performance and accelerate key deployment, we have adjusted the cache’s Time To Live (TTL) to update every 5 minutes. This change will significantly improve deployment efficiency.


## Post 10 by @Cecile.Richetta (2026-07-02T11:52:58.686Z)

Dear HPC admin,
I can’t connect via SSH, and I am bamboozled (to say the least) as to why that would be the case.
- With PuTTYgen I created a private key, with a passphrase, saved locally on my .ssh folder
- I copy/pasted the ssh-rsa AAAA […] public key into my UNIGE account more than 50 minutes ago. No line breaks or copy/paste mistake. In my UNIGE account, it shows the key + “rsa-key-20260702 Cecile.Richetta@unige.ch 2026 registration” at the end
- In PuTTY, trying to connect to bamboo, I have the login node “login1.bamboo.hpc.unige.ch” on Port 22. In Data, my auto login username is “cecile.richetta@unige.ch” (although i’ve tried also with my short username, richetta). In Connection > SSH > Auth > Credentials > Private key file is pointing correctly to the .ppk file saved in my .ssh folder locally
- Trying to log in “Server refused our key” / "No supported authentication methods available (server sent:publickey,hostbased).
I truly do not see where I went wrong, and if you could help me in any way, I would appreciate it.
*EDIT: By logging through On Demand, and typing /usr/bin/sss_ssh_authorizedkeys richetta I see that an old SSH key is registered that I can no longer use. Any chance you could force a synch on your side?
Thanks, Cecile.


## Post 11 by @Adrien.Albert (2026-07-02T13:20:42.886Z)

Hello @Cecile.Richetta[@Cecile.Richetta](https://hpc-community.unige.ch/u/cecile.richetta)
Could you please try again and let me know it’s working ? Sometimes the synchronization can take a little longer than expected.


## Post 12 by @Cecile.Richetta (2026-07-02T13:25:54.709Z)

Hi Adrien, thanks for taking the time. Unfortunately it still hasn’t synched (I see the correct SSH key on my-account.unige.ch but still the oudated one with /usr/bin/sss_ssh_authorizedkeys richetta by logging through On Demand). It’s been a little over 2 hours now since I uploaded the key. Do you think it’s a CAD or HPC matter?


## Post 13 by @Adrien.Albert (2026-07-03T07:32:40.588Z)

Hi @Cecile.Richetta[@Cecile.Richetta](https://hpc-community.unige.ch/u/cecile.richetta),
We have encountered a case where a user attempted to register the same public key twice in My Account, resulting in a silent blocking issue in LDAP. As far as we can tell, this situation can only be resolved by an LDAP administrator.
Does this case sound familiar to you?
I have also sent you an email with additional details for further support.
Best


## Post 14 by @Cecile.Richetta (2026-07-03T07:55:44.389Z)

Hi Adrien, thanks again. Yes it’s possible I did this, because I struggled to set up the SSH access so I might very well have copy/pasted the same one twice. I sent a request through the CAD, hopefully that doesn’t end on your side but on the side of the LDAP admin – there are no specific options for SSH issues on the CAD.
Best


## Post 15 by @Jaspreet.Saini (2026-07-06T08:53:23.284Z)

Hi everyone,
I seem to be experiencing the same issue.
I uploaded my ED25519 public key to the UNIGE portal and waited for propagation, but I still cannot log in to Baobab.
Username: `sainij`
Host: `login1.baobab.hpc.unige.ch`
My SSH client is correctly offering the key:
```
Offering public key: /Users/sainij/.ssh/id_ed25519 ED25519 SHA256:0oXBJe2lT8xhy2xoUp5FNAV/Sf6QFatkOaJed0fSZoI explicit
```
but the server responds with:
```
Permission denied (publickey,hostbased).
```
I have also tried:
```
ssh -i ~/.ssh/id_ed25519 -o IdentitiesOnly=yes sainij@login1.baobab.hpc.unige.ch
```
with the same result.
Is anyone else still experiencing this issue, or could HPC support please check whether my public key has been correctly propagated/mapped to my account?
Thanks in advance,
Jaspreet


## Post 16 by @Yann.Sagon (2026-07-08T14:09:58.583Z)

Dear @Jaspreet.Saini[@Jaspreet.Saini](https://hpc-community.unige.ch/u/jaspreet.saini)
Please check this post:
[tutorial] locked_with_key SSH login troubleshooting (key mismatch)[[tutorial]  SSH login troubleshooting (key mismatch)](https://hpc-community.unige.ch/t/tutorial-ssh-login-troubleshooting-key-mismatch/4339) HPC Support[HPC Support](https://hpc-community.unige.ch/c/hpc-support/13)
> locked_with_key SSH login troubleshooting (key mismatch)
If you cannot connect to the Baobab cluster via SSH, it is often due to a mismatch between: 

the private key on your machine
the public key stored in LDAP

You can easily verify this by comparing their fingerprints. 

white_check_mark Step 1 — Check your local key
On your computer, run: 
ssh-keygen -lf ~/.ssh/id_rsa

(Replace id_rsa with id_ed25519 or your key file if different.) 
This will display a fingerprint like: 
2048 SHA256:XXX…
