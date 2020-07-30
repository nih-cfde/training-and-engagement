---
layout: page
title: Access the AWS instance
---

Log in to your AWS instance on a Mac terminal
==============================================

OK, so you've created a running computer on the cloud. How do you get to it? We will use the Mac terminal window to access the cloud computer.

Start this tutorial by opening up a terminal window by searching (type cmd+space_bar) for "terminal" on your Mac.

## Getting to the AWS instance

To access the cloud computer, you need the network name of our new computer. This can be found at the bottom of the instance log page shown here:

![](images/publicDNS.png)

Copy this name, connect to the cloud computer with ssh under the username ‘ubuntu’, as follows:

* Find the private key file; it’s the `.pem` file you downloaded when starting up the EC2 instance. Remember, you named it `amazon.pem` and saved in on the Desktop. Select the terminal window, then type:

```
chmod og-rwx ~/Desktop/amazon.pem
```

* This sets the permissions on the private key file to “closed to all evildoers”.

* Finally, log in to the cloud computer:

```
ssh -i ~/Desktop/amazon.pem ubuntu@ec2-???-???-???-???.compute-1.amazonaws.com
```

!!! Important
    Replace the stuff after the ‘@’ sign with the name of the host; see the red circle on your own instance:
    ![](images/publicDNS.png)


!!! Error
    If you see this message:

    The authenticity of host 'ecc2-???-???-???-???.compute-1.amazonaws.com (3.129.57.169)' can't be established.
    ECDSA key fingerprint is XXX.
    Are you sure you want to continue connecting (yes/no/[fingerprint])? **yes**

    Type "yes" and press enter.



* If everything works ok, you should see something like this:

![](images/AWS_Connected.png)

!!! Note
    My terminal window is yellow, but yours may not be!

* You have now successfully logged in as user ‘ubuntu’ to the machine ‘ec2-18-216-20-166.us-east-2.compute.amazonaws.com’ using the authentication key located in ‘amazon.pem’ on the Desktop.
