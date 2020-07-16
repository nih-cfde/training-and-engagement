---
layout: page
title: Accessing the AWS instance
---

Log in to your AWS instance on a Mac terminal
==============================================

OK, so you've created a running computer on the cloud. How do you get to it?

## Getting to the AWS instance

To access the cloud computer, you need the network name of our new computer. This can be found at the bottom of the instance log page shown here:

![](images/publicDNS.png)

Copy this name, and connect to the cloud computer with ssh under the username ‘ubuntu’, as follows:

* Find the private key file; it’s the `.pem` file you downloaded when starting up the EC2 instance. Remember, you named it `amazon.pem` and saved in on the Desktop. Open up the terminal, then type:

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

* If everything works ok, you should see something like this:


![](images/AWS_Connected.png)
!!! Note
    My terminal window is yellow, but yours may not be!
* You have now successfully logged in as user ‘ubuntu’ to the machine ‘ec2-18-216-20-166.us-east-2.compute.amazonaws.com’ using the authentication key located in ‘amazon.pem’ on the Desktop.
