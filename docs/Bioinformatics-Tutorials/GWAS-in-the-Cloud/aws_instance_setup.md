---
layout: page
title: Setting up an AWS Instance
---

Setting up an AWS Instance
==========================


Amazon offers a cloud computing platform called Amazon Web Services (AWS). AWS is not free, however, you receive the benefits of the Free Tier automatically for 12 months after you sign up for an AWS account. If you are no longer eligible for the Free Tier, you're charged at the [standard billing rate](https://docs.aws.amazon.com/awsaccountbilling/latest/aboutv2/free-tier-eligibility.html) for your AWS usage.

## Free Tier AWS Ubuntu Instance

Let's set up a Free Tier AWS Ubuntu instance:

* Go to http://aws.amazon.com/ and click on the create an AWS account button on the top right. Proceed to create your account. If you have an existing AWS account, click the `Sign in to an existing AWS account` below the `continue` option on the sign-up page and log in to your account as a root user.

* Next, click on the 'Launch a virtual machine' option as shown in the image:

![](images/Launch.png)

* On the left navbar, check the Free tier only option, then select an Ubuntu machine

![](images/Ubuntu.png)

* Make sure the free tier version is selected. Then click review and launch, then launch. You should see a pop-up window like this:

![](images/KeyPair.png)

* Choose the 'Create a new key pair' option from the drop down menu. Under key pair name, type 'amazon'. Save it on your Desktop. Check the acknowledgement box and click "Launch Instance". Next time you launch an instance, you can reuse the key pair you just generated.

* After you click 'Launch Instance', you should see this:

![](images/launching.png)

* Click on this first hyperlink: i-038c58bfbe9612c57

* You have now successfully launched your AWS instance. DO NOT CLOSE THIS PAGE YET.
