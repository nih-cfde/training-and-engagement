# Introduction to Amazon Web Services

Summary:

* Go to [http://aws.amazon.com/](http://aws.amazon.com/) and log in, then select "EC2" (upper left)
* Select "Launch instance"
* Select "Ubuntu Pro 20.04 LTS" from the list
* Select "t2.micro" from the list 
* Click "Review and launch"
* Select "Launch"
* If your first time through, create a key pair; otherwise select existing
* Click "launch instance"

# Details

## What is Amazon Web Services?

Amazon Web Services is a subsidiary of Amazon that provides on-demand cloud computing platforms and APIs to individuals, companies and governments, on a metered pay-as-you-go basis. Subscribers can pay for a single virtual AWS computer, a dedicated physical computer or clusters of either.

The aim of this tutorial is to expose researchers to cloud computing resources and their benefits.

Amazon Web Services will require you to provide contact information (e-mail, etc.) and billing information (as you will have to pay for some of their services).
However, the _AWS Free Tier program_ allows you to use more than 60 products to start building with limited amount of data storage, limited amount of time **or** are always free.

To start we will learn how to set-up Ubuntu 20.04 Pro LTS open source software operating system. Ubuntu 20.04 Pro LTS is one of the programs offered in Amazon Free Tier as well as one of the most popular open source operating systems. 

## Let's Start

<iframe id="kaltura_player" src="https://cdnapisec.kaltura.com/p/1770401/sp/177040100/embedIframeJs/uiconf_id/29032722/partner_id/1770401?iframeembed=true&playerId=kaltura_player&entry_id=0_1a78cmhv&flashvars[mediaProtocol]=rtmp&amp;flashvars[streamerType]=rtmp&amp;flashvars[streamerUrl]=rtmp://www.kaltura.com:1935&amp;flashvars[rtmpFlavors]=1&amp;flashvars[localizationCode]=en&amp;flashvars[leadWithHTML5]=true&amp;flashvars[sideBarContainer.plugin]=true&amp;flashvars[sideBarContainer.position]=left&amp;flashvars[sideBarContainer.clickToClose]=true&amp;flashvars[chapters.plugin]=true&amp;flashvars[chapters.layout]=vertical&amp;flashvars[chapters.thumbnailRotator]=false&amp;flashvars[streamSelector.plugin]=true&amp;flashvars[EmbedPlayer.SpinnerTarget]=videoHolder&amp;flashvars[dualScreen.plugin]=true&amp;flashvars[Kaltura.addCrossoriginToIframe]=true&amp;&wid=0_0rywna7h" width="1000" height="500" allowfullscreen webkitallowfullscreen mozAllowFullScreen allow="autoplay *; fullscreen *; encrypted-media *" sandbox="allow-forms allow-same-origin allow-scripts allow-top-navigation allow-pointer-lock allow-popups allow-modals allow-orientation-lock allow-popups-to-escape-sandbox allow-presentation allow-top-navigation-by-user-activation" frameborder="0" title="Kaltura Player"></iframe>


Go to 'https://aws.amazon.com' in a Web browser.

Select 'My Account' menu option 'AWS Management Console."

![AWS Management Console](../../images/aws_1.PNG)

Log in with your username & password.

!!! Note

          If you need to create an account, please follow these instructions. (https://aws.amazon.com/premiumsupport/knowledge-center/create-and-activate-aws-account/)

Check to make sure the location corresponds your location is correct. 

For me it is US West (N. California)


![AWS Dashboard](../../images/aws_2.PNG)

Click on Services (upper left).


![AWS Services](../../images/aws_3.png)

Click on EC2.

![EC2](../../images/aws_4.png)

Click on Launch Instance

!!! Note
         Amazon Elastic Cloud Computing features virtual computing environments called instances. These instances can vary in configurations of CPU, memory, storage, networking capacity.
         For the purposes of future tutorials, we will launch Ubuntu 20.04 Pro LTS. LTS releases are the ‘enterprise grade’ releases of Ubuntu and are utilised the most. 


![Launch Instance](../../images/aws_5.png)



Select AWS Marketplace on the left hand side



![AWS Marketplace](../../images/aws_6.png)

Type Ubuntu 20.04 Pro LTS on the search bar underneath

Step 1: Choose an Amazon Machine Image (AMI)


An Amazon Machine Image is a special type of virtual appliance that is used to create a virtual machine within the Amazon Elastic Compute Cloud. It is a template for the root volume of an instance (operating system, application server, and applications.)

Select Ubuntu Server 20.04 Pro LTS (HVM)

!!! Note
  	The Free Tier Eligible tag lets us know that this particular operating system is covered by the Free Tier program where you use the service with limits without being charged. Limits could be based on how much storage you can use and how many hours you can use in one month.



![AMI](../../images/aws_7.png)

Step 2: Choose an Instance Type

Amazon EC2 provides a wide selection of instance types optimized to fit different use cases. Instances are virtual servers that can run applications. They have varying combinations of CPU, memory, storage, and networking capacity, and give you the flexibility to choose the appropriate mix of resources for your applications. Learn more about instance types and how they can meet your computing needs.


Select t2 micro

![t2.micro](../../images/aws_8.png)

Click Continue

![Ubuntu Pro](../..images/aws_9.png)

Steps 3: Configure Instance

Configure the instance to suit your requirements. You can launch multiple instances from the same AMI, request Spot instances to take advantage of the lower pricing, assign an access management role to the instance, and more.

!!! Note
	A Spot Instance is an unused EC2 instance that is available for less than the On-Demand price. Because Spot Instances enable you to request unused EC2 instances at steep discounts, you can lower your Amazon EC2 costs significantly.

Step 4: Add Storage

Your instance will be launched with the following storage device settings. You can attach additional EBS volumes and instance store volumes to your instance, or edit the settings of the root volume. You can also attach additional EBS volumes after launching an instance, but not instance store volumes. Learn more about storage options in Amazon EC2.

Step 5: Add Tags

A tag consists of a case-sensitive key-value pair. For example, you could define a tag with key = Name and value = Webserver.
A copy of a tag can be applied to volumes, instances or both.
Tags will be applied to all instances and volumes. Learn more about tagging your Amazon EC2 resources.

Step 6: Configure Security Groups

### First Time Through (generate a new key pair)

If you do not have any key pairs, enter a key pair name and then download a key pair. Then click Launch Instance

!!! Note

	**Why do I need a key pair?**
    Easy Answer: *Security*. 

    The SSH (Secure Shell) protocol uses encryption to secure the connection between a client and a server. All user authentication, commands, output, and file transfers are encrypted to protect against attacks in the network. With SSH protocol (secure Shell) public key authenticantion improves security as it frees users from remembering complicated passwords.

![mobaxterm](../../images/aws_10.PNG)

### Next times through (select an existing key pair)

Select a key pair and click 'Launch Instances'.

## Click on View Instances

* After you click 'Launch Instance', you should see this:

![SSH](../../images/aws_11.PNG)

* Click on this first hyperlink: i-038c58bfbe9612c57. Your hyperlink may not be exactly the same. 

## Select the public DNS name for later use

Highlight and copy the Public DNS address and save for future steps.

![Remote Host](../../images/aws_12.PNG)

# Next Steps

## Logging into your new instance "in the cloud" (Windows Version)

Ok, so you've created a running computer. How do you get to it?

The main thing you'll need is the network name of your new computer. To retrieve this, go to the instance view and click on the instance, and find the "Public DNS." This is the public name of your computer on the Internet.

Copy this name, and connect to that computer with ssh under the username 'ubuntu,' as follows.

## Install mobaxterm

MobaXterm is a terminal for Windows with an X11 server, a tabbed SSH client and several other network tools for remote computing (VNC, RDP, telnet,rlogin)
MobaXterm brings all the essential Unix commands to Windows desktop, in a single portable exe file which works out of the box.

First, download mobaxterm and run it.

Download Link: https://download.mobatek.net/2032020060430358/MobaXterm_Portable_v20.3.zip

## Start a new session

![mobaxterm1](../../images/mobaxterm_1.PNG)

## Fill in session settings

Put in your hostname (should be ec2-XXX-YYY-AAA.compute-1.amazon.aws.com), select 'specify username', and enter 'ubuntu'

![Remote Host](../../images/mobaxterm_2.PNG)

## Specify the session key

Copy the downloaded.pem file into your primary hard disk (generally C:) and then put the full path into it.

![Private Key](../../images/mobaxterm_3.PNG)

## Click OK

If you see this screen and the line ubuntu@ip-###-##-#-##:~$,
this means your instance computer is ready.

![Ubuntu Terminal](../../images/mobaxterm_3.PNG)

## Terminating the Instance

Once you have completed your tasks and are sure you do not need the instance any longer, you may terminate the instance by returning to AWS Management Console

- Click on Services 
- Click EC2 
- Click Instance on the left hand side bar menu and it should bring you to the list of running instances on your account.
- Click on the instance you would like to terminate
- Click actions
- Click instance state
- Select terminate



![Terminate](../../images/Terminate.png)

!!! Warning
        
 	Terminating an instance will erase all the work you have saved on the instance.


















