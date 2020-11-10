# Setting up an AWS instance

To start, we will learn how to set-up Ubuntu 20.04 Pro LTS open source software operating system. Ubuntu 20.04 Pro LTS is one of the programs offered in the Amazon Free Tier as well as one of the most popular open source operating systems.

Follow along with these steps and/or watch our [walk-through tutorial](./introtoaws2.md) to get started!

### Step 1: Log in to an AWS account

Go to [Amazon Web Services](https://aws.amazon.com) in a web browser. Select the "My Account" menu option "AWS Management Console". Log in with your username & password.

![AWS Management Console](../../images/aws_1.PNG "AWS my account button")

!!! Note

          If you need to create an account, please follow the [AWS instructions for creating an account](https://aws.amazon.com/premiumsupport/knowledge-center/create-and-activate-aws-account/).

!!! Warning

    If you are creating a new account, it could take up to 24 hours to be activated. You'll need a credit card to set up the account.

### Step 2: Choose virtual machine

For this tutorial, select the AWS region that is closest to your current geographic location. The AWS region of your remote machine is displayed on the top right of this page. Click on it and choose the location that best describes the region you are currently located.

![AWS Dashboard](../../images/aws_2.PNG "AWS amazon machine selection")

!!! note "AWS Region"

    The default region is automatically displayed in the AWS Dashboard. The [choice of region](https://docs.aws.amazon.com/emr/latest/ManagementGuide/emr-plan-region.html) has implications on fees, speed, and performance.

Click on "Services" (upper left):

![AWS Services](../../images/aws_3.png "AWS Services button")

Click on "EC2":

![EC2](../../images/aws_4.png "AWS EC2 button")

!!! Note

         Amazon Elastic Cloud Computing features virtual computing environments called instances. These instances can vary in configurations of CPU, memory, storage, networking capacity. For the purposes of future tutorials, we will launch Ubuntu 20.04 Pro LTS. LTS releases are the ‘enterprise grade’ releases of Ubuntu and are utilized the most.

Click on "Launch Instance":

![Launch Instance](../../images/aws_5.png "AWS launch button")

Select "AWS Marketplace" on the left hand side tab:

![AWS Marketplace](../../images/aws_6.png "AWS marketplace button")

### Step 3: Choose an Amazon Machine Image (AMI)

An Amazon Machine Image is a special type of virtual appliance that is used to create a virtual machine within the Amazon Elastic Compute Cloud. It is a template for the root volume of an instance (operating system, application server, and applications).

Type `Ubuntu 20.04 Pro LTS` in the search bar. Click "Select":

![AMI](../../images/aws_7.png "AWS Ubuntu AMI")

Click "Continue":

![Ubuntu Pro](../../images/aws_9.PNG "Ubuntu Pro information")

### Step 4: Choose an instance type

Amazon EC2 provides a wide selection of instance types optimized to fit different use cases. Instances are virtual servers that can run applications. They have varying combinations of CPU, memory, storage, and networking capacity, and give you the flexibility to choose the appropriate mix of resources for your applications. Learn more about instance types and how they can meet your computing needs.

Select the row with `t2.micro`, the free tier eligible option:

![t2.micro](../../images/aws_8.png "t2 micro instance type")

!!! Note

    The Free Tier Eligible tag lets us know that this particular operating system is covered by the Free Tier program where you use (limited) services without being charged. Limits could be based on how much storage you have access to and/or how many hours of compute you can perform in a one month.

### Step 5: Optional Configurations

There are several optional set up configurations. You can either click "Review and Launch" now to start the instance we've configured thus far in the tutorial without these additional configurations or as necessary, click on the following tabs to continue configuring. Start the first option by clicking "Next: Configure Instance Details" on the AWS page.

=== "Configure Instance Details"

    Configure the instance to suit your requirements. You can launch multiple instances from the same AMI, request Spot instances to take advantage of the lower pricing, assign an access management role to the instance, and more.

    A Spot Instance is an unused EC2 instance that is available for less than the On-Demand price. Because Spot Instances enable you to request unused EC2 instances at steep discounts, you can lower your Amazon EC2 costs significantly.

=== "Add Storage"

    Your instance will be launched with the following storage device settings. You can attach additional EBS volumes and instance store volumes to your instance, or edit the settings of the root volume. You can also attach additional EBS volumes after launching an instance, but not instance store volumes. Learn more about storage options in Amazon EC2.

=== "Add Tags"

    A tag consists of a case-sensitive key-value pair. For example, you could define a tag with key = Name and value = Webserver. A copy of a tag can be applied to volumes, instances or both. Tags will be applied to all instances and volumes. Learn more about tagging your Amazon EC2 resources.

=== "Configure Security Group"

    A security group is a set of firewall rules that control the traffic for your instance. On this page you can add rules to allow specific traffic to reach your instance. You can create a new security group or select from an existing one.

### Step 6: Review and Launch instance

After configuration settings are complete, click "Review and Launch" and "Launch". If you are launching an AWS instance for the first time, you will need to generate a key pair.

Choose the "Create a new key pair" option from the drop down menu. Under key pair name, type "amazon" and click "save". The default location for saving files on a Mac is the "Downloads" folder -- that's where your key pair can be found. Next time you launch an instance, you can reuse the key pair you just generated.

If you have a previously generated key pair, you can reuse it to launch an instance. For this tutorial, we are calling the key pair "amazon.pem".

![mobaxterm](../../images/aws_10.png "key pair set up")

!!! Note "Why do I need a key pair?"

    For security purposes, the SSH (Secure Shell) protocol uses encryption to secure the connection between a client and a server. All user authentication, commands, output, and file transfers are encrypted to protect against attacks in the network. With SSH protocol (secure Shell) public key authentication improves security as it frees users from remembering complicated passwords.

Then select your key pair, check the acknowledgement box, and click "Launch Instance". Now you should see:

![SSH](../../images/aws_11.png "Instance ID link")

Click on this first hyperlink, in the image above, "i-038c58bfbe9612c57". Your hyperlink may be different.

![Remote Host](../../images/aws_12.PNG "AWS instance running page")

This page shows you a list of all your active instances. Users may launch as many instances as they wish. Just remember that every instance costs money if you don't qualify for the Free Tier. On this page, there is a "Public DNS" address, with the format `ec2-XXX-YYY-AAA.compute-1.amazon.aws.com`. You'll need this address to connect to your AWS computer.


You have now successfully launched your AWS instance! You will need the Public DNS address from this amazon webpage to access your AWS computer, so do not close the page yet.

If you happen to close the webpage on accident, [click on this link](https://us-west-1.console.aws.amazon.com/ec2/v2/home?region=us-west-1#Instances:sort=instanceId)

Continue on to the next lesson to learn how to connect to your AWS computer!

