# A Hands On Introduction to AWS

**When:**  

**Instructors:**  

**Moderator:**  

**Helpers:**  

 - [ ] Please fill out our [pre-workshop survey](https://forms.gle/TNhgV1KYCHgRKr1v6)!


## Description

This 2 hour hands-on tutorial will introduce you to creating a computer "in the cloud" and logging into it, via Amazon Web Services. We'll create a small general-purpose Linux computer, connect to it, and run a small job while discussing the concepts and technologies involved.

Links for moderator to share with participants:

- [ ] [Pre-workshop survey](https://forms.gle/TNhgV1KYCHgRKr1v6)
- [ ] [Workshop notes link](https://hackmd.io/jTqsKMSZRU-L6NkJM4fZew?view)
- [ ] [mobaXterm install link](https://mobaxterm.mobatek.net/)
- [ ] [AWS lessons](https://training.nih-cfde.org/en/latest/Bioinformatics-Skills/Introduction_to_Amazon_Web_Services/introtoaws1/) 
- [ ] [Post-workshop survey](https://forms.gle/2tGthbddRCQ72vQG9)


## Hello!

Your instructors are part of the training and engagement team for the [NIH Common Fund Data Ecosystem](https://nih-cfde.org/), a project supported by the NIH to increase data reuse and cloud computing for biomedical research.

This is our third workshop for our Amazon Web services lesson :slightly_smiling_face:, and we have the following goals:
* run the material by you all!
* gather questions and refine the tutorial materials!
* help you think about if and how to use cloud computers for your work!

So, please ask lots of questions, and even the ones we can't answer yet we'll figure out for you!

### Costs and payment

Today, everything you do will be paid for by us (well, the NIH). In the future, if you create your own AWS account, you'll have to put your own credit card on it. We'd be happy to answer questions about how to pay for AWS - we've used invoicing as well as credit cards.


 _**Your free login credentials will work for the next 24 hours !**_


### Workshop structure and plan

* Brief introduction to AWS and the cloud
* Set up an instance and connect to it
* Install and run things in the cloud computer
* Learn how to download output files to local machine
* Take your questions

### How to ask questions
If you have questions at any point, 
- Drop them in the chat, or
- Direct messages to the moderator (Marisa Lim) are welcome, or
- Unmute yourself and ask during the workshop, or
- Use the raise hand reaction in zoom

We're going to use the green check mark reaction in zoom to make sure people are on board during the hands-on activities.

## Some background

What is cloud computing? 
- Renting and use of IT services over the internet.
- Amazon and Google, among others, rent compute resources over the internet for money.

Why might you want to use a cloud computer?

There are lots of reasons, but basically "you need a kind of compute or network access that you don't have."

- More memory than you have available otherwise
- An operating system you don't have access to (Windows? Mac?)
- Installation privileges for software
- May not want to install brand new software on your local computer


## Amazon, terminology, and logging in!

- Amazon web services is one of the most broadly adopted cloud platforms
- It is a hosting provider that gives you a lot of services including cloud storage and cloud compute. 


Terminology:
* Instance - a computer that is running ...somewhere in "the cloud". The important thing is that someone else is worrying about the hardware etc, so you're just renting what you need!
* Cloud computer - same as an "instance".
* Image - the basic computer install from which an instance is constructed. The configuration of your instance at launch is a copy of the Amazon Machine Image (AMI)
* [EC2](https://en.wikipedia.org/wiki/Elasticity_(cloud_computing)#:~:text=In%20cloud%20computing%2C%20elasticity%20is,demand%20as%20closely%20as%20possible%22.) - elastic compute cloud.


**Amazon's main compute rental service is called Elastic Compute Cloud (or EC2) and that's what we'll be showing you today.**



### EC2

- Amazon Elastic Compute Cloud (Amazon EC2) is a web service that provides secure, resizable compute capacity in the cloud.
- Basically, you rent virtual computers that are configured according to your needs and run applications and analyses on that computer. 
- Best suited for analyses that could crash your local computer. E.g. those that generate or use large output files or take too long


### Advantages of using AWS

- Sign up process is relatively easy (you need a credit card and some patience to deal with delays in two-factor authentication)
- Simple billing
- Stable services with only 3-4 major outages that only lasted 2-3 hours and did not affect all customers (region-specific). A large team of employees who are on top of any problems that arise!
- Lots of people use it, so there are a ton of resources
- Spot instances (unused EC2 instances) - you can bid for a price. It is cheap, but your services might be terminated if someone outbids you. 

---

## Let's get started!


We will create a cloud computer - an "instance" - and then log in to it.

**Log in at**: https://cfde-training-workshop.signin.aws.amazon.com/console

Use your registration e-mail (see bottom of this page if you forgot!) and password `CFDErocks!`

### "Spinning up" instances

Checklist for hands-on walk-through
- [ ] Select a region: geographic area where AWS has data centers
- [ ] Pick the AMI (OS)
- [ ] Pick an instance (T2 micro free tier!) 
- [ ] Edit security groups
- [ ] Launch

[Link to tutorial](https://training.nih-cfde.org/en/latest/Bioinformatics-Skills/Introduction_to_Amazon_Web_Services/introtoaws3/)

### Connecting to instances

- [ ] [Connect to the instance via the web browser](https://training.nih-cfde.org/en/latest/Bioinformatics-Skills/Introduction_to_Amazon_Web_Services/introtoaws4/)

Other ways to connect to the instance:

We have tutorials on connecting to an instance for **Windows** Users using MobaXterm and for **Mac Users** using MacOS Terminal. Please visit our ["Connect to an Instance"](https://training.nih-cfde.org/en/latest/Bioinformatics-Skills/Introduction_to_Amazon_Web_Services/introtoaws4/) webpage and select your OS using the tabs on the top of the page.


## Installing programs and running them in the cloud

- Install a simple bioinformatic software (FastQC)
- Download fastq (raw RNA Sequence) data
- Brief overview of the FastQC HTML report
- Demo - How to terminate an instance

FastQC Documentation: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
https://hackmd.io/4zmi-AQBQIelVp3_WFaxyQ?view#What-is-FASTQC

### What is FastQC?

FastQC aims to provide a simple way to do some quality control checks on raw sequence data coming from high throughput sequencing pipelines. 

- It provides a modular set of analyses which you can use to give a quick impression of whether your data has any problems of which you should be aware before doing any further analysis.
- The aim of this tool is to spot issues that originate from the sequencer or in the starting library material.
- Output of fastqc is an HTML based permanent report

FastQC functions include:

- Import of Data from BAM, SAM or FastQ files (any variant)
- Providing a quick overview to tell you in which areas of your data there may be problems
- Export of results to an HTML based permanent report
- Offline operation to allow automated generation of reports without running the interactive application



### Commands to run 
(explain commands)


1) Update system packages:

```
sudo apt update
```

2) Make a directory
```
mkdir fastq
```

3) Change into the directory
```
cd fastq
```

4) Download a fastq data file from osf.io
```
curl -L https://osf.io/8rvh5/download -o ERR458494.fastq.gz
```

5) Check if your file has been downloaded
```
ls -l
```

6) Install FastQC
```
sudo apt install fastqc -y
```

7) Run FastQC on the dowloaded file
```
fastqc ERR458494.fastq.gz
```

8) view files
```
ls
```

Learn more about the individual commands :

```apt-cache search [search term 1]```
- search available software for installation

```sudo apt update```
- download packaged information from all configured sources from the internet
- This will update the package lists from all repositories in one go. Remember to do this after every added repository!

```sudo apt install <program1> -y```
- install package
- Other programs ("ncbi blast+", )

```mkdir <directory name>```
- make a new directory
- equivalent to making a new folder in Windows

```cd <directory name>```
- change directory
- equivalent to double clicking a folder

```curl -L <url> <filename> -o <file.html>```

- curl stands for "Client URL"
    - transfers data to or from a network server
    - "-L" or location 
    - "-o" output

```fastqc ERR458494.fastq.gz```
- Run FastQC on ERR458494.fastq.gz
- ERR458494.fastq.gz - "Yeast" Sample

***

Analysis Modules Documentation: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/

What a good data file looks like
https://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html

What bad data looks like
https://www.bioinformatics.babraham.ac.uk/projects/fastqc/bad_sequence_fastqc.html

Video Walkthrough:

[FastQC tool for read data quality evaluation](https://www.youtube.com/watch?v=lUk5Ju3vCDM)

[Using FastQC to check the quality of high throughput sequence](https://www.youtube.com/watch?v=bz93ReOv87Y&t=116s)
:::

***

## Using `screen`

So far in this workshop, we have only encountered programs that install quickly. The analysis we ran was also pretty quick because we only ran it on one file! 

In your own work, you may encounter programs that have lengthy installations, and/or you may need to analyze a large number of files.

While performing a long-running task on a remote machine, a sudden drop in your internet connection would terminate the SSH session and your work would be lost!

The `screen` utility provides a work-around to this problem. `screen` is a terminal multiplexer i.e. you can open many virtual terminals. Processes running in `screen` will continue to run even when the terminal is not visible, or if you get disconnected from the internet.

### Commands to run 

1) Install `screen`:

```
sudo apt-get install screen
```

2) Running `screen`
```
screen
```
Press space (twice) or enter to get the command prompt

3) Run a code inside `screen` session
```
top
```

4) Detaching screen
```
Press ctrl + a + d keys
```

5) List screen sessions
```
screen -ls
```

6) Reattach screen
```
screen -r <screen_ID>
```

7) Repeat step 4 to detach



## Downloading data from AWS instance onto local computer

### WindowsOS

#### MobaXterm installation

1. Go to the MobaXterm website to [download](https://mobaxterm.mobatek.net/)
2. Click on "GET MOBAXTERM NOW!"
3. The Home Edition is perfect for normal use and it is free! Click "Download now"
4. Click on "MobaXterm Home Edition v20.6 (Portable edition)" and save as in your Downloads folder
5. Go to Downloads folder, click on the zipped folder, click "Extract all", click "Extract"
6. The MobaXterm application is now in the unzipped folder
7. Click on the MobaXterm application to open it!

#### Connecting to instance
1. Go back to your [instance page](https://us-west-1.console.aws.amazon.com/ec2/v2/home?region=us-west-1#Instances:), select it and click on "Connect". The Public DNS information you need to connect to your instance via ssh can be found in the "SSH client" tab:

![](https://i.imgur.com/EilADhq.png)

2. In MobaXterm, click on "Session"
3. Click on "SSH"
4. Enter the Public DNS as the "Remote host"
5. Check box next to "Specify username" and enter "ubuntu" as the username
6. Click the "Advanced SSH settings" tab
7. Check box by "Use private key"
8. Use the document icon to navigate to where you saved the private key (e.g., "amazon.pem") from AWS on your computer. It is likely on your Desktop or Downloads folder
9. Click "OK"
10. A terminal session should open up with a left-side panel showing the file system of our AWS instance! You can click on the FastQC html file and view in browser to open. There are also options in the panel to download files.


#### MacOS

- Start Terminal 
- Change the permissions on the .pem file for security purposes (removes read, write, and execute permissions for all users except the owner (you)
```
chmod og-rwx ~/Desktop/amzon.pem
```

- Change directory to Desktop. Your `.pem` file is on your Desktop

```
cd ~/Desktop
```
Go back to your [instance page](https://us-west-1.console.aws.amazon.com/ec2/v2/home?region=us-west-1#Instances:), select it and click on "Connect". The information you need to connect to your instance via ssh can be found in the "SSH client" tab:

![](https://i.imgur.com/EilADhq.png)


- Use the `scp` command on your local terminal to copy your `.html` file!

```
scp -i <your-.pem> ubuntu@???-??-??-???-??.us-west-1.compute.amazonaws.com:/home/ubuntu/fastq/ERR458494_fastqc.html ./
```

Don't forget to change the stuff after `ubuntu@` to match your instance!


## Shutting down instances

When you shutdown your instance any data that is on a non-persistent disk goes away permanently. But you also stop being charged for any compute and data, too!

**Stopping vs hibernation vs termination**

- Stopping: 
    - saves data to EBS root volume 
    - only EBS data storage charges apply 
    - No data transfer charges or instance usage charges 
    - RAM contents not stored

- Hibernation: 
    - charged for storage of any EBS volumes 
    - stores the RAM contents 
    - it's like closing the lid of your laptop

- Termination: 
    - complete shutdown 
    - EBS volume is detached 
    - data stored in EBS root volume is lost forever
    - instance cannot be relaunched





## Checklist of things you learned today!

- [x] A little bit about AWS and cloud computing
- [x] How to launch an instance 
- [x] How to connect to the instance
- [x] How to install and run a software program on the instance 
- [x] How to terminate your instance 


## Post-workshop survey


[Please fill out our post workshop survey!](https://forms.gle/2tGthbddRCQ72vQG9)


## Questions and comments?

We'll send around this link via e-mail -- please do fill it out, thank you!

## Upcoming CFDE workshops
- [Intro to Conda](https://registration.genomecenter.ucdavis.edu/events/intro_to_conda_march4/)


## Appendix 

### Additional Resources

- Understanding data transfer costs in AWS: https://github.com/open-guides/og-aws#aws-data-transfer-costs
- Useful tips: https://wblinks.com/notes/aws-tips-i-wish-id-known-before-i-started/
- Consolidated billing: https://docs.aws.amazon.com/awsaccountbilling/latest/aboutv2/consolidated-billing.html


## FAQs

**What are the advantages of using AWS over an academic HPC?**
- Most universities don't have a HPC
- No queues!
- Can set up as many instances as you want (as long as you are willing to pay for it)
- Can install anything without needing admin permissions
- Almost no scheduled or unscheduled outages
- Easier to set up 
- Easier to learn and get help on the internet
- Costs more over time, but someone is paying for the HPC too! 

*But if you have a good HPC, please use it!*

**Can you set up multiple instances at once**
- Yes!
- There is a limit per account but it is a very large number and doesn't apply to most people

**Can you launch more than one instance with the same configurations?**
- Yes, there is an option to do this on the instance set up page.
- Look in the second tab!

**Can you copy an instance or share an instance with collaborators?** 
- Yes, but this is not as straightforward as it seems.
- The way to clone an instance is via [snapshots](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ebs-creating-snapshot.html)

