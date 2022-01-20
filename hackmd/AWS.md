# A Hands On Introduction to AWS:

**When:**  

**Instructors:**  

**Moderator:**  

**Helpers:**  


:speech_balloon: Please fill out our [pre-workshop survey](https://docs.google.com/forms/d/e/1FAIpQLSd7SuaESZ0g8MZswdgxlnjC6m8r5fQUcyH6QYMzN2RW4mjCew/viewform)!


## Description

This 2 hour hands-on tutorial will introduce you to creating a computer "in the cloud" and logging into it, via Amazon Web Services. We'll create a small general-purpose Linux computer, connect to it, and run a small job while discussing the concepts and technologies involved.

:computer: Links for moderator to share with participants:

- [x] [Pre-workshop survey](https://docs.google.com/forms/d/e/1FAIpQLSd7SuaESZ0g8MZswdgxlnjC6m8r5fQUcyH6QYMzN2RW4mjCew/viewform)
- [x] [Workshop notes link](https://hackmd.io/4zmi-AQBQIelVp3_WFaxyQ?view)
- [x] [Post-workshop survey](https://docs.google.com/forms/d/e/1FAIpQLSe2qR2t3PIJlRWoOERPK-lRz3JvQv6z_nHGCCQ1771-EeZ0EQ/viewform)

## Hello!

Your instructors are part of the training and engagement team for the [NIH Common Fund Data Ecosystem](https://nih-cfde.org/), a project supported by the NIH to increase data reuse and cloud computing for biomedical research.


This is only our second workshop for our Amazon Web services lesson :slightly_smiling_face:	, and we have the following goals:
* run the material by you all!
* gather questions and refine the tutorial materials!
* help you think about if and how to use cloud computers for your work!

So, please ask lots of questions, and even the ones we can't answer yet we'll figure out for you!

### Costs and payment

Today, everything you do will be paid for by us (well, the NIH). In the future, if you create your own AWS account, you'll have to put your own credit card on it. We'd be happy to answer questions about how to pay for AWS - we've used invoicing as well as credit cards.

:smiley_cat: **Your free login credentials will work for the next 24 hours**

### Workshop structure and plan

* Brief introduction to AWS and the cloud
* An instructor will do a demo, taking questions
* Everyone will do it themselves!
* Install and run things in the cloud computer
* Take your questions

### How to ask questions
If you have questions at any point, 
- Drop them in the chat, or
- Unmute yourself and ask during the workshop, or
- Direct messages to the moderator, Saranya Canchi, are also welcome.

We're going to use the "raise hand" reaction in zoom to make sure people are on board during the hands-on activities. So it's probably not the best way to get our attention if you have questions.

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
- Best suited for analyses that could crash your local computer. E.g. those that generates or uses large output files or takes too long


### Advantages of using AWS

- Sign up process is relatively easy (you need a credit card and some patience to deal with delays in two-factor authentication)
- Simple billing
- Stable services with only 3-4 major outages that only lasted 2-3 hours and did not affect all customers (region-specific). A large team of employees who are on top of any problems that arise!
- Lots of people use it, so there are a ton of resources
- Spot instances (unused EC2 instances) - you can bid for a price. It is cheap, but your services might be terminated if someone outbids you. 



## Let's get started!


We will create a cloud computer - an "instance" - and then log in to it.


**Log in at**: https://cfde-training-workshop.signin.aws.amazon.com/console

Use your registration e-mail (see bottom of this page if you forgot!) and password `CFDErocks!`

### "Spinning up" instances

Checklist:
- [ ] Demo only, then hands-on walk-through
- [ ] Select a region: geographic area where AWS has data centers
- [ ] New EC2 experience vs. old
- [ ] Pick the AMI (OS)
- [ ] Pick an instance (T2 micro free tier!) 
- [ ] Launch

[Link to tutorial](https://training.nih-cfde.org/en/latest/Bioinformatics-Skills/Introduction_to_Amazon_Web_Services/introtoaws3/)

### Connecting to instances

- [ ] Connect to the instance via the web browser
- [Link to tutorial pages](https://hackmd.io/_r6ENz34QbWQLeKah0D6Ug?view)

Other ways to connect to the instance:
- [Connect to an Instance for Windows Users using MobaxTerm](https://training.nih-cfde.org/en/latest/Bioinformatics-Skills/Introduction_to_Amazon_Web_Services/introtoaws4/)
- [Connect to an Instance for Mac Users using MacOS Terminal](https://training.nih-cfde.org/en/latest/Bioinformatics-Skills/Introduction_to_Amazon_Web_Services/introtoaws4/)


## Installing programs and running them in the cloud

- Install simple Bioinformatic Software (FASTQC)
- Download FastQ (raw RNA Sequence) data
- Brief Overview of the FASTQC HTML report
- Demo - How to Terminate an Instance

FASTQC Documentation: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
https://hackmd.io/4zmi-AQBQIelVp3_WFaxyQ?view#What-is-FASTQC

### What is FASTQC?

FastQC aims to provide a simple way to do some quality control checks on raw sequence data coming from high throughput sequencing pipelines. 

- It provides a modular set of analyses which you can use to give a quick impression of whether your data has any problems of which you should be aware before doing any further analysis.
- The aim of this tool is to spot issues that originate from the sequencer or in the starting library material.
- Output of fastqc is an HTML based permanent report



FastQC Functions include:

- Import of Data from BAM, SAM or FastQ files (any variant)
- Providing a quick overview to tell you in which areas of your data there may be problems
- Export of results to an HTML based permanent report
- Offline operation to allow automated generation of reports without running the interactive application



Commands to Run (Explain commands)

```
# Update system packages
sudo apt update

# make a directory fastqc
mkdir fastq

# change into the directory fastqc
cd fastq

# install fastqc
sudo apt install fastqc -y

# download a FASTQ data file from osf.io
curl -L https://osf.io/8rvh5/download -o ERR458494.fastq.gz

# run fastqc on the dowloaded file
fastqc ERR458494.fastq.gz

# view files
ls
```

Additional details or individual commands:

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
- equivalent to making a new folder in  Windows

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

***




## Shutting down instances

When you shutdown your instance any data that is on a non-persistent disk goes away permanently. But you also stop being charged for any compute and data, too!

:bulb: **Stopping vs hibernation vs termination**

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

- [ ] A little bit about AWS and cloud computing
- [ ] How to launch an instance 
- [ ] How to connect to the instance
- [ ] How to install and run a software program on the instance 


## Post-workshop survey


[Please fill out our post workshop survey!](https://forms.gle/RUELsnBPnCM6AbTWA)

## Questions and comments?

We'll send around this link via e-mail -- please do fill it out, thank you!


## Appendix 

### Additional Resources

- Understanding data transfer costs in AWS: https://github.com/open-guides/og-aws#aws-data-transfer-costs
- Useful tips: https://wblinks.com/notes/aws-tips-i-wish-id-known-before-i-started/
- Consolidated billing: https://docs.aws.amazon.com/awsaccountbilling/latest/aboutv2/consolidated-billing.html
