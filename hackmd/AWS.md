# A Hands On Introduction to AWS 

**When:**  

**Instructors:**  

**Helpers:**  

## Description
This 2 hour hands on tutorial will introduce you to creating a computer "in the cloud" and logging into it, via Amazon Web Services. We'll create a small general-purpose Linux computer connect to it, and run a small job, while discussing the concepts and technologies involved.

## Agenda

### Introduction and workshop goals

Hello!
I'm Titus Brown, a professor here at UC Davis. I'm joined today by Dr. Abhijna Parigi. We are both part of the training and engagement team for the [NIH Common Fund Data Ecosystem](https://nih-cfde.org/), a project supported by the NIH to increase data reuse and cloud computing for biomedical research.
You can contact us at ctbrown@ucdavis.edu, and abhijna.parigi@gmail.com.
This is our first workshop for our Amazon Web services lesson, and we have the following goals:
* run the material by you all!
* gather questions and refine the tutorial materials!
* help you think about if and how to use cloud computers for your work!
So, please ask lots of questions, and even the ones we can't answer yet we'll figure out for you!
Today, everything you do will be paid for by us (well, the NIH). In the future, if you create your own AWS account, you'll have to put your own credit card on it. We'd be happy to answer questions about how to pay for AWS - we've used invoicing as well as credit cards.
How we will do this!
* first, an instructor will walk through process, taking questions
* then everyone will do it themselves!
* then an instructor will walk through things again, taking more questions.
* if you have questions, drop them in chat!

### Some background
Amazon and Google, among others, rent compute resources for money.
Why might you want to use a cloud computer?
There are lots of reasons, but it basically boils down to "you need a kind of compute or network access that you don't have."
- GPU nodes
- more memory than you have available otherwise
- an operating system you don't have access to (Windows? Mac?)
- installation privileges for software
Amazon's main compute rental service is called Elastic Cloud Computing (or EC2) and that's what we'll be showing you today.

### Amazon, terminology, and logging in!
What we're going to start with is creating a cloud computer - an "instance" - and then logging in to it.
Terminology:
* instance - a computer that is running ...somewhere in "the cloud". The important thing is that someone else is worrying about the hardware etc, so you're just renting what you need!
* cloud computer - same as an "instance".
* image - the basic computer install from which an instance is constructed
* "EC2" - elastic cloud computer.
Log in at https://cfde-training-workshop.signin.aws.amazon.com/console
Use your registration e-mail and password `CFDErocks!`

### "Spinning up" instances

* Select a region
* New EC2 experince vs. old
* Pick the AMI (OS)
* Pick an instance (T2 micro free tier!) 
* Launch
[Link to tutorial](https://training.nih-cfde.org/en/latest/Bioinformatics-Skills/Introduction_to_Amazon_Web_Services/introtoaws3/)

### Connecting to instances
[Link to tutorial pages](https://hackmd.io/_r6ENz34QbWQLeKah0D6Ug?view)

### Installing software on them and running BLAST
[Link to tutorial pages](https://training.nih-cfde.org/en/latest/Bioinformatics-Skills/Command-Line-BLAST/BLAST1/)
Points to cover:
* you are the sysadmin and can do anything with `sudo`!
* the memory and disk space you have is configured at instances startup
* this is just like any command line instance, BUT no queuing software
* you have to shut it down yourself when you're done.

### Shutting down instances
Points to cover:
* any data that is on a non-persistent disk ...goes away permanently.
* but you also stop being charged for any compute and data, too!

### Questions and comments?
We'll send around this link via e-mail -- please do fill it out, thank you!
[Survey link](https://docs.google.com/forms/d/e/1FAIpQLSeVK_iwccWGSzsMvzsZyHxmkklz1cP5GbsQbEVoVVvLhoeWcQ/viewform)
