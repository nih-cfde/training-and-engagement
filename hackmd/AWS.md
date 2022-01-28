# A Hands-on Introduction to AWS for Cloud Computing

**When:**   January 26, 2022, from 10 am PST - 12 pm PST 

**Instructors:** Dr. Rayna Harris

**Helpers:** Dr. Amanda Charbonneau  Jessica Lumian and Jeremy Walter 

Your instructors are part of the training and engagement team for the [NIH Common Fund Data Ecosystem](https://nih-cfde.org/), a project supported by the NIH to increase data reuse and cloud computing for biomedical research.

#### Description

This 2-hour hands-on tutorial will introduce you to creating a computer "in the cloud" and logging into it, via Amazon Web Services. We will launch a small general-purpose Linux instance, connect to it, and run a small job while discussing the concepts and technologies involved.

#### Overview

[TOC]

#### Before we start

:pencil: Please fill out our [pre-workshop survey](https://forms.gle/hCtbXrDyH9omAYyB6) if you have not already done so! 

:heavy_check_mark: Windows users should install Mobaxterm. Read our [quick installation guide.](https://hackmd.io/YivVGhznRwqVO3fxdGpkOQ?view)

#### Questions?

If you have questions,  
- Type them in the group chat
- Direct message the moderator 
- Unmute and ask them outloud

We're going to use the raised hand :raised_hand: reaction in zoom to make sure people are on board during the hands-on activities.

## 1. Terminology and Sign-On

Cloud computing is the on-demand use of data storage and compute power without direct active management by the user. Amazon Web Services (AWS) is one of the most broadly adopted cloud platforms.

Some advantages of using AWS include:

- Easy sign-on 
- Simple billing
- Stable services
- Customizable images
- Customer support
- Online resources


![](https://uploads-ssl.webflow.com/5e1f17bab0dc6527c1ecc801/5e55f0ab6725fd082d2ea435_amazon-hosting.jpeg)

Amazon's Elastic Compute Cloud (**EC2**) is a web service that provides secure, resizable compute capacity in the cloud. Amazon's Simple Storage Service (**S3**) is widely used for storing and sharing data. 

An **instance** is a virtual machine that runs in the cloud. An **image** (or AMI for Amazon Machine Image) is a template that contains the software configuration (including operating system and applications) required to launch your instance. You can select an image provided by the AWS Marketplace, the AWS community, or you can select one of your own images. When you launch an instance, you specify the type of image to use. 

Today, everything you do will be paid for by us. Your free login credentials will work for the next 24 hours. In the future, if you create an AWS account, you will have to add a credit card for billing. We'd be happy to answer questions about how to pay for AWS.

Log in to your account by going to this web address:  https://cfde-training-workshop.signin.aws.amazon.com/console. 

![](https://hackmd.io/_uploads/SJfyT66pt.png)

Find your first name in the table below and log in with that as your IAM user name and the password `cfde2022!!` 


|   IAM |   IAM  |   IAM  |   IAM  |   Password  |
|---|---|---|---|---|
|   Abdullahi  |   Jasleen  |   Minoo  |   Somi  |   `cfde2022!!`  |
|   Brett  |   Jenna  |   Nathan  |   Sophia  |   `cfde2022!!`  |
|   Connor  |   Jesus  |   Nicole  |   Stefan  |   `cfde2022!!`   |
|   Daulet  |   Jiefei  |   Nina  |   Stephen  |   `cfde2022!!`   |
|   David  |   Ketaki  |   Owen  |   Triveni  |   `cfde2022!!`   |
|   Francisco  |   Kirtan   |   Paul  |   Uma  |   `cfde2022!!`   |
|   Harsh  |   Layla  |   Pornlada  |   Vincent  |   `cfde2022!!`   |
|   Hasim  |   Lessa  |   ramakrishna  |   student2  |  `cfde2022!!`   |
|   Ian  |   Michaelangelo  |   Sai  |   student3  |   `cfde2022!!`   |
|   Jackie  |   Minji  |   sajjad  |   student4  |   `cfde2022!!`   |


> :raised_hand: Raise your hand in Zoom when you've successfully logged in with the workshop user credentials.

## 2. Launching an EC2 Instance 

You can launch an instance using the AWS launch instance wizard. The launch instance wizard specifies all the launch parameters required for launching an instance. Where the launch instance wizard provides a default value, you can accept the default or specify your own value. At the very least, you need to select an AMI and a key pair to launch an instance. Let's walk through the following steps.

1. Open the Amazon EC2 console at https://console.aws.amazon.com/ec2/

2. AWS has servers all over the world. In the top right corner, click the drop-down menu to select a global region. For this workshop choose **US West (N. California) us-west-1**. _In the future, you should pick a region near your or one that contains your data._

![](https://hackmd.io/_uploads/Hk6dhoapY.png)

3. Now, click the  [![  - Launch instances](https://img.shields.io/badge/_-Launch_instances-ec7211)](https://us-west-1.console.aws.amazon.com/ec2/v2/home?region=us-west-1#LaunchInstanceWizard:) button. 


<!---
You should see a page that looks like this:
--->


4. AWS is beta testing a new version of the launch version of the wizard that goes through all the steps in one page instead of many. It is awesome. Click the [![  - Try it now!](https://img.shields.io/badge/_-Try_it_now!-276fc4)](https://us-west-1.console.aws.amazon.com/ec2/v2/home?region=us-west-1#LaunchInstances:) button at the top to get started. _If you accidentally close the banner with the beta button, refresh the page to bring it up again._

![](https://hackmd.io/_uploads/BkDh0opaF.png)

5. First, give your instance a name (such as your first name) so that you can distinguish your instances from your classmates'. This is optional but very useful for keeping track of multiple instances on the same account.  

6. The next step is to pick an image. Our preferred image is not listed in the Quick Start list, so we must find it in the Marketplace. Type **Ubuntu 20.04 LTS - Focal** in the search bar. Then click AWS Marketplace AMIs. Once you see, Ubuntu 20.04 LTS - Focal, click [![Select](https://img.shields.io/badge/Select-ec7211)](https://).  

7. Next, we must specify how much memory and ram we need by specifying an **instance type**. The **t2.micro** instance is "Free tier eligible" and provides 1CPU and 1GB of memory. This is perfect for our class. 

9. The final step is to **create a new key pair**. This will be used in the next section to connect to your instance via `ssh`. Give your key pair a name (without spaces). Use the default settings of RSA type and .pem format. Save this file locally (e.g. in your downloads or your desktop).  

10. For this workshop, we will choose the default network, security, and storage settings, so there is nothing else to change. 

11. Scroll down to the bottom of the page and click [![Launch instance](https://img.shields.io/badge/Launch_instance-ec7211)](https://). 

12. Once your instance launches, click the [![View all instances](https://img.shields.io/badge/View_all_instances-ec7211)](https://) button at the bottom of the page.

> :raised_hand: Raise your hand in Zoom when you've successfully launched an instance.

Congratulations! You have successfully launched an instance. The next step is to connect to your instance.  

## 3. Connecting to AWS instances

There are three ways to connect an AWS instance:
1. with a web browser 
2. using `ssh` from the Terminal
3. using an ssh client such as MobaXterm

Let's connect to our instances using a web browser.

1. Find your instance in the list of running instances. 
2. Click the empty check box next to your name. 
3. Then click "Connect" in the top center of your browser.


![](https://hackmd.io/_uploads/BkcjPnTTt.png)

4. This will open a window that provides details about your instances. Click the [![Connect](https://img.shields.io/badge/Connect-ec7211)](https://) button at the bottom of your screen.

![](https://hackmd.io/_uploads/BkEUdnT6Y.png)

After you click connect, a new tab will open in your browser with a Terminal window that looks something like this. 

![](https://hackmd.io/_uploads/By7fqaaaK.png)

> :raised_hand: Raise your hand in Zoom when you've successfully connected to your instance.

_If at any time, your instance stops responding, hit the "refresh" button and functionality should be restored, right where you left off._

## 4. Running programs at the command line


Now that you have successfully launched a terminal in your browser, you can run programs at the command line. If you attended last week's [Intro to UNIX for Cloud Computing](https://hackmd.io/O2BDsqNWRhSdia26WcCyoQ) workshop, we used a variety of commands to navigate the file system and work with files. Let's revisit a few. 


First, print your working directory with `pwd`.

```
pwd
```
 
 You should see:
 
 
 ```
 /home/ubuntu
 ```


Now if you type `ls`, it may look like you do not have any files, but remember some files are hidden. Let's use the `-a` option to list all files, `-l` for long listing format, and the `-F` option to append a classifier.

```
ls -alF
```

From this, we can see a few directories. The `.ssh` directory contains the ssh key you created a moment ago. 

```
drwxr-xr-x 4 ubuntu ubuntu 4096 Jan 25 18:59 ./
drwxr-xr-x 3 root   root   4096 Jan 25 18:59 ../
-rw-r--r-- 1 ubuntu ubuntu  220 Feb 25  2020 .bash_logout
-rw-r--r-- 1 ubuntu ubuntu 3771 Feb 25  2020 .bashrc
drwx------ 2 ubuntu ubuntu 4096 Jan 25 18:59 .cache/
-rw-r--r-- 1 ubuntu ubuntu  807 Feb 25  2020 .profile
drwx------ 2 ubuntu ubuntu 4096 Jan 25 18:59 .ssh/
```

Your instance comes pre-configured with a number of computer programs. These are stored in your **root directory (`/`)** in the `bin` directory. You can list the programs installed on your instance by providing `ls` with the full path `/bin`. 


```
ls /bin
```

This will print all the installed programs to your screen. Here are a few of the programs. You may recognize a few of the programs we used last week, such as `gunzip`, `gzip`, and `head`. 


```
...
grub-render-label         rvim                               zcat
grub-script-check         savelog                            zcmp
grub-syslinux2cfg         sbattach                           zdiff
gsettings                 sbkeysync                          zdump
gtbl                      sbsiglist                          zegrep
gunzip                    sbsign                             zfgrep
gzexe                     sbvarsign                          zforce
gzip                      sbverify                           zgrep
h2ph                      scp                                zipdetails
h2xs                      screen                             zless
hd                        screendump                         zmore
head                      script                             znew
helpztags                 scriptreplay
hexdump                   scsi_logging_level
...
```

To practice working with some of these command-line programs, we need some files to work on. Let's use the `curl` command to download the same files we used in last week's workshop, which are stored in a .zip file in an Amazon S3 bucket. The `-O` option says to use the same filename that is specified in the web address. 

```
curl -O https://s3.us-west-1.amazonaws.com/dib-training.ucdavis.edu/shell-data2.zip
```

You should see the following message. 

```
% Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
100 38.4M  100 38.4M    0     0  67.8M      0 --:--:-- --:--:-- --:--:-- 67.8M
```

You can type `ls` to confirm the download was completed successfully. 

```
ls
```

You should now have one file in your working directory. 

```
shell-data2.zip
```

Now, we need to uncompress this file with the command `unzip`.


```
unzip shell-data2.zip 
```

However, instead of uncompressing the file, we get the following error and help messages. 

```
Command 'unzip' not found, but can be installed with:

sudo apt install unzip
```

Unfortunately, this tells us that the command we want to use is not installed. Fortunately, it tells us how to install it. Let's try it. 


```
sudo apt install unzip
```

Once the program finishes installing, we can now unzip our files. Remember, you can use the up arrow to scroll through your history to a previous command, then select enter to run it. 


```
unzip shell-data2.zip 
```

After the files are inflated, you can remove .zip file if you wish with the `rm` command.

```
rm shell-data2.zip
```

Now, when we type `ls` we see 6 directories and a README.md file. 

```
ls
```

```
MiSeq  README.md  binder  books  images  seattle  southpark
```

Now that we have these files on a cloud computer, we can run bioinformatic programs that are typically too computationally expensive to run on our local computer. 


#### What is FastQC?

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is a bioinformatic program performs quality control checks on raw sequence data coming from high throughput sequencing pipelines. It runs a set of analyses to help you identify problems in the quality of your samples or sequence. The output of fastqc is an HTML document. 

Here are some links if you would like to learn more about FASTQ on your own time. 

- [Analysis Modules Documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/)
- [What a good data file looks like](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html)
- [What bad data looks like](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/bad_sequence_fastqc.html)
- [Video: FastQC tool for read data quality evaluation](https://www.youtube.com/watch?v=lUk5Ju3vCDM)
- [Video: Using FastQC to check the quality of high throughput sequence](https://www.youtube.com/watch?v=bz93ReOv87Y&t=116s)


To install FASTQC, we first need to update the `apt` package. Then we can use it to install the program. Add the `-y` option to say "yes" to all the installing prompts automatically. 

```
sudo apt update
sudo apt install fastqc -y
```

To double check it was successful, type `fastqc --version`. If it returns 0.11.9, that means installation was successful. You can also type `fastqc --help` to view the manual.


A fastqc command looks like this: `fasqtc -o <output directory> <file>`


The output directory must exist! Before we begin, let's navigate to the `MiSeq` directory with `cd` and create a `results` directory with a subdirectory called `fastqc`. Remember, the `-p` argument will tell `mkdir` to create any missing parent directories.


```
cd MiSeq
mkdir -p results/fastqc
ls
```

We now see the following files and directories. 

```
F3D0_S188_L001_R1_001.fastq    F3D148_S214_L001_R1_001.fastq  F3D7_S195_L001_R1_001.fastq
F3D0_S188_L001_R2_001.fastq    F3D148_S214_L001_R2_001.fastq  F3D7_S195_L001_R2_001.fastq
F3D141_S207_L001_R1_001.fastq  F3D149_S215_L001_R1_001.fastq  F3D8_S196_L001_R1_001.fastq
F3D141_S207_L001_R2_001.fastq  F3D149_S215_L001_R2_001.fastq  F3D8_S196_L001_R2_001.fastq
F3D142_S208_L001_R1_001.fastq  F3D150_S216_L001_R1_001.fastq  F3D9_S197_L001_R1_001.fastq
F3D142_S208_L001_R2_001.fastq  F3D150_S216_L001_R2_001.fastq  F3D9_S197_L001_R2_001.fastq
F3D143_S209_L001_R1_001.fastq  F3D1_S189_L001_R1_001.fastq    HMP_MOCK.v35.fasta
F3D143_S209_L001_R2_001.fastq  F3D1_S189_L001_R2_001.fastq    Mock_S280_L001_R1_001.fastq
F3D144_S210_L001_R1_001.fastq  F3D2_S190_L001_R1_001.fastq    Mock_S280_L001_R2_001.fastq
F3D144_S210_L001_R2_001.fastq  F3D2_S190_L001_R2_001.fastq    README.md
F3D145_S211_L001_R1_001.fastq  F3D3_S191_L001_R1_001.fastq    mouse.dpw.metadata
F3D145_S211_L001_R2_001.fastq  F3D3_S191_L001_R2_001.fastq    mouse.time.design
F3D146_S212_L001_R1_001.fastq  F3D5_S193_L001_R1_001.fastq    results
F3D146_S212_L001_R2_001.fastq  F3D5_S193_L001_R2_001.fastq    stability.batch
F3D147_S213_L001_R1_001.fastq  F3D6_S194_L001_R1_001.fastq    stability.files
F3D147_S213_L001_R2_001.fastq  F3D6_S194_L001_R2_001.fastq
```

We can run fastq on one file at a time or on all of them at once. Both of the following commands work. 


```
fastqc -o results/fastqc F3D0_S188_L001_R1_001.fastq 
fastqc -o results/fastqc *fastq
```

The standard output looks like this:

```
Started analysis of F3D0_S188_L001_R1_001.fastq
Approx 10% complete for F3D0_S188_L001_R1_001.fastq
Approx 25% complete for F3D0_S188_L001_R1_001.fastq
Approx 35% complete for F3D0_S188_L001_R1_001.fastq
Approx 50% complete for F3D0_S188_L001_R1_001.fastq
Approx 60% complete for F3D0_S188_L001_R1_001.fastq
Approx 75% complete for F3D0_S188_L001_R1_001.fastq
Approx 85% complete for F3D0_S188_L001_R1_001.fastq
```

Now, we can navigate to our results directory to view the results. 

```
cd results/fastqc
ls
```

For every input, there are two outputs: an html file and a ziped folder. The html files are of interest. 

```
F3D0_S188_L001_R1_001_fastqc.html    F3D150_S216_L001_R1_001_fastqc.html
F3D0_S188_L001_R1_001_fastqc.zip     F3D150_S216_L001_R1_001_fastqc.zip
F3D0_S188_L001_R2_001_fastqc.html    F3D150_S216_L001_R2_001_fastqc.html
F3D0_S188_L001_R2_001_fastqc.zip     F3D150_S216_L001_R2_001_fastqc.zip
F3D141_S207_L001_R1_001_fastqc.html  F3D1_S189_L001_R1_001_fastqc.html
F3D141_S207_L001_R1_001_fastqc.zip   F3D1_S189_L001_R1_001_fastqc.zip
F3D141_S207_L001_R2_001_fastqc.html  F3D1_S189_L001_R2_001_fastqc.html
F3D141_S207_L001_R2_001_fastqc.zip   F3D1_S189_L001_R2_001_fastqc.zip
F3D142_S208_L001_R1_001_fastqc.html  F3D2_S190_L001_R1_001_fastqc.html
F3D142_S208_L001_R1_001_fastqc.zip   F3D2_S190_L001_R1_001_fastqc.zip
F3D142_S208_L001_R2_001_fastqc.html  F3D2_S190_L001_R2_001_fastqc.html
```


> Click the raised hand :hand: if you have html files in a results directory.


Congratulations, you have successfully, launched and connected to an instance,
navigated the file system, downloaded data, installed programs, and executed programs at the command line. Here's an overview of some of the commands we used. 

:::info
#### Summary of commands

|Command |Description|
|-|-| 
|`pwd`| print name of current/working directory|
| `ls` [options] [path] | list directory contents | 
|`cd` [path]| change the working directory |
|curl -O [URL] | download a file from a URL and save it using the original file name | 
| `sudo apt install [program]` | install a program |
| `unzip [filename]` | uncompress filename
|rm [path] | removes (deletes) a file |
|mkdir -p [path/to/files] | creates a hierarchy of directories |
:::


:heavy_check_mark:  Let's take a 3 min break before moving on to the next section.  


## 5. Copying data from AWS instance onto your local computer

After processing your data in the cloud, you most likely need to copy some of your files to your local computer for viewing and sharing. In this section, we will use our ssh keys and public Domain Name System (DNS) to securely copy files from the cloud to our local computer using either a secure shell (`ssh`) or secure copy (`scp`).


### Windows Users

If you have a Windows machine, you will need to download a Terminal program. We recommend MobaXterm which is both a Terminal and an SSH client.

Read the following steps and/or watch [this short video tutorial](https://us06web.zoom.us/rec/play/1GfdPKpeJ5CVd8L6aPdYnOYuU3VmRmeoIHmdChyTNtUvpPIbezzdxdAGghAmsDPzhrGdi2SgkGSa_RqZ.7Ml7d475z1S9cItV).

**MobaXterm installation**

1. Go to the MobaXterm website to [download](https://mobaxterm.mobatek.net/)
2. Click on "GET MOBAXTERM NOW!"
3. The Home Edition works great and is free. Click "Download now".
4. Click on "MobaXterm Home Edition v20.6 (Portable edition)" and save as in your Downloads folder.
5. Go to your Downloads folder, click on the zipped folder, click "Extract all", click "Extract"
6. The MobaXterm application is now in the unzipped folder
7. Click on the MobaXterm application to open it!


Now that you have MobaXTerm installed you need to find the name and the address of your instance. To do so, let's reconnect to our instances. 

#### (Re)Connect to your EC2 instance

1. In a new browser tab or winder, navigate to the [instances page](https://us-west-1.console.aws.amazon.com/ec2/v2/home?region=us-west-1#Instances:).
2. Check the empty box next to your instance. 
3. Click the "Connect" button. 
4. Click the **SSH client** tab.
5. Find the "Example:" ssh command. Copy the last piece of information, which contains the public DNS for your instance and the computer name. It will look something like "ec2-54-193-121-227.us-west-1.compute.amazonaws.com"

![](https://i.imgur.com/EilADhq.png)

6. In MobaXterm, click on "Session"
7. Click on "SSH"
8. Enter the Public DNS as the "Remote host"
95. Check the box next to "Specify username" and enter "ubuntu" as the username
6. Click the "Advanced SSH settings" tab
7. Check box by "Use private key"
8. Use the document icon to navigate to where you saved the private key (e.g., "amazon.pem") from AWS on your computer. It is likely on your Desktop or Downloads folder
9. Click "OK"
10. A terminal session should open up with a left-side panel showing the file system of our AWS instance! 
11. Click on one of the FastQC html files to view it in a browser.

> Click the raised hand :hand: in zoom once you have viewed opened an html file.


### MacOS

Mac users do not need to install any additional programs to transfer files. You do however need to locate the ssh key file you saved at the beginning of the workshop. 

1. Open a Terminal window
2. Navigate your private key file and change the permissions using `chmod 400` to ensure your key is not publicly viewable. _Note: your .pem file may be in a different directory and have a different name. Modify the following commands accordingly._

```
cd ~/Desktop/
chmod 400 aws-jan-2022.pem
```


3. In a new browser tab or window, navigate to the [instances page](https://us-west-1.console.aws.amazon.com/ec2/v2/home?region=us-west-1#Instances:)
4. Check the empty box next to your instance. 
5. Click the "Connect" button. 
6. Click the **SSH client** tab.
7. Find the "Example:" ssh command. Copy the last piece of information, which contains the public DNS for your instance and the computer name. It will look something like "ubuntu@ec2-54-193-121-227.us-west-1.compute.amazonaws.com"

![](https://i.imgur.com/EilADhq.png)


8. Use the `scp` command on your local terminal to copy all the `.html` files. The `-i` option is used to specify the ssh key file. As with the copy (`cp`) command, you must specify both the location of the source file and the location of the copied file.  When specifying the source file, you must first include the Public DNS link (ec2-.....amazon.com) and the name of the user (@ubuntu).  To specify the path to the file, add a `:` after the DNS and the paste the path to the file.

Your command will look something like this. Remember to use your .pem file and your DNS. You can specify the current directory on your local computer with `.`

```
scp -i keys.pem  ubuntu@ec2-54-193-121-227.us-west-1.compute.amazonaws.com:~/MiSeq/results/fastqc/F3D141_S207_L001_R1_001_fastqc.html .
```

If this is your first time connecting to an instance, you may be prompted with the following question" Are you sure you want to continue connecting (yes/no/[fingerprint])?". Type "yes".

If you want to copy all the html files, you will need to put the path the files in single quotes to escape the wildcard. 

```
scp -i keys.pem 'ubuntu@ec2-54-193-121-227.us-west-1.compute.amazonaws.com:~/MiSeq/results/fastqc/*fastqc.html' .
```

> Click the raised hand :hand: in zoom once you have viewed opened an html file.


Congratulations! You have now successfully downloaded files from the cloud to your local computer. 

:::info
#### Summary of commands

|Command |Description|
|-|-| 
|`ssh -i keys.pem user@publicDNS`| connect to secure shell with ssh keys  |
| `scp -i keys.pem user@publicDNS:~/path/to/directory/file .` | download files with ssh keys | 
| `scp -i keys.pem file user@publicDNS:~/path/to/directory` | upload files with ssh keys | 
:::



## 6. Shutting down instances

The AWS Free Tier of services only remains free if you stay within the usage limits. If your instance is running in the cloud, you may be charged even if you aren't using it for computer power or storage. It is therefore good practice to shut down your instances when not in use. 

There are three options for shutting down instances. 

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

These accounts will remain available for 24 hours before your instructor deletes them. If you wish to return to your instance within the next 24 hours, stopping it is a good idea. If you are done practicing, terminating the instance is the best idea. 

To shut down an instance:

1. Navigate to the [instances page](https://us-west-1.console.aws.amazon.com/ec2/v2/home?region=us-west-1#Instances:)
2. Check the empty box next to your instance
3. Click the "Instance state" button
4. Select "Stop instance" or "Terminate instance" as appropriate




#### Exercise

Launch a t2.nano, Ubuntu 20.04 LTS - Focal instance in the **East US (Ohio) region**. Change the root storage volume to 16 GiB and add an additional EBS volume (8 GiB). 

> Hint
> - Go to Amazon Marketplace and search for the "Ubuntu 20.04 LTS - Focal". Should be the first result.
> - Look in tab 4 called "Add Storage" to add additional storage volumes.



## 7. Summary

In today's workshop, we covered the following topics:

- [x] AWS terminology and login
- [x] How to launch an instance 
- [x] How to connect to the instance
- [x] How to install and run a software programs on the instance 
- [x] How to terminate your instance 

We hope this workshop was helpful. Please complete the [post-workshop survey](https://docs.google.com/forms/d/e/1FAIpQLSe_IgvBrX_lr3B8Z_wWMTTt_qTfFwCR1lZLvEVe-BkCbOsKGw/viewform?usp=sf_link) to let us know if the workshop was useful or if you have any suggestions for improvement. 

We will send a follow email with links to the resources used today. 

If you have any questions, comments, or concerns, feel free to contact us at training@cfde.atlassian.net 

Check our [Events page](https://www.nih-cfde.org/events/) for information on upcoming workshops!

