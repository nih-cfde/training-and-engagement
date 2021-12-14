# An Introduction to UNIX for Cloud Computing - December 8, 2021

**When:** Wednesday, December 8th from 10:00 am-12:00 pm PDT

**Instructors:** Dr. Rayna Harris

**Helpers:**  Dr. Saranya Canchi and Jeremy Walter

**Zoom:** https://zoom.us/j/7575820324?pwd=d2UyMEhYZGNiV3kyUFpUL1EwQmthQT09

**Slides:** https://osf.io/dw3he/

**Computing environment:** [![Binder](https://binder.pangeo.io/badge_logo.svg)](https://mybinder.org/v2/gh/nih-cfde/training-rstudio-binder/data?urlpath=rstudio)

**Cheatsheet:** https://training.nih-cfde.org/en/latest/General-Tools/Cheat-Sheets/bash_cheatsheet/


#### Workshop Overview  

This free 2-hour workshop introduces the UNIX command line. It was designed for scientists and clinicians who need to use cloud-based and remote computers for basic and biomedical research. 

:::success
#### Today's Learning Goals

* Understand basic UNIX command structure 
* Navigate through hierarchical directory structures
* Read, write, create, copy, move, and remove files
* Understand wildcard expression 
* Combine commands into workflows
:::



The lesson materials were adapted from the UC Davis Data Lab's [Intro to Cloud Computing](https://ngs-docs.github.io/2021-august-remote-computing/) workshop, Data Carpentry's [Introduction to the Command Line for Genomics](https://datacarpentry.org/shell-genomics/) lesson, and the Lab for Data Intensive Biology's [Advanced Beginner/Intermediate Shell](https://dib-training.readthedocs.io/en/pub/2016-01-13-adv-beg-shell.html) workshop. Today's workshop is designed to prepare you for taking subsequent cloud computing and bioinformatic workshops. 


:::info
While we wait to get started --

1. :heavy_check_mark: Have you looked at the [pre-workshop resources page](https://github.com/nih-cfde/training-and-engagement/wiki/Resources-for-Workshop-Attendees)?

2. :pencil: Please fill out our pre-workshop survey if you have not already done so:
[Click here!](https://forms.gle/kf87XXWa7cZFDXUA9)


3. :hand: Raise your hand on Zoom when you've completed the pre-workshop survey.
:::


#### Table of Contents
[TOC]


---

## 1. The shell and terminal

The **shell** is a computer program that uses a command line interface (CLI) to give commands made by your keyboard to your operating system. Most people are used to interacting with a graphic user interface (GUI), where you can use a combination of your mouse and keyboard to carry out commands on your computer. 

We can use the shell through a **terminal** program. From the terminal. We can open programs, run analyses, create documents, delete files and create folders. 

For this remote workshop, we will be using a custom created create custom computing environment using [Binder](https://mybinder.org/). Click the **launch binder** button below, wait for it to launch, then open a new terminal window by clicking **Terminal**. 

[![Binder](https://binder.pangeo.io/badge_logo.svg)](https://mybinder.org/v2/gh/nih-cfde/training-rstudio-binder/data?urlpath=rstudio)

After the server launches and you select terminal, you should see an RStudio environment and the UNIX **prompt** `$`.

![](https://i.imgur.com/c3sQ3K4.png)

React with a :heavy_check_mark: when you have a terminal window open.
React with a :negative_squared_cross_mark: if you have issues opening a terminal.


:::danger
#### Getting help

When you open up the terminal in Binder, you may see a line of text or a **prompt statement** that tells us useful things such as the name of the directory we are currently in, our username, or what computer we are currently running terminal on. Let your instructors know if you see any warning or error messages. 

:::

:::success
#### Key points about the UNIX shell and terminal

* A shell is a program that reads commands and runs programs.
* We are using a remote terminal provided by myBinder.org  [![Binder](https://binder.pangeo.io/badge_logo.svg)](https://mybinder.org/v2/gh/nih-cfde/training-rstudio-binder/data?urlpath=rstudio)
:::


 
---

## 2. Navigating with `pwd`, `ls`, and `cd`

UNIX commands are like sentences that can be very simple or complex. The simplest commands consist of only the command name. Many require the name of a file or directory and allow specially formatted arguments, known as flags or options, which modify the default behavior of the program. The grammar of a shell allows you to combine existing tools into powerful pipelines and handle large volumes of data automatically. Sequences of commands can be written into a script, improving the reproducibility of workflows. The ease of getting things done via the shell will increase with your exposure to the program.


We should note that _folders_ are called **directories** at the command line. For all intents and purposes, they can be used interchangeably, but if you'd like more information please read about ["the folder metaphor"](https://en.wikipedia.org/wiki/Directory_%28computing%29#Folder_metaphor).


This Binder comes preloaded with data provided by your instructors.  _If you want to do these exercises locally, on your own computer, you can [download the data here](https://s3-us-west-1.amazonaws.com/dib-training.ucdavis.edu/shell-data.zip)._

For today's lesson, we will focus on three different sets of data. `books` contains ebooks such as A Tale of Two Cities and The Wizard of Oz.`southpark` contains a compressed csv file containing all the lines spoken by each character across 14 seasons. `seattle` contains data from the Open Seattle Data Portal, including a csv file with names of pets.  The `MiSeq` directory contains FASTQ and FASTA files that are commonly used in next-generation sequencing experiments. These data are useful for practicing commonly used UNIX commands to explore genome-scale data. 


The commands `pwd` and `ls` are two simple commands that can be used to answer the two commonly asked questions "where am I?" and "what files are here?". We will use these frequently through the next sections.

To answer the question "where am I?", we can use the **print working directory** or `pwd` command to see what directory we are currently located in. 


```
pwd
```


This will print **absolute path** to the directory where we are located. An absolute path shows the complete series of directories you need to locate either a directory or a file starting from the **root directory** of your computer. The absolute path to the root directory is `/`. A useful way to start thinking about directories and files is through levels. At the highest level of your computer, you have the root directory. Everything that is contained in your computer is located in directories below your root directory. 


The **home directory** is typically two levels down. For many mac users, the home directory is `/Users/USERNAME`. If you are using the Binder provided for this workshop, the home directory is `/home/jovyan`. Because the absolute path to the home directory is different for every user, you can refer to the home directory with the tilde symbol `~`.


```
/home/jovyan
```


Who or what is `Jovyan`? According to [Project Juypter](https://jupyter.readthedocs.io/en/latest/community/content-community.html#what-is-a-jovyan) the word “Jovian” describes several planets that share Jupiter-like properties. Much like the planet Jupiter and our solar system, the Jupyter community is large, distributed, and nebulous, so the word "Jovyan" is used to describe members of the community. Thus, the name of the User for this remote computer is "jovan". 



This **list** or `ls` command is a simple yet powerful command that is used to list the contents of your computer. It can be executed with or without optional flags and directories or files. Let's look at the contents in our working directory by using the `ls`.

```
ls
```

We can see the following files:

```
binder  books  images  MiSeq  README.md  seattle  southpark
```

If we want more information about the files, such as the date they were created and their file size, we can add "flags" `-l` for long listing format.

Flags (sometimes called options) allow us to finely control the behavior of the command. But how did we know to add `-l` after ls? The [`ls` manual ](https://man7.org/linux/man-pages/man1/ls.1.html) describes the command and all its options in the details. Like most commands, you can type the command followed `--help` to view the manual in your terminal.

```
ls --help
```


:::warning

#### Challenge: Multiple flags

You can use multiple flags or options at the same time to modify the behavior of a command. What does the command `ls` do when used with:
1. the `-h` option? 
1. the `-l` and the `-h` option?
1. the `-l`, the `-h`, and the `-F` option?

:::spoiler

1. The `-h` option makes the file size human readable, but it is only noticeable if you are printing the file size. 
1. If you use both the `-h` option and the `-l` option (with`ls -lh` or `ls -l -h`), this makes the file size ‘human readable’, _i.e._ displaying something like 5.3K instead of 5369.
1. The -F flag will class the file types by appending an identifier. This works best if there are directories present. 
:::

Now we have seen how to list around our computers and what is located in the directory we are. But some of the beauty of the shell is that we can execute activities in locations that we are not currently in. To do this we can either use an absolute path or a relative path. A **relative path** is the path to another directory from the one you are currently in. 


To move from one directory to the other, we use the `cd` command to **change directories**. 
Let's navigate into the `books` using the `cd` command followed by

```
cd books
```

Now, we can use the `pwd` and/or `ls` commands to confirm that we did indeed change directories.  

```
pwd
ls
```


Because you can change directories using either the relative or absolute path, there are multiple ways to successfully move up or down in the directory hierarchy.
Let's return to our home directory using the `cd` command and a relative path, then print the working directory to confirm.  


:::warning
#### Challenge: Navigating with relative and absolute paths

Starting from `books`, which of the following commands could Jovyan use to navigate to the `MiSeq` directory? 


1. `cd MiSeq`
2. `cd ./MiSeq`
3. `cd ~/MiSeq`
4. `cd ../MiSeq`
5. `cd /MiSeq`
6. `cd ../../Miseq`
7. `cd /home/jovyan/MiSeq`


:::spoiler

1. No, MiSeq does not exist in the current working directory.
2. No, MiSeq does not exist in the current working directory.
3. Yes, MiSeq is in the home directory.
4. Yes, MiSeq is in the directory one level above.
5. No, MiSeq is not in the root directory.
6. No, MiSeq is not in the directory two levels above.
7. Yes, this is the full path to MiSeq.
:::

 
Let's practice using the cd and ls commands to explore files in different directories.  
 
#### books

What files are in the `books` directory?

```
ls books
cd ~/books
ls 
```

We can see the following files:

```
A-tale-of-two-cities.txt	TheWonderfulWizardOfOz.txt.gz
Alice_in_wonderland.txt		book.txt
PeterPan.txt			    get_books.sh
README.md
```

#### southpark

How large are all the files in the `southparkdata` directory?

```
cd ~/southpark
ls -l 
```

We will see the following:

```
-rw-r--r-- 1 jovyan jovyan 1866609 Aug  4 04:34 All-seasons.csv.gz
-rw-r--r-- 1 jovyan jovyan     252 Aug  4 04:34 License.md
-rw-r--r-- 1 jovyan jovyan     345 Aug  4 04:34 README.md
```

#### MiSeq


How large are all the `.fastq` files?


```
cd ~/MiSeq
ls *fastq
```




:::success
#### Key Points

|Command |Description|
|-|-| 
|`pwd`| print name of current/working directory|
| `ls` [options] | list directory contents | 
|`cd` [path]| change the working directory |

|Path |Description|
|-|-| 
|`/`| root directory|
| `~/` | home directory | 
|`./` | current or working directory |
|`../` | directory one level up |


:::


## 3. Reading with `head`, `cat`, `less` and `more`

Now that we know what files exist on our computer, it's time to look at the contents of the file. There are multiple ways to look at the contents of a file. 

The `cat` command prints the entirety of a file to the stdout of our computer. We can scroll through files using the `less` command. Less is a safe way of looking at the contents of a file without the ability to change it. `head` prints, by default the first 10 lines of a file.


All three of the commands use the same syntax:

```
head [filename]
cat [filename]
less [filename]
```

:::info

#### Tab completion
You can use TAB to do filename completion, so if you type `cat R` and then press your Tab key once, it will autocomplete if there is a unique match. If there is more than one match, the first Tab will do nothing, and the second will show all the possible matches.
:::


Let's navigate to the `books` directory and use the `head` command to view the `README.md` file. 

```
cd ~/books/
head README.md
```

You should see an output that looks like this. The `README.md` and `License.md` files are written in [Markdown](https://en.wikipedia.org/wiki/Markdown). To learn more about Markdown syntax, read this excellent [Markdown guide](https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet).

```
# Tale of Two Cities

downloaded from [Project Gutenberg](https://www.gutenberg.org/ebooks/98)

```


Now we can view the file with `head`, `cat`, or `less` and `tail`. 

```
head book.txt
tail book.txt
cat book.txt
less book.txt
```


We can there are a lot more books, and we can look at the first few lines of all the txt files with the *. 

```
head *.txt
```

Notice, there is one book that is compressed. We can uncompress it with the command `gunzip`.

```
gunzip TheWonderfulWizardofOz.txt.gz
```


::: success
#### Key UNIX commands for viewing files
| Command [OPTION] | Description |
| -------- | -------- | 
|`head [filename]` | print first 10 lines of  `FILENAME` | 
|`cat [filename]`| print `FILENAME`'s contents to stdout|
|`less [filename]`|view `FILENAME` without printing  to stdout |
| gunzip -k [filename] | uncompress a file and keep the original
:::


## 4. Working with files and directories using `cp`, gunzip, `mv`, `mkdir` 

We are quite used to copying and moving files using a GUI. These functions can be also carried out at the command line.

The `cp` and `mv` commands can be used to copy and move (or rename) files and directories respectively. For both commands, you must specify the old and new names. Specifying the path is necessary if you want to move files out of the current working directory.


```
cp [original-filename] [copy-filename]
mv [original-filename] [new-filename]
```

Let's make a copy of some raw data before we start modifying it. 

```
cp book.txt book-copy.txt
ls
```


The `mv` command can be used to either move files to a new location or to rename them (which is essentially moving the contents from the old filename to the new file name. Let's use the `mv` command to rename the copied and compressed file back to the original name.

```
mv book-copy.txt book-2cities.txt

```

Take care when renaming files. It is good practice to keep track of changes in file names and links to the source data. The commands used to get these books are stored in `get_books.sh`.  We will talk about how to execute this script later.

```
head get_books.sh
```


:::warning

### Challenge: Creating and deleting hierarchy of directories

Now you know how to copy and move files, but you may encounter errors if you try to move files to a directory that doesn't exist. But, have no fear, we can create new directories at the command line with the command `mkdir` followed by the path to the directories you want to create. 

What happens when you run the following commands?

```
1. mkdir data results images/
2. mkdir -p data/results/images
```



:::spoiler Answer

1. data, results, and images all created in the working directory. 2. images is a subdirectory of results, which is a subdirectory of data 

:::

::: warning
Now, imagine you could create the perfect directory hierarchy for a project. What would it look like? Type a command or series of commands to create your ideal directory structure. Share with your group. 

:::spoiler Hint

An example project directory could be setup like this:

```
mkdir awesome/
cd awesome/
mkdir -p data/ results/2020/ results/2021 images/ notes/
ls *
```
:::

::: warning
If you created some files or directories that you don't want, you can remove them with the `rm` and `rmdir` commands. How could you remove `data/results/`

:::spoiler A solutions


```
rmdir data/results/images
rmdir data/results
```

or 

```
rm -r data/results
```
:::



:::success

#### Key Commands: Working with files and directories
| Command | Description |
| -------- | -------- | 
|cp [old] [new] | copies a file | 
|mv [old] [new]| moves or renames a file or directory |
|rm [path] | removes (deletes) a file |
|mkdir [path] | creates a new directory |
|rmdir [path] | removees an empty  directory |
:::



## 5. Finding things with `grep`, `find`, `cut`




A big part of data science is making sure what you _expect_ in a particular file is what you _have_ in that file. This is fairly easy when your files are small but is challenging when the files are much larger than your screen.


To explore this topic in more detail, navigate to the `data/MiSeq/` directory.


```
cd ~/MiSeq
ls
```

This directory contains multiple [**FASTQ** files]((https://en.wikipedia.org/wiki/FASTQ_format)). A FASTQ file normally uses four lines per sequence.

* Line 1 begins with a '@' character and is followed by a sequence identifier and an optional description (like a FASTA title line).
* Line 2 is the raw sequence letters.
* Line 3 begins with a '+' character and is optionally followed by the same sequence identifier (and any description) again.
* Line 4 encodes the quality values for the sequence in Line 2, and must contain the same number of symbols as letters in the sequence.

A FASTQ file containing a single sequence might look like this:


::: info
An example FASTQ file
```
@SEQ_ID
GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65
```
:::



We can use the `cat` command to print fastq files to the screen, but thousands of lines of text would crowd your screen. Instead, let's use the `head` command to view the first 8 lines file. You can copy the file name and paste it into the console or you can type and use tab complete to pick a particular file. 


```
head -n 4 F3D0_S188_L001_R1_001.fastq
```

```
@M00967:43:000000000-A3JHG:1:1101:18327:1699 1:N:0:188
NACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGCCTGCCAAGTCAGCGGTAAAATTGCGGGGCTCAACCCCGTACAGCCGTTGAAACTGCCGGGCTCGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACCCCGATTGCGAAGGCAGCATACCGGCGCCCTACTGACGCTGAGGCACGAAAGTGCGGGGATCAAACAG
+
AABABBFFFGGGGGGGGGGGGGGGGHHHHHHHGGGHHHHHGHGGGGGGGHGGGGGGHHHHHHHHHHGGGGGHHHHGHGGGGGGHHBGHGDGGGGGHHHGGGGHHHHHHHHGGGGGHG@DHHGHEGGGGGGBFGGEGGGGGGGG.DFEFFFFFFFDCFFFFFFFFFFFFFFFFFFFFFFFFFFDFDFFFEFFCFF?FDFFFFFFFFAFFFFFFFFFFFBDDFFFFFEFADFFFFFBAFFFA?EFFFBFF
@M00967:43:000000000-A3JHG:1:1101:14069:1827 1:N:0:188
TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGCCTGCCAAGTCAGCGGTAAAATTGCGGGGCTCAACCCCGTACAGCCGTTGAAACTGCCGGGCTCGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACCCCGATTGCGAAGGCAGCATACCGGCGCCCTACTGACGCTGAGGCACGAAAGTGCGGGGATCAAACAG
+
3AA?ABBDBFFBEGGEGGGGAFFGGGGGHHHCGGGGGGHFGHGGCFDEFGGGHGGGEGF1GGFGHHHHHGGEGGHHHHHFGGGGGGHHHHHGGGGCDDGHHGGGFHHHHHHHHCD@CCHGGGGHEHGGG@GFGGGGGGG@BGGGEGCEBFFFBFFB;9@EFFFEFFFFFFFFFFFFAFBBBFFFFFBBBFFFFBBBFFFFFFFFFFFBBBBBBBFFFFFFFFFDDFAFFFFF.AF9/FBBBBB.EAFFE?F
```

`head` prints the first ten lines of a file out onto your screen. Similarly, the `tail` command prints the last 10 lines of a file. 

```
tail -4 F3D0_S188_L001_R1_001.fastq
```

```
@M00967:43:000000000-A3JHG:1:1105:19125:28016 1:N:0:188
TACGTAGGGGGCAAGCGTTATCCGGAATTACTGGGTGTAAAGGGAGCGTAGACGGTAATGCAAGTCTGGAGTGAAAGGCGGGGGCCCAACCCCCGGACTGCTCTGGAAACTGTGTAACTGGAGTGCAGGAGAGGCAGGCGGAATTCCTAGTGTAGCGGTGAAATGCGTAGATATTAGGAGGAACACCAGTGGCGAAGGCGGCCTGCTGGACTGTAACTGACGTTGAGGCTCGAAAGCGGGGGGGGCAAAAA
+
>AAAAFFBBBBDGFGGGEGFGGAEEEGGHHGHHHHFFEHHHHGGHGGGGFGGHGGGEGHHGHHHHHHHHGHGHFFHHHHGGGGG@DGDGDGGGGGGDDGHHHHGHHCFHGHHHH0GGFFHHHGFFEFGGGGGGGGGG@.BB-@EF/BFFF9//9FF;B.;.:F//9A.;9.9;B/////BFDB9FF.AAB/:?B-;-99@B-9;EA/;B99BE99B///B:9.;:9F..;9=-FD.9A-;@BB----?F..
@M00967:43:000000000-A3JHG:1:1105:20429:28046 1:N:0:188
TACGGAGGATTCAAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGTTCGATAAGTTAGAGGTGAAATCCCGGGGCTCAACTCCGGCACTGCCTCTGATACTGTCGGGCTAGAGTTTAGTTGCGGTAGGCGGAATGTATGGTGTAGCGGTGAAATGCATAGAGATCATACAGAACACCGATTGCGAAGGCAGCTTACCAAACTACGACTGACGTTGAGGCACGAAAGCGTGGGGAGCAAACA
+
BBBBBBBBBFFFGGGFFGGGGGGGGGGGHHHHHHHGGGHGHHHGHGGGGGGGHGGEGGFGGHGHHHHHHHHHGHHHHEHHGGGGGGHHHHHHHGGGGGHHHGFHHHHGGHHHHHGGGGGHFGH?GHHGHHHHGG<CDFHGGGGGHCFGHHHHEHFHCGCGEFFGGGGGGEGGEFGGGFFFFGGG/FDCFFFFFADFFFFFDFFFFFFFFFABFFFFF?DFFFE.EEFFFEFFFFAAFADFFFFECDEFB..
```

**FASTQ** files should not be confused with **FASTA** files. FASTQ files contain information about the quality of the sequence, but FASTA files only contain the sequence and an identifier.

::: info
An Example FASTA file
```
> SEQUENCE_1
MTEITAAMVKELRESTGAGMMDCKNALSETNGDFDKAVQLLREKGLGKAAKKADRLAAEG
LVSVKVSDDFTIAAMRPSYLSYEDLDMTFVENEYKALVAELEKENEERRRLKDPNKPEHK
IPQFASRKQLSDAILKEAEEKIKEELKAQGKPEKIWDNIIPGKMNSFIADNSQLDSKLTL
MGQFYVMDDKKTVEQVIAEKEKEFGGKIKIVEFICFEVGEGLEKKTEDFAAEVAAQL
>SEQUENCE_2
SATVSEINSETDFVAKNDQFIALTKDTTAHIQSNSLQSVEELHSSTINGVKFEEYLKSQI
ATIGENLVVRRFATLKAGANGVVNGYIHTNGRVGVVIAAACDSAEVASKSRDLLRQICMH
```
:::

Let's look at a synthetic FASTA file. 

```
head -4 HMP_MOCK.v35.fasta 
```


```
> A.baumannii.1 
TGGGGAATATTGGACAATGGGGGGAACCCTGATCCAGCCATGCCGCGTGTGTGAAGAAGGCCTTATGGTTGTAAAGCACTTTAAGCGAGGAGGAGGCTACTTTAGTTAATACCTAGAGATAGTGGACGTTACTCGCAGAATAAGCACCGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGGTGCGAGCGTTAATCGGATTTACTGGGCGTAAAGCGTGCGTAGGCGGCTTATTAAGTCGGATGTGAAATCCCCGAGCTTAACTTGGGAATTGCATTCGATACTGGTGAGCTAGAGTATGGGAGAGGATGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGATGGCGAAGGCAGCCATCTGGCCTAATACTGACGCTGAGGTACGAAAGCATGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCATGCCGTAAACGATGTCTACTAGCCGTTGGGGCCTTTGAGGCTTTAGTGGCGCAGCTAACGCGATAAGTAGACCGCCTGGGGAGTACGGTC
> A.odontolyticus.1
TGGGGAATATTGCACAATGGGCGAAAGCCTGATGCAGCGACGCCGCGTGAGGGATGGAGGCCTTCGGGTTGTAAACCTCTTTCGCTCATGGTCAAGCCGCAACTCAAGGTTGTGGTGAGGGTAGTGGGTAAAGAAGCGCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGCGCGAGCGTTGTCCGGAATTATTGGGCGTAAAGGGCTTGTAGGCGGTTGGTCGCGTCTGCCGTGAAATCCTCTGGCTTAACTGGGGGCGTGCGGTGGGTACGGGCTGACTTGAGTGCGGTAGGGGAGACTGGAACTCCTGGTGTAGCGGTGGAATGCGCAGATATCAGGAAGAACACCGGTGGCGAAGGCGGGTCTCTGGGCCGTTACTGACGCTGAGGAGCGAAAGCGTGGGGAGCGAACAGGATTAGATACCCTGGTAGTCCACGCTGTAAACGTTGGGCACTAGGTGTGGGGGCCACCCGTGGTTTCTGCGCCGTAGCTAACGCTTTAAGTGCCCCGCCTGGGGAGTACGGCC
```

FASTQ and FASTA files are often used in combination to map reads to a genome or transcription. You've seen three ways to read large files. Next, we will learn how to use and manipulate the files.  


Sometimes you know a file or directory exists, but you can't find it. Sometimes you want to find many files with similar properties. This is where the wildcard (`*`) comes in handy.  What do the following commands do?

 
1. `ls *` 
1. `ls MiSeq/F3D*`
1. `ls MiSeq/*fasta`
 

:::spoiler

1. `ls *` lists files in the working directory and 1 level down. 
1. `ls MiSeq/F3D*` lists files in the data/MiSeq directory that start with "F3D".
1. `ls MiSeq/*fasta` lists files in the data/MiSeq directory that end with "fasta".

:::

A lot of the time we want to know if a file contains what we expect. A useful thing to do is to be able to **search the contents of files** for a particular string of characters you would like to find.  We can use the file pattern searcher `grep` to find things.



The `MiSeq/` directory contains many of the sequence files ending in`.fastq`. We expect these files to contain information in a particular format throughout the file with four lines of information for each sequence string. Looking through a million line file using less will take a long time. Rather than manually looking at the whole file, we can print only a portion of the file's contents to standard output. 

Let's say you'd like to find the sequence `CATTAG` in your MiSeq files. We can also use the **wildcard** regular expression to search `CATTAG` in all of the fastq files located in our current working directory:

```
cd ../data/MiSeq/
grep CATTAG F3D0_S188_L001_R2_001.fastq
grep CATTAG *.fastq
```


:::warning
#### CHALLENGE: grep
What line does `CATTAG` occur on in `F3D141_S207_L001_R1_001.fastq`? 

:::spoiler Hint
Use `grep --help` to search for `grep` options related to line number.
`grep -n [filename]`` will print the line number.
:::

How many lines are in the file total? The word count (`wc`) command will give us that answer.

```
wc F3D0_S188_L001_R2_001.fastq
```

By default, `wc` prints the characters, words, and lines in a file. To extract just the line numbers, we use `wc -l`. We will use this more later. 


```
wc -l F3D0_S188_L001_R2_001.fastq
```

:::success
#### Key Points:
|Command|Description|
|-|-
|`find [filename]` | finds files with specific properties that match patterns|
| `grep [option] [filename]`  | selects lines in files that match patterns|
| `wc` [filename] | will print the total characters, words, and lines in a file.
:::





## 6. Redirection, appending, and piping with `>`, `>>`, and `|`

By default, many UNIX commands like `cat` send output to something called
standard out, or "stdout". This is a catch-all phrase for "the basic
place we send regular output." (There's also standard error, or "stderr",
which is where errors are printed; and standard input, or "stdin", which
is where input comes from.)

Much of the power of the UNIX command line comes from working with
stdout output, and if you work with UNIX a lot, you'll see characters
like the `>` (redirect), `>>` (append) , and `|` (pipe) thrown around. These
are redirection commands that say, respectively, "send stdout to a new
file", "append stdout to an existing file", and "send stdout from one
program to another program's stdin."

If you know you want to save an output file, you can use the redirect symbol `>`. Note, if you want to save a file in a different directory, that directory must exist.

```
mkdir results
grep CATTAG *.fastq > results/files-with-CATTAG.txt
```


Let's pipe (`|`) the output of a grep search to the command head to only show the first 10 lines. This is very handy when you only want to view snapshot of the data. 



```
grep CATTAG *.fastq 
grep CATTAG *.fastq | wc -l
```

The pipe is a useful way to combine multiple commands on line line into a pipe. 

However, the two lines above told us, how many times CATTAG appeared in all of the .fastq files in total. If we want to know who many times it occurs in each, we need a for loop. 

A for loop looks like this

```
for [thing] in [list of things]
do
command $[thing]
done
```


We could use this with `grep` on our fastq files like so. 


grep CATTAG *.fastq 

```
for file in *fastq
do
grep CATTAG $file | wc -l
done
```

If we cant to know which file these results came from, we can add an echo statement to echo the variable.


```
for file in *fastq
do
echo $file
grep CATTAG $file | wc -l
done
```



When working with `csv` files, sometimes you only want to look for patterns in 1 column versus all the columns. So, you need to filter the data first to only have the column of interest.  

Notice, if you search for characters' names in the South Park TV series, grep will return both instances where the character spoke a line and where the line mentions their name. 

```
cd ~/SouthParkData
gunzip -k All-seasons.csv
grep Kenny All-seasons.csv
```

Let's say you want to know which of your Southpark charcters had the most spoken lines. You could do like by counting how often their name appeared in the 3rd column of the csv file. We need the command `cut`.  We use the `-f3` to specify the third column and `-d,` to specify that it is a csv. 

```
head All-seasons.csv
cut -d, -f3 All-seasons.csv
```

Now, after we extract only the column with character names, we can search for our favorite characters and count how many spoken lines they had. 

```

cut -d, -f3 All-seasons.csv | grep Kenny | wc -l
```

If we want to do this on more than 1 character, we can use a for loop. 


```
for character in Kenny Cartman Chef Kyle
do
echo $character
cut -d, -f3 All-seasons.csv | grep $character | wc -l
done
```

:::success
#### Key points: redirecting outputs

|Command|Description|
|-|-
|`\|` | pipes the standard output to a new command|
| `>`  | redirects the standard output to a new file |
| `>>`  | append the standard output to a new or existing file|

:::



## Concluding thoughts

This lesson focused on file and directory exploration because that's
something everyone needs to know, and all these commands will work on
pretty much any computer that is running a UNIX compatible shell (including
Mac OS X and Windows Subsystem for Linux). 

We've shown you a whole plethora of hopefully
not-too-confusing options for editing and working with text files.


The redirection and compression stuff is really useful, but again,
you just need to know it exists and that there's this tutorial on it.


These skills may seem confusing, but they will become second nature if
you use them regularly.

If you want to save all the commands we used today, you can use the `history` command to print out all the commands you typed.

```
history
```

You can save the file in your home directory with: 

```
history > ~/history-2021-nov-17.txt
```


The binder and this documentation page will stay working for the foreseeable
future, so please feel free to come back and revisit some of these commands!

Google (and especially stackoverflow) is your friend! Use Internet
search whenever you have questions about what a command does, or what
commands to use to achieve a particular task. The wedbsite [Explain Shell](https://explainshell.com/) is great for defining what each command and flag does.




::: success
#### Key points 
This workshop teaches a dozen commonly used UNIX commands that can be combined to perform power, reproducible bioinformatic workflows. The commands taught `pwd` `ls`  `cd` `cat` `head` `less` `cp` `mv` `rm` `mkdir` `grep` `wc` `cut` `gunzip` and `gzip` (and probably a few others).
:::



### Today's Commands


```
clear
pwd
PS1="\w $ "

ls
ls --help
ls -l books
ls -l MiSeq
ls -lh MiSeq
ls -lh *

cd books
pwd
ls
cd ../Miseq
pwd
ls
ls *fastq
ls F3*

cd ~books

head README.md

head book.txt
tail book.txt
cat book.txt
less book.txt

head *txt

wc -l *txt
wc -l *txt | sort
wc -l *txt | sort -nr

gunzip TheWonderfulWizardofOz.txt
wc -l *txt | sort -nr

rm book.txt
mkdir results

wc -l *txt | sort -nr > results/book_lengths.txt
cat results/book_lengths.txt

grep "Chapter" *txt
grep -i "Chapter" *txt
grep "CHAPTER" *txt > results/chapter_titles.txt
cat results/chapter_titles.txt

grep "The" *txt
grep "^The" *txt
grep -w "^The" *txt
grep -w -A 1 "^The" *txt
grep -w  "^The" *txt > results/The.txt
cat results/The.txt

head README.md
head -4 F3D0_S188_L001_R1_001.fastq
tail -4 F3D0_S188_L001_R1_001.fastq
less F3D0_S188_L001_R1_001.fastq

wc README.md
wc F3D0_S188_L001_R1_001.fastq
wc -l F3D0_S188_L001_R1_001.fastq	
wc -l *fastq
wc -l *R1*fastq
wc -l *R1*fastq | sort -nr

head -4 F3D0_S188_L001_R1_001.fastq
grep "^@M" F3D0_S188_L001_R1_001.fastq
grep "^@M" F3D0_S188_L001_R1_001.fastq | wc -l
grep "^@M" *R1*.fastq | wc -l

for file in *R1*.fastq
do
echo $file
grep "^@M" $file | wc -l
done

mkdir results


for file in *R1*.fastq
do
echo $file >> results/read_count.txt
grep "^@M" $file | wc -l >> results/read_count.txt
done

head results/read_count.txt

for file in *R1*.fastq
do
echo $file >> results/samples.csv
grep "^@M" $file | wc -l >> results/count.csv
paste -d , results/samples.csv results/count.csv > results/read_count.csv
done

head results/read_count.csv

cd
history
history > notes_dec_7.txt

```


