---
title: Naviagate
---

### 

UNIX commands are like sentences that can be very simple or complex. The simplest commands consist of only the command name. Many require the name of a file or directory and allow specially formatted arguments, known as flags or options, which modify the default behavior of the program. The grammar of a shell allows you to combine existing tools into powerful pipelines and handle large volumes of data automatically. Sequences of commands can be written into a script, improving the reproducibility of workflows. The ease of getting things done via the shell will increase with your exposure to the program.

We should note that _folders_ are called **directories** at the command line. For all intents and purposes, they can be used interchangeably, but if you'd like more information please read about ["the folder metaphor"](https://en.wikipedia.org/wiki/Directory_%28computing%29#Folder_metaphor).

This Binder comes preloaded with data provided by your instructors.  _If you want to do these exercises locally, on your own computer, you can [download the data here](https://s3.us-west-1.amazonaws.com/dib-training.ucdavis.edu/shell-data2.zip)._

The commands `pwd` and `ls` are two simple commands that can be used to answer the two commonly asked questions "where am I?" and "what files are here?". We will use these frequently through the next sections.

### `pwd`

To answer the question "where am I?", we can use the **print working directory** or `pwd` command to see what directory we are currently located in. 

```bash
pwd
```

This will print the **absolute path** to the directory where we are located. An absolute path shows the complete series of directories you need to locate either a directory or a file starting from the **root directory** of your computer. The absolute path to the root directory is `/`. A useful way to start thinking about directories and files is through levels. At the highest level of your computer, you have the root directory. Everything that is contained in your computer is located in directories below your root directory. 

The **home directory** is typically two levels down. For many Mac users, the home directory is `/Users/USERNAME`. If you are using the Binder provided for this workshop, the home directory is `/home/jovyan`. Because the absolute path to the home directory is different for every user, you can refer to the home directory with the tilde symbol `~`.

```bash
/home/jovyan
```

### `ls`

The **list** or `ls` command is a simple yet powerful command that is used to list the contents of your computer. It can be executed with or without optional flags and directories or files. Let's look at the contents in our working directory by using the `ls`.

```bash
ls
```

We can see the following files:

```bash
books          images  README.md             seattle
CFDE-logo.png  MiSeq   rstudio-terminal.png  
```

If we want more information about the files, such as the date they were created and their file size, we can add "flags" `-l` for long listing format.

```bash
ls -l
```

```bash
drwxr-xr-x 2 jovyan jovyan   4096 Jan 18 21:13 books
-rw-r--r-- 1 jovyan jovyan  71154 Jan 18 21:13 CFDE-logo.png
drwxr-xr-x 2 jovyan jovyan   4096 Jan 18 21:13 images
drwxr-xr-x 2 jovyan jovyan   4096 Jan 18 21:13 MiSeq
-rw-r--r-- 1 jovyan jovyan   2089 Jan 18 21:13 README.md
-rw-r--r-- 1 jovyan jovyan 188705 Jan 18 21:13 rstudio-terminal.png
drwxr-xr-x 2 jovyan jovyan   4096 Jan 18 21:13 seattle
```

Flags (sometimes called options) allow us to finely control the behavior of the command. But how did we know to add `-l` after ls? The [`ls` manual ](https://man7.org/linux/man-pages/man1/ls.1.html) describes the command and all its options in detail. Like most commands, you can type the command followed `--help` to view the manual in your terminal.

```bash
ls --help
```

=== "Challenge"

You can use multiple flags, wildcards, and specify directories to modify the behavior of a command. What does the command `ls` do when used with the following option:

1. `ls -a`
2. `ls -F`
3. `ls -aF`

=== "Answer"

1. The `-a` flag will list hidden files and directories.  
2. The `-F` flag will class the file types by appending an identifier. This works best if there are directories present. 
3.  We can combine `-a` and `-F` to be `-aF` to use both options.

### `cd`

Now we have seen how to list the contents of folders on our computers and what is located in the directory we are presently in. But some of the beauty of the shell is that we can perform activities in locations that we are not currently in. To do this we can either use an absolute path or a relative path. A **relative path** is the path to another directory from the one you are currently in. An **absolute path** starts from the root and ends in the apropriate subdirectory. 

To move from one directory to the other, we use the `cd` command to **change directories**. We can use the `pwd` and/or `ls` commands to confirm that we did indeed change directories.  Because you can change directories using either the relative or absolute path, there are multiple ways to successfully move up or down in the directory hierarchy.

Let's return to our home directory using the `cd` command and a relative path, then print the working directory to confirm.  
 
Let's practice using the cd and ls commands to explore files in different directories.  

Because books/ is in our working directory, we can navigate there with a relative path. What files are in the `books` directory and how large are they?
bash
```
cd books/
pwd
ls -lh
```

We can see the following files:

```bash
-rw-r--r-- 1 jovyan jovyan 171K Jan 18 21:13 Alice_in_wonderland.txt
-rw-r--r-- 1 jovyan jovyan 789K Jan 18 21:13 A-tale-of-two-cities.txt
-rw-r--r-- 1 jovyan jovyan 789K Jan 18 21:13 book.txt
-rw-r--r-- 1 jovyan jovyan 282K Jan 18 21:13 PeterPan.txt
-rw-r--r-- 1 jovyan jovyan 1.1K Jan 18 21:13 README.md
-rw-r--r-- 1 jovyan jovyan  80K Jan 18 21:13 WizardOfOz.txt.gz
-rw-r--r-- 1 jovyan jovyan  12M Jan 18 21:13 yeast.fasta
```

=== "Challenge"

Starting from `books`, which of the following commands could Jovyan use to navigate to the `MiSeq` directory? 


1. `cd MiSeq`
2. `cd ./MiSeq`
3. `cd ~/MiSeq`
4. `cd /home/jovyan/MiSeq`
5. `cd ../MiSeq`
6. `cd ../../MiSeq`
7. `cd /MiSeq`

=== "Answer"

1. No, MiSeq does not exist in the current working directory.
2. No, MiSeq does not exist in the current working directory.
3. Yes, MiSeq is in the home directory.
4. Yes, this is the full path to MiSeq.
5. Yes, MiSeq is in the directory one level above.
6. No, MiSeq is not in the directory two levels above.
7. No, MiSeq is not in the root directory.
:::


Most, but not all of the files in the MiSeq directory are .fastq files. Which .fastq files are the largest? We can use the wildcard `*` to list only files that end in .fastq. We can use the `-S` option to sort by size.


```bash
cd ~/MiSeq
pwd
ls -lhS *.fastq
```

```bash
-rwxr-xr-x 1 jovyan jovyan  11M Jan 18 21:13 F3D2_S190_L001_R1_001.fastq
-rwxr-xr-x 1 jovyan jovyan  11M Jan 18 21:13 F3D2_S190_L001_R2_001.fastq
-rwxr-xr-x 1 jovyan jovyan 9.2M Jan 18 21:13 F3D147_S213_L001_R1_001.fastq
-rwxr-xr-x 1 jovyan jovyan 9.2M Jan 18 21:13 F3D147_S213_L001_R2_001.fastq
-rwxr-xr-x 1 jovyan jovyan 7.1M Jan 18 21:13 F3D149_S215_L001_R1_001.fastq
-rwxr-xr-x 1 jovyan jovyan 7.0M Jan 18 21:13 F3D149_S215_L001_R2_001.fastq
-rwxr-xr-x 1 jovyan jovyan 6.7M Jan 18 21:13 F3D148_S214_L001_R1_001.fastq
-rwxr-xr-x 1 jovyan jovyan 6.7M Jan 18 21:13 F3D148_S214_L001_R2_001.fastq
-rwxr-xr-x 1 jovyan jovyan 4.3M Jan 18 21:13 F3D6_S194_L001_R1_001.fastq
-rwxr-xr-x 1 jovyan jovyan 4.3M Jan 18 21:13 F3D6_S194_L001_R2_001.fastq
```


### Key points

|Command |Description|
|-|-| 
|`pwd`| print name of current/working directory|
| `ls` [options] [path] | list directory contents | 
|`cd` [path]| change the working directory |

|Path |Description|
|-|-| 
|`/`| root directory|
| `~/` | home directory | 
|`./` | current or working directory |
|`../` | directory one level up |
