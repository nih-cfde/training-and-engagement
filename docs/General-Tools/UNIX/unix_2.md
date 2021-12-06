# Navigating directories

UNIX commands are like sentences that can be very simple or complex. The simplest commands consist of only the command name. Many require the name of a file or directory and allow specially formatted arguments, known as flags or options, which modify the default behavior of the program. The grammar of a shell allows you to combine existing tools into powerful pipelines and handle large volumes of data automatically. Sequences of commands can be written into a script, improving the reproducibility of workflows. The ease of getting things done via the shell will increase with your exposure to the program.


We should note that _folders_ are called **directories** at the command line. For all intents and purposes, they can be used interchangeably, but if you'd like more information please read about ["the folder metaphor"](https://en.wikipedia.org/wiki/Directory_%28computing%29#Folder_metaphor).

![](https://i.imgur.com/tS4uw77.png)


This Binder comes preloaded with data provided by your instructors.  _If you want to do these exercises locally, on your own computer, you can [download the data here](https://s3-us-west-1.amazonaws.com/dib-training.ucdavis.edu/shell-data.zip)._

For today's lesson, we will focus on three different sets of data. `2cities` contains a compressed text file containing the book A Tale of Two Cities, `SouthParkData` contains a compressed csv file containing all the lines spoken by each character across 14 seasons. This dataset is useful for teaching UNIX commans on medium sized text data. The `data/MiSeq` directory contains FASTQ and FASTA files that are commonly used in next-generation sequencing experiments. These data are useful for teaching commonly used UNIX commands for exploring genome-scale data. 



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
MiSeq		binder		images		southpark
README.md	books		seattle
```

If we want more information about the files, such as the date they were created and their file size, we can add "flags" `-l` for long listing format.

Flags (sometimes called options) allow us to finely control the behavior of the command. But how did we know to add `-l` after ls? The [`ls` manual ](https://man7.org/linux/man-pages/man1/ls.1.html) describes the command and all its options in the details. Like most commands, you can type the command followed `--help` to view the manual in your terminal.

```
ls --help
```

### Exercise 1

===  Question


	You can use multiple flags or options at the same time to modify the behavior of a command. What does the command `ls` do when used with:

	1. the `-h` option? 
	1. the `-l` and the `-h` option?
	1. the `-l`, the `-h`, and the `-F` option?

=== Answer

	1. The `-h` option makes the file size human readable, but it only noticable if you are printing the file size. 
	1. If you use both the `-h` option and the `-l` option (with`ls -lh` or `ls -l -h`), this makes the file size ‘human readable’, _i.e._ displaying something like 5.3K instead of 5369.
	1. The -F flag will class the file types by appending an identifier. This works best if there are directories present. 



What files are in the each of the different sub-directories? To list files in a different directory, you must specify the path or change directories. 

Because you can change directories using either the relative or absolute path, there multiple ways to successfully move up or down in the directory hierarchy.
Let's return to our home directory using the `cd` command and a relative path, then print the working directory to confirm.  



### Exercise 1

===  Question

	Starting from `books`, which of the following commands could Jovyan use to navigate to the `MiSeq` directory? 

	1. `cd MiSeq`
	2. `cd ./MiSeq`
	3. `cd ~/MiSeq`
	4. `cd ../MiSeq`
	5. `cd /MiSeq`
	6. `cd ../../Miseq`
	7. `cd /home/jovyan/MiSeq


=== Answer

1. No, MiSeq does not exist in the current working directory.
2. No, MiSeq does not exist in the current working directory.
3. Yes, MiSeq is in the home directory.
4. Yes, MiSeq is in the directory one level above.
5. No, MiSeq is not in the root directory.
6. No, MiSeq is not in the directory two levels above.
7. Yes, this is the full path to MiSeq.
:::

 
Let's practice using the cd and ls commands to explore files in different directories.  
 
### books

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

### southpark

How large are all the files in the `southparkdata` directory?

```
cd ~/books
ls -l 
```

We will see the following:

```
-rw-r--r-- 1 jovyan jovyan 1866609 Aug  4 04:34 All-seasons.csv.gz
-rw-r--r-- 1 jovyan jovyan     252 Aug  4 04:34 License.md
-rw-r--r-- 1 jovyan jovyan     345 Aug  4 04:34 README.md
```

### MiSeq


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