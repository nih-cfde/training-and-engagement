# Concluding thoughts

This lesson focused on file and directory exploration because that’s something everyone needs to know, and all these commands will work on pretty much any computer that is running a UNIX compatible shell (including Mac OS X and Windows Subsystem for Linux).

We have shown you multiple options for editing and working with text files. These tools may seem confusing at first, but they will become second nature if you use them regularly.

If you want to save all the commands we used today, you can use the history command to print out all the commands you typed.

If you want to save all the commands we used today, you can use the `history` command to print out all the commands you typed.

```
history
```

You can save the file in your home directory with: 

```
history > ~/history.txt
```

You can also pipe your history to `grep` to focus in on frequently used commands such as:

```
history | grep "ls"
history | grep "grep"
history | grep "for"
``` 

=== "Tutorial Resources"

The binder [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/nih-cfde/training-rstudio-binder/data?urlpath=rstudio) and the lesson documentation pages will stay working for the foreseeable
future, so please feel free to come back and revisit some of these commands!

Google (and especially stackoverflow) is your friend! Use Internet
search whenever you have questions about what a command does, or what
commands to use to achieve a particular task.

Please refer to our [UNIX Cheatsheet](https://training.nih-cfde.org/en/latest/General-Tools/Cheat-Sheets/bash_cheatsheet/) for a list of commonly used commands.

The website [Explain Shell](https://explainshell.com/) searches the help text and prints the information for a command line entered. 

The lesson materials were adapted from the UC Davis Data Lab's [Intro to Cloud Computing](https://ngs-docs.github.io/2021-august-remote-computing/) workshop, Data Carpentry's [Introduction to the Command Line for Genomics](https://datacarpentry.org/shell-genomics/) lesson, and the Lab for Data Intensive Biology's [Advanced Beginner/Intermediate Shell](https://dib-training.readthedocs.io/en/pub/2016-01-13-adv-beg-shell.html) workshop. 


!!! note "Key Points"
	This workshop teaches a dozen commonly used UNIX commands that can be combined to perform power, reproducible bioinformatic workflows. The commands taught `pwd` `ls`  `cd` `cat` `head` `less` `cp` `mv` `rm` `mkdir` `grep` `wc` `cut` `gunzip` and `gzip` (and probably a few others).

