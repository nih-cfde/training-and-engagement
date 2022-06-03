---
title: Work
---



We are quite used to working with files using a graphical user interface (or GUI). In this section, you will learn how to copy, move, create, and delete directories and files.

The `cp` and `mv` commands can be used to copy and move (or rename) files and directories respectively. For both commands, you must specify the old and new names. Specifying the path is necessary if you want to move files out of the current working directory.


```bash
cp [original-filename] [copy-filename]
mv [original-filename] [new-filename]
```

### `cp` 

Let's make a copy of some raw data before we start modifying it. We will use `ls` to check our work.

```bash
cp book.txt book-copy.txt
ls
```

###  `mv`

The `mv` command can be used to either move files to a new location or to rename them (which is essentially moving the contents from the old filename to the new file name. Let's use the `mv` command to rename the copied and compressed file back to the original name.

```bash
mv book-copy.txt book-2cities.txt

```

!!! warning

Take care when renaming files. It is good practice to keep track of changes in file names and links to the source data. 

### `mkdir` 

Now you know how to copy and move files, but you may encounter errors if you try to move files to a directory that does not exist. But, have no fear, we can create new directories at the command line with the command `mkdir` followed by the path to the directories you want to create. 

What happens when you run the following commands?

```bash
mkdir data results images/
mkdir -p data/results/images
```
The first line creates two directories, results and images. The second line also create two directories, but they are nested with the parent directory data. 
The `-p` argument creates parent directories if they do not already exist.
 
### `rmdir` 

If you created some files or directories that you do not want, you can remove them with the `rm` and `rmdir` commands. 
`rmdir` will only remove empty directories, but `rm -r` will remove recursively.

```bash
rmdir data/results/images
rmdir data/results
```

or 

```bash
rm -r data/results
```

### Key points

| Command | Description |
| -------- | -------- | 
|cp [old] [new] | copies a file | 
|mv [old] [new]| moves or renames a file or directory |
|rm [path] | removes (deletes) a file |
|mkdir [path] | creates a new directory |
|rmdir [path] | removes an empty  directory |
