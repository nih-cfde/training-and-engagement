# Creating and modifying files

We are used to copying and moving files using a GUI. These functions can be also carried out at the command line.

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

The `mv` command can be used to either move files to a new location or to rename them (which is essentially moving the contents from the old filename to the new file name. Let's use the `mv` command to rename one of the books.

```
mv A-tale-of-two-cities.txt Two-Cities.txt
```

If you want to move your book-copy.txt to a new folder, you first have to create that folder using the command `mkdir`. Then you move the copy. We can use `ls *` to list files in this directory and subdirectories.

```
mkdir book-copies
mv book-copy.txt book-copies
ls *
```

Take care when naming and renaming files. File names should not contain spaces or slashes. The use of capital letters (e.g. CamelCase) and underscores (e.g. snake_case) are often preferred over periods or spaces. It is good practice to keep track of where you got files. The commands used to get these books are stored in the `README.md`. The `README.md`  file is written in [Markdown](https://en.wikipedia.org/wiki/Markdown). (To learn more about Markdown syntax, read this excellent [Markdown guide](https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet).)

=== "Input"

	```
	head README.md
	```

=== "Expected Output"

	```
	# Books
	
	These books were downloaded from [Project 	Gutenberg](https://www.gutenberg.org/ebooks/) using the following commands. 
	
	
	curl https://www.gutenberg.org/files/98/98-0.txt -o book.txt
	curl https://www.gutenberg.org/files/98/98-0.txt -o A-tale-of-two-cities.txt
	curl https://www.gutenberg.org/files/11/11-0.txt -o Alice_in_wonderland.txt
	curl https://www.gutenberg.org/files/16/16-0.txt -o PeterPan.txt
	curl https://www.gutenberg.org/files/55/55-0.txt -o WizardOfOz.txt
	```

If the books are deleted or modified, they can easily be downloaded again from the source with `curl` command followed by the path to the file. Use the `-o` option to specify the file name. 

```
curl https://www.gutenberg.org/files/98/98-0.txt -o A-tale-of-two-cities.txt
```

=== "Exercise"

	Now you know how to copy and move files, but you may encounter errors if you try to move files to a directory that does not exist. But, have no fear, we can create new directories at the command line with the command `mkdir` followed by the path to the directories you want to create. 

	What happens when you run the following commands?


	1. `mkdir data results images/`
	2. `mkdir data/results/images`
	3. `mkdir -p data/results/images`


=== "Answer"

	1. 3 subfolders (data, results, and images) all created in the working directory. 
	2. An error message (`mkdir: cannot create directory ‘data/results/images’: No such file or directory`) appears because the results directory does not exist
	3. The `-p` option tells `mkdir` to create parent directories if they do not already exist, so the subdirectory results and its subdirectory, images, are successfully created.


!!! note "Key Points"

	| Command | Description |
	| -------- | -------- | 
	|cp [old] [new] | copies a file | 
	|mv [old] [new]| moves or renames a file or directory |
	|rm [path] | removes (deletes) a file |
	|mkdir -p [path/to/files] | creates a hierarchy of directories |
	|rmdir [path] | removes an empty  directory |
	|curl [webaddress] -o [filename] | downloads from URL named [filename]