---
title: Read
---


### `head`, `cat`, `less`, and `tail`

Now that we know what files exist on our computer, it's time to look at the contents of the file. There are multiple ways to look at the contents of a file. 

The `cat` command prints the entirety of a file to the stdout of our computer. We can scroll through files using the `less` command. Less is a safe way of looking at the contents of a file without the ability to change it. `head` prints, by default the first 10 lines of a file and `tail` prints the last 10 lines.


All four of the commands use the same syntax:

```bash
head [filename]
cat [filename]
less [filename]
tail [filename]
```



!!! Tip

    You can use TAB to do filename completion, so if you type `cat R` and then press your TAB key once, it will autocomplete if there is a unique match. If there is more than one match, the first TAB will do nothing, and the second will show all the possible matches.



Let's navigate to the `books` directory and use the `head` command to view the `README.md` file. 

```bash
cd ~/books/
head README.md
```

You should see an output that looks like this. The `README.md` file is written in [Markdown](https://en.wikipedia.org/wiki/Markdown). To learn more about Markdown syntax, read this excellent [Markdown guide](https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet).

``` 
# Books

These books were downloaded from [Project Gutenberg](https://www.gutenberg.org/ebooks/) using the following commands. 

curl https://www.gutenberg.org/files/98/98-0.txt -o book.txt
curl https://www.gutenberg.org/files/98/98-0.txt -o A-tale-of-two-cities.txt
curl https://www.gutenberg.org/files/11/11-0.txt -o Alice_in_wonderland.txt
curl https://www.gutenberg.org/files/16/16-0.txt -o PeterPan.txt
curl https://www.gutenberg.org/files/55/55-0.txt -o WizardOfOz.txt
```


Now we can view the file with `head`, `cat`, or `less` and `tail`. 

```bash
head book.txt
cat book.txt
less book.txt
tail book.txt
```


We can see there are several more books in the directory, and we can look at the first few lines of all the txt files with the *. 

```bash
head *.txt
```

### `gunzip`

Notice, there is one book that is compressed. We can uncompress it with the command `gunzip`.

```bash
gunzip WizardOfOz.txt.gz
```

Now the `ls` command will show that the `WizardOfOz.txt.gz` has been replaced with the unzipped `WizardOfOz.txt` file. `gunzip` also has a number of flags you can use including `-k` which will allow you to unzip the file and keep the original.

### Key points

| Command [OPTION] | Description |
| -------- | -------- | 
|`head [filename]` | print first 10 lines of `FILENAME` | 
|`cat [filename]`| print `FILENAME`'s contents to stdout|
|`less [filename]`|view `FILENAME` without printing  to stdout |
|`tail [filename]` | print last 10 lines of `FILENAME` |
|`gunzip -k [filename]` | uncompress a file and keep the original |

