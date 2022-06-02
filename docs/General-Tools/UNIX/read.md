---
title: Read
---


### `head`, `cat`, `less` and `more`

Now that we know what files exist on our computer, it's time to look at the contents of the file. There are multiple ways to look at the contents of a file. 

The `cat` command prints the entirety of a file to the stdout of our computer. We can scroll through files using the `less` command. Less is a safe way of looking at the contents of a file without the ability to change it. `head` prints, by default the first 10 lines of a file.


All three of the commands use the same syntax:

```bash
head [filename]
cat [filename]
less [filename]
```



!!! Tip

    You can use TAB to do filename completion, so if you type `cat R` and then press your Tab key once, it will autocomplete if there is a unique match. If there is more than one match, the first Tab will do nothing, and the second will show all the possible matches.



Let's navigate to the `books` directory and use the `head` command to view the `README.md` file. 

```bash
cd ~/books/
head README.md
```

You should see an output that looks like this. The `README.md` and `License.md` files are written in [Markdown](https://en.wikipedia.org/wiki/Markdown). To learn more about Markdown syntax, read this excellent [Markdown guide](https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet).

``` 
# Tale of Two Cities
# downloaded from [Project Gutenberg](https://www.gutenberg.org/ebooks/98)
```


Now we can view the file with `head`, `cat`, or `less` and `tail`. 

```bash
head book.txt
tail book.txt
cat book.txt
less book.txt
```


We can there are a lot more books, and we can look at the first few lines of all the txt files with the *. 

```bash
head *.txt
```

### `gunzip`

Notice, there is one book that is compressed. We can uncompress it with the command `gunzip`.

```bash
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