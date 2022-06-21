---
title: Redirect
---

In this section, we will explore redirecting, appending, piping, and looping with `>`, `>>`,  `|`, and `for`.

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

```bash
mkdir results
grep CATTAG *.fastq > results/files-with-CATTAG.txt
```

### `|`

Let's pipe (`|`) the output of a grep search to the command head to only show the first 10 lines. This is very handy when you only want to view a snapshot of the data. 


```bash
grep CATTAG *.fastq 
grep CATTAG *.fastq | wc -l
```

The pipe is a useful way to combine multiple commands on line line into a pipe. 

### For loops

However, the two lines above tell us how many times CATTAG appeared in all of the .fastq files in total. 
If we want to know who many times it occurs in each, we need a for loop. 

A for loop looks like this

```bash
for [thing] in [list of things]
do
command $[thing]
done
```


We could use this with `grep` on our fastq files like so. 


grep CATTAG *.fastq 

```bash
for file in *fastq
do
grep CATTAG $file | wc -l
done
```

If we cant to know which file these results came from, we can add an echo statement to echo the variable.


```bash
for file in *fastq
do
echo $file
grep CATTAG $file | wc -l
done
```

### history

If you want to save all the commands we used today, you can use the `history` command to print out all the commands you typed.

```bash
history
```

### `>`

You can redirect the output from the screen to a file using `>`. Note that `>` will overright existing conent, but `>>` will append. 

```bash
history > ~/history-2021-nov-17.txt
```

### Key points 

|Command|Description|
|-|-
|`\|` | pipes the standard output to a new command|
| `>`  | redirects the standard output to a new file |
| `>>`  | append the standard output to a new or existing file|
| `for` | initiates a for loop |
