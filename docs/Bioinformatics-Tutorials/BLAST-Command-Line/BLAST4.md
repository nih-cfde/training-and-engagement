---
layout: page
title: How to Run BLAST+
---

## Download Data


First! We need some data.  Let's grab the mouse and zebrafish RefSeq
protein data sets from NCBI, and put them in /mnt, which is the
scratch disk or storage space for temp user files for Amazon machines.
```
   sudo chmod a+rwxt /mnt
   cd /mnt
```


## Running BLAST

We need some data!  Let's grab the mouse and zebrafish RefSeq
protein data sets from NCBI, and put them in our home directory.

Now, we'll use `curl` to download the files from a Web site onto our
computer; note, these files originally came from the
[NCBI FTP site](ftp://ftp.ncbi.nih.gov/refseq/M_musculus/mRNA_Prot)

```
curl -o mouse.1.protein.faa.gz -L https://osf.io/v6j9x/download
curl -o mouse.2.protein.faa.gz -L https://osf.io/j2qxk/download
curl -o zebrafish.1.protein.faa.gz -L https://osf.io/68mgf/download
```

If you look at the files in the current directory:

```
ls -l
```

You should now see these 3:

```
total 29876
-rw-rw-r-- 1 ubuntu ubuntu 12553742 Jul 16 17:06 mouse.1.protein.faa.gz
-rw-rw-r-- 1 ubuntu ubuntu  4074490 Jul 16 17:06 mouse.2.protein.faa.gz
-rw-rw-r-- 1 ubuntu ubuntu 13963093 Jul 16 17:06 zebrafish.1.protein.faa.gz
```

The three files you just downloaded are the last three on the list - the
`.faa.gz` files.

All three of the files are FASTA protein files (that's what the .faa
suggests) that are compressed with `gzip` (that's what the .gz means).

Uncompress them:

```
gunzip *.faa.gz
```

and let's look at the first few sequences in the file:

```
head mouse.1.protein.faa 
```

These are protein sequences in FASTA format.  FASTA format is something
many of you have probably seen in one form or another -- it's pretty ubiquitous.  

It's a text file, containing records; each record
starts with a line beginning with a '>', and then contains one or more lines of sequence text.

Let's take those first two sequences and save them to a file.  We'll
do this using output redirection with '>', which says "take
all the output and put it into this file here."

```
head -n 11 mouse.1.protein.faa > mm-first.faa
```

So now, for example, you can do `cat mm-first.faa` to see the contents of
that file (or `less mm-first.faa`). TIP: if you try `less mm-first.faa` you will need to exit by pressing the `q` key in your keyboard.

Now let's BLAST these two sequences against the entire zebrafish
protein data set. First, we need to tell BLAST that the zebrafish
sequences are (a) a database, and (b) a protein database.  That's done by calling 'makeblastdb':

```
makeblastdb -in zebrafish.1.protein.faa -dbtype prot
```

Next, we call BLAST to do the search:

```
blastp -query mm-first.faa -db zebrafish.1.protein.faa
```

This should run pretty quickly, but you're going to get a lot of output!!
To save it to a file instead of watching it go past on the screen,
ask BLAST to save the output to a file that we'll name `mm-first.x.zebrafish.txt`:

```
blastp -query mm-first.faa -db zebrafish.1.protein.faa -out mm-first.x.zebrafish.txt
```

and then you can 'page' through this file at your leisure by typing:

```
less mm-first.x.zebrafish.txt
```

(Type spacebar to move down, and 'q' to get out of paging mode.)

-----

Let's do some more sequences (this one will take a little longer to run):

```
head -n 498 mouse.1.protein.faa > mm-second.faa
blastp -query mm-second.faa -db zebrafish.1.protein.faa -out mm-second.x.zebrafish.txt
```

will compare the first 96 sequences.  You can look at the output file with:

```
less mm-second.x.zebrafish.txt
```

(and again, type 'q' to get out of paging mode.)

Notes:

* you can copy/paste multiple commands at a time, and they will execute in order;

* why did it take longer to BLAST ``mm-second.faa`` than ``mm-first.faa``?


----

Last, but not least, let's generate a more machine-readable version of that
last file --

```
blastp -query mm-second.faa -db zebrafish.1.protein.faa -out mm-second.x.zebrafish.tsv -outfmt 6
```

You can open the file with `less mm-second.x.zebrafish.tsv` to see how the file looks like.

See [this link](http://www.metagenomics.wiki/tools/blast/blastn-output-format-6) for a description of the possible BLAST output formats.

!!! Warning

    Remember to terminate your instance when you have completed your tasks with AWS.