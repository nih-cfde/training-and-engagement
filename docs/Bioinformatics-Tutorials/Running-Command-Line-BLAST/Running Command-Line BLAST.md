# Running Command-Line BLAST

Learning Objectives

- Gain hands-on exposure to the linux command line
- Understand how data is turned into results by programs run at the command line


## Getting Started 

<iframe id="kaltura_player" src="https://cdnapisec.kaltura.com/p/1770401/sp/177040100/embedIframeJs/uiconf_id/29032722/partner_id/1770401?iframeembed=true&playerId=kaltura_player&entry_id=0_q6gmtglm&flashvars[mediaProtocol]=rtmp&amp;flashvars[streamerType]=rtmp&amp;flashvars[streamerUrl]=rtmp://www.kaltura.com:1935&amp;flashvars[rtmpFlavors]=1&amp;flashvars[localizationCode]=en&amp;flashvars[leadWithHTML5]=true&amp;flashvars[sideBarContainer.plugin]=true&amp;flashvars[sideBarContainer.position]=left&amp;flashvars[sideBarContainer.clickToClose]=true&amp;flashvars[chapters.plugin]=true&amp;flashvars[chapters.layout]=vertical&amp;flashvars[chapters.thumbnailRotator]=false&amp;flashvars[streamSelector.plugin]=true&amp;flashvars[EmbedPlayer.SpinnerTarget]=videoHolder&amp;flashvars[dualScreen.plugin]=true&amp;flashvars[Kaltura.addCrossoriginToIframe]=true&amp;&wid=0_qy8cnqw9" width="1000" height="500" allowfullscreen webkitallowfullscreen mozAllowFullScreen allow="autoplay *; fullscreen *; encrypted-media *" sandbox="allow-forms allow-same-origin allow-scripts allow-top-navigation allow-pointer-lock allow-popups allow-modals allow-orientation-lock allow-popups-to-escape-sandbox allow-presentation allow-top-navigation-by-user-activation" frameborder="0" title="Kaltura Player"></iframe>

## What is BLAST?

BLAST is the **B**asic **L**ocal **A**lignment **S**earch large sequence databases; It starts by finding small matches between the two sequences and extending those matches.  For in-depth information on how BLAST works and the different BLAST functionality, check out the resources page [here](https://blast.ncbi.nlm.nih.gov/Blast.cgi).

BLAST can be helpful for identifying the source of a sequence or finding a similar sequence in another organism.  In this lesson, we will use BLAST to find zebrafish proteins that are similar to a small set of mouse proteins.

## Why use the command line?
BLAST has a very nice graphical interface for searching sequences in NCBI's database.
However, running BLAST through the commmand line has many benefits:


  * It's much easier to run many BLAST queries using the command line than the GUI
  * Running BLAST with the command line is reproducible and can be documented in a script
  * The results can be saved in a machine-readable format that can be analyzed later on
  * You can create your own databases to search rather than using NCBI's pre-built databases
  * It allows the queries to be automated
  * It allows you to use a remote computer to run the BLAST queries
  
Later on in the workshop we will talk more about these advantages and have a more in-depth explanation of the shell.

## Set-up an Instance

Boot a t2.micro instance on AWS and connect your shell prompt.

!!! Tip 
	
	If you need assistance setting up an instance check out the tutorial Intro to AWS!


## Install BLAST software using Command Line
Now, install some software. We will need NCBI BLAST for the below tutorial


Copy and paste the following commands
```
   sudo apt-get update && sudo apt-get -y install python ncbi-blast+
```

This updates the software list and installs the Python programming
language and NCBI BLAST+.

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
