# How to run a BLAST search

## Step 1: Download data

First, we need some data!  For this tutorial, we're using mouse and zebrafish RefSeq protein data sets from NCBI. We'll put them in our home directory:
```
sudo chmod a+rwxt /mnt
cd /mnt
```

!!! Note

	`chmod` command stands for **change mode** and it is used to define or change permissions or modes on files and limit access to only those who are allowed access. It's the same as using your mouse to right-click a file or folder and selecting the permission tabs and defining who can access the resource.

Now, we'll use `curl` to download the files from an online OSF repository to the AWS remote instance. These are copies of the files from the [NCBI FTP site](ftp://ftp.ncbi.nih.gov/refseq/M_musculus/mRNA_Prot).

!!! Note
	You can Copy/Paste multiple commands at a time, and they will execute in order.

```
curl -o mouse.1.protein.faa.gz -L https://osf.io/v6j9x/download
curl -o zebrafish.1.protein.faa.gz -L https://osf.io/68mgf/download
```

!!! Note
	`curl` is a command line tool that allows you to transfer data from or to a remote server. The lowercase `-o` allows you to specifiy the name of the saved file, while the uppercase `-O` saves the file with its original filename. The `-L` is the location of the file.

There should now be two protein files in the current directory:

=== "Input"
	```
	ls -lht
	```
	Use the `ls` command to list the files in the directory. The `-l` flag adds file permissions, ownership, size, and creation time information. The `-h` flag converts the file sizes into human-readable units and the `-t` flag sorts the files starting with the most recently created file.

=== "Expected Output"

	Your output should look similar, though the download dates will be the current date and time when you downloaded the files:
	```
	total 30M
	-rw-rw-r-- 1 ubuntu ubuntu  14M Oct  2 20:30 zebrafish.1.protein.faa.gz
	-rw-rw-r-- 1 ubuntu ubuntu  12M Oct  2 20:28 mouse.1.protein.faa.gz
	```

These files are FASTA protein files (that's what the ".faa"
means) that are compressed with `gzip` (that's what the ".gz" means).

Uncompress them:

```
gunzip *.faa.gz
```

Let's look at the first few sequences in the file:

```
head mouse.1.protein.faa
```

These are protein sequences in FASTA format.  FASTA format is something
many of you have probably seen in one form or another -- it's ubiquitous!

It's a text file, containing sequence records. Each record
starts with a `>`, and then contains one or more lines of sequence text. The information associated with each
sequence is in the header line: accession number, gene, and species.

## Step 2: Subset data for demo

Let's take those first two sequences and save them to a file.  We'll
do this using output redirection with `>`, which says "take
all the output from the preceding command and put it into this new file."

Select the first 11 lines of the file and save as a new file:

=== "Input"

	```
	head -n 11 mouse.1.protein.faa > mm-first.faa
	```

=== "Expected Output"

	So now, for example, you can do `cat mm-first.faa` to see the contents of
	that file (or `less mm-first.faa`):

	```
	ubuntu@ip-172-31-17-217:~$ cat mm-first.faa
	>YP_220550.1 NADH dehydrogenase subunit 1 (mitochondrion) [Mus musculus domesticus]
	MFFINILTLLVPILIAMAFLTLVERKILGYMQLRKGPNIVGPYGILQPFADAMKLFMKEPMRPLTTSMSLFIIAPTLSLT
	LALSLWVPLPMPHPLINLNLGILFILATSSLSVYSILWSGWASNSKYSLFGALRAVAQTISYEVTMAIILLSVLLMNGSY
	SLQTLITTQEHMWLLLPAWPMAMMWFISTLAETNRAPFDLTEGESELVSGFNVEYAAGPFALFFMAEYTNIILMNALTTI
	IFLGPLYYINLPELYSTNFMMEALLLSSTFLWIRASYPRFRYDQLMHLLWKNFLPLTLALCMWHISLPIFTAGVPPYM
	>YP_220551.1 NADH dehydrogenase subunit 2 (mitochondrion) [Mus musculus domesticus]
	MNPITLAIIYFTIFLGPVITMSSTNLMLMWVGLEFSLLAIIPMLINKKNPRSTEAATKYFVTQATASMIILLAIVLNYKQ
	LGTWMFQQQTNGLILNMTLMALSMKLGLAPFHFWLPEVTQGIPLHMGLILLTWQKIAPLSILIQIYPLLNSTIILMLAIT
	SIFMGAWGGLNQTQMRKIMAYSSIAHMGWMLAILPYNPSLTLLNLMIYIILTAPMFMALMLNNSMTINSISLLWNKTPAM
	LTMISLMLLSLGGLPPLTGFLPKWIIITELMKNNCLIMATLMAMMALLNLFFYTRLIYSTSLTMFPTNNNSKMMTHQTKT
	KPNLMFSTLAIMSTMTLPLAPQLIT
	```

!!! tip

	if you try `less mm-first.faa`, exit back to the terminal by pressing the ++q++ key in your keyboard.

## Step 3: Create blast database

Now let's BLAST these two mouse sequences against the entire zebrafish
protein data set. First, we need to tell BLAST that the zebrafish
sequences are (a) a database, and (b) a protein database.  That's done by calling `makeblastdb`:

=== "Input"

	```
	makeblastdb -in zebrafish.1.protein.faa -dbtype prot
	```

=== "Expected Output"

	```
	Building a new DB, current time: 10/02/2020 20:38:17
	New DB name:   /home/ubuntu/zebrafish.1.protein.faa
	New DB title:  zebrafish.1.protein.faa
	Sequence type: Protein
	Keep MBits: T
	Maximum file size: 1000000000B
	Adding sequences from FASTA; added 53088 sequences in 1.69959 seconds.
	```

## Step 4: Run a BLAST search

We'll use the `blastp` command to conduct a protein sequence search comparing the mouse sequences (`-query`) with the zebrafish database (`-db`):

```
blastp -query mm-first.faa -db zebrafish.1.protein.faa
```

This should run pretty quickly, but you're going to get a lot of output!!
To save it to a file instead of watching it go past on the screen,
ask BLAST to save the output to a file that we'll name `mm-first.x.zebrafish.txt`:


```
blastp -query mm-first.faa -db zebrafish.1.protein.faa -out mm-first.x.zebrafish.txt
```

Now you can page through this file by typing:

```
less mm-first.x.zebrafish.txt
```

Use the ++up++ and ++down++ keys to move up/down, and ++q++ to exit the paging mode.


Let's run some more sequences (this search will take a little longer to run).
Using the `head -n` command, we'll select the first 498 lines of the mouse protein FASTA file.
Since these FASTA files are divided over five lines per sequence, this subset will include 96 sequences.

```
head -n 498 mouse.1.protein.faa > mm-second.faa
blastp -query mm-second.faa -db zebrafish.1.protein.faa -out mm-second.x.zebrafish.txt
```

Check the output file:
```
less mm-second.x.zebrafish.txt
```

!!! tip

	We can verify the number of records by searching for all of the lines that include the `>` symbol.
	Type
	```
	grep '>' mm-second.faa | wc -l
	```

	`grep` is a powerful tool for finding patterns and words on Unix and Linux Systems. `wc -l` tells `wc` to count the number of lines that have the matching pattern `>`. The command should tell you there are 96 lines in this file.


=== "Exercise"

	Why did it take longer to BLAST `mm-second.faa` than `mm-first.faa`?

=== "Answer"

	`mm-second.faa` has 96 sequences in comparison to `mm-first.faa` which only had 2 sequences.


Last, but not least, let's generate a more machine-readable version of that
last file. There are several BLAST output formats. The [blast output format 6](http://www.metagenomics.wiki/tools/blast/blastn-output-format-6)
is the most commonly used.

=== "Input"

	```
	blastp -query mm-second.faa -db zebrafish.1.protein.faa -out mm-second.x.zebrafish.tsv -outfmt 6
	```

=== "Expected Output"

	```
	less mm-second.x.zebrafish.tsv
	```

	The first few lines look like this:
	```
	YP_220550.1     NP_059331.1     69.010  313     97      0       4       316     10      322     1.24e-150       426
	YP_220551.1     NP_059332.1     44.509  346     188     3       1       344     1       344     8.62e-92        279
	YP_220551.1     NP_059341.1     24.540  163     112     3       112     263     231     393     5.15e-06        49.7
	```

For this particular example, the blast search results show that the sequences we compared are quite different. The top match from the mouse data was only 69% similar to the zebrafish reference database. Of course, this makes sense since mouse and zebrafish species are very divergent!

!!! Warning

	Remember to download any files you want to save and terminate your instance when you have completed your tasks with AWS. See the AWS [tutorial](../Introduction_to_Amazon_Web_Services/introtoaws1.md) for help.

!!! note "Key Points"

	BLAST searches can be conducted between DNA and protein sequences, including translating DNA to protein or protein to DNA. For this tutorial, we used the `blastp` command to compare protein sequences. Other commands include `blastn` for comparing DNA sequences and `blastx` or `tblastn` for comparing translated sequences. This flexibility and ease of use has made BLAST a very popular tool for quickly comparing sequences. Command line BLAST searches can also be scaled up to compare sequences for multiple samples to large databases. For more information about BLAST, check out the [NCBI webpage](https://blast.ncbi.nlm.nih.gov/Blast.cgi).
