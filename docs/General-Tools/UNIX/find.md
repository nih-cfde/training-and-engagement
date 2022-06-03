---
title: Find
---


A big part of data science is making sure what you _expect_ in a particular file is what you _have_ in that file. This is fairly easy when your files are small but is challenging when the files are much larger than your screen.


To explore this topic in more detail, navigate to the `data/MiSeq/` directory.


```bash
cd ~/MiSeq
ls
```

### FASTQ format

This directory contains multiple [**FASTQ** files](https://en.wikipedia.org/wiki/FASTQ_format). A FASTQ file normally uses four lines per sequence.

* Line 1 begins with a '@' character and is followed by a sequence identifier and an optional description (like a FASTA title line).
* Line 2 is the raw sequence letters.
* Line 3 begins with a '+' character and is optionally followed by the same sequence identifier (and any description) again.
* Line 4 encodes the quality values for the sequence in Line 2, and must contain the same number of symbols as letters in the sequence.

A FASTQ file containing a single sequence might look like this:

!!! info

     An example FASTQ file
```
@SEQ_ID
GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65
```

We can use the `cat` command to print fastq files to the screen, but thousands of lines of text would crowd your screen. Instead, let's use the `head` command to view the first 8 lines file. You can copy the file name and paste it into the console or you can type and use tab complete to pick a particular file. 

```bash
head -n 4 F3D0_S188_L001_R1_001.fastq
```

```
@M00967:43:000000000-A3JHG:1:1101:18327:1699 1:N:0:188
NACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGCCTGCCAAGTCAGCGGTAAAATTGCGGGGCTCAACCCCGTACAGCCGTTGAAACTGCCGGGCTCGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACCCCGATTGCGAAGGCAGCATACCGGCGCCCTACTGACGCTGAGGCACGAAAGTGCGGGGATCAAACAG
+
AABABBFFFGGGGGGGGGGGGGGGGHHHHHHHGGGHHHHHGHGGGGGGGHGGGGGGHHHHHHHHHHGGGGGHHHHGHGGGGGGHHBGHGDGGGGGHHHGGGGHHHHHHHHGGGGGHG@DHHGHEGGGGGGBFGGEGGGGGGGG.DFEFFFFFFFDCFFFFFFFFFFFFFFFFFFFFFFFFFFDFDFFFEFFCFF?FDFFFFFFFFAFFFFFFFFFFFBDDFFFFFEFADFFFFFBAFFFA?EFFFBFF
@M00967:43:000000000-A3JHG:1:1101:14069:1827 1:N:0:188
TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGCCTGCCAAGTCAGCGGTAAAATTGCGGGGCTCAACCCCGTACAGCCGTTGAAACTGCCGGGCTCGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACCCCGATTGCGAAGGCAGCATACCGGCGCCCTACTGACGCTGAGGCACGAAAGTGCGGGGATCAAACAG
+
3AA?ABBDBFFBEGGEGGGGAFFGGGGGHHHCGGGGGGHFGHGGCFDEFGGGHGGGEGF1GGFGHHHHHGGEGGHHHHHFGGGGGGHHHHHGGGGCDDGHHGGGFHHHHHHHHCD@CCHGGGGHEHGGG@GFGGGGGGG@BGGGEGCEBFFFBFFB;9@EFFFEFFFFFFFFFFFFAFBBBFFFFFBBBFFFFBBBFFFFFFFFFFFBBBBBBBFFFFFFFFFDDFAFFFFF.AF9/FBBBBB.EAFFE?F
```

`head` prints the first ten lines of a file out onto your screen. Similarly, the `tail` command prints the last 10 lines of a file. 

```bash
tail -4 F3D0_S188_L001_R1_001.fastq
```

```
@M00967:43:000000000-A3JHG:1:1105:19125:28016 1:N:0:188
TACGTAGGGGGCAAGCGTTATCCGGAATTACTGGGTGTAAAGGGAGCGTAGACGGTAATGCAAGTCTGGAGTGAAAGGCGGGGGCCCAACCCCCGGACTGCTCTGGAAACTGTGTAACTGGAGTGCAGGAGAGGCAGGCGGAATTCCTAGTGTAGCGGTGAAATGCGTAGATATTAGGAGGAACACCAGTGGCGAAGGCGGCCTGCTGGACTGTAACTGACGTTGAGGCTCGAAAGCGGGGGGGGCAAAAA
+
>AAAAFFBBBBDGFGGGEGFGGAEEEGGHHGHHHHFFEHHHHGGHGGGGFGGHGGGEGHHGHHHHHHHHGHGHFFHHHHGGGGG@DGDGDGGGGGGDDGHHHHGHHCFHGHHHH0GGFFHHHGFFEFGGGGGGGGGG@.BB-@EF/BFFF9//9FF;B.;.:F//9A.;9.9;B/////BFDB9FF.AAB/:?B-;-99@B-9;EA/;B99BE99B///B:9.;:9F..;9=-FD.9A-;@BB----?F..
@M00967:43:000000000-A3JHG:1:1105:20429:28046 1:N:0:188
TACGGAGGATTCAAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGTTCGATAAGTTAGAGGTGAAATCCCGGGGCTCAACTCCGGCACTGCCTCTGATACTGTCGGGCTAGAGTTTAGTTGCGGTAGGCGGAATGTATGGTGTAGCGGTGAAATGCATAGAGATCATACAGAACACCGATTGCGAAGGCAGCTTACCAAACTACGACTGACGTTGAGGCACGAAAGCGTGGGGAGCAAACA
+
BBBBBBBBBFFFGGGFFGGGGGGGGGGGHHHHHHHGGGHGHHHGHGGGGGGGHGGEGGFGGHGHHHHHHHHHGHHHHEHHGGGGGGHHHHHHHGGGGGHHHGFHHHHGGHHHHHGGGGGHFGH?GHHGHHHHGG<CDFHGGGGGHCFGHHHHEHFHCGCGEFFGGGGGGEGGEFGGGFFFFGGG/FDCFFFFFADFFFFFDFFFFFFFFFABFFFFF?DFFFE.EEFFFEFFFFAAFADFFFFECDEFB..
```

### FASTA format

**FASTQ** files should not be confused with **FASTA** files. FASTQ files contain information about the quality of the sequence, but FASTA files only contain the sequence and an identifier.

!!! info

     An Example FASTA file

```
> SEQUENCE_1
MTEITAAMVKELRESTGAGMMDCKNALSETNGDFDKAVQLLREKGLGKAAKKADRLAAEG
LVSVKVSDDFTIAAMRPSYLSYEDLDMTFVENEYKALVAELEKENEERRRLKDPNKPEHK
IPQFASRKQLSDAILKEAEEKIKEELKAQGKPEKIWDNIIPGKMNSFIADNSQLDSKLTL
MGQFYVMDDKKTVEQVIAEKEKEFGGKIKIVEFICFEVGEGLEKKTEDFAAEVAAQL
>SEQUENCE_2
SATVSEINSETDFVAKNDQFIALTKDTTAHIQSNSLQSVEELHSSTINGVKFEEYLKSQI
ATIGENLVVRRFATLKAGANGVVNGYIHTNGRVGVVIAAACDSAEVASKSRDLLRQICMH
```


Let's look at a synthetic FASTA file. 

```bash
head -4 HMP_MOCK.v35.fasta 
```


```
> A.baumannii.1 
TGGGGAATATTGGACAATGGGGGGAACCCTGATCCAGCCATGCCGCGTGTGTGAAGAAGGCCTTATGGTTGTAAAGCACTTTAAGCGAGGAGGAGGCTACTTTAGTTAATACCTAGAGATAGTGGACGTTACTCGCAGAATAAGCACCGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGGTGCGAGCGTTAATCGGATTTACTGGGCGTAAAGCGTGCGTAGGCGGCTTATTAAGTCGGATGTGAAATCCCCGAGCTTAACTTGGGAATTGCATTCGATACTGGTGAGCTAGAGTATGGGAGAGGATGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGATGGCGAAGGCAGCCATCTGGCCTAATACTGACGCTGAGGTACGAAAGCATGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCATGCCGTAAACGATGTCTACTAGCCGTTGGGGCCTTTGAGGCTTTAGTGGCGCAGCTAACGCGATAAGTAGACCGCCTGGGGAGTACGGTC
> A.odontolyticus.1
TGGGGAATATTGCACAATGGGCGAAAGCCTGATGCAGCGACGCCGCGTGAGGGATGGAGGCCTTCGGGTTGTAAACCTCTTTCGCTCATGGTCAAGCCGCAACTCAAGGTTGTGGTGAGGGTAGTGGGTAAAGAAGCGCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGCGCGAGCGTTGTCCGGAATTATTGGGCGTAAAGGGCTTGTAGGCGGTTGGTCGCGTCTGCCGTGAAATCCTCTGGCTTAACTGGGGGCGTGCGGTGGGTACGGGCTGACTTGAGTGCGGTAGGGGAGACTGGAACTCCTGGTGTAGCGGTGGAATGCGCAGATATCAGGAAGAACACCGGTGGCGAAGGCGGGTCTCTGGGCCGTTACTGACGCTGAGGAGCGAAAGCGTGGGGAGCGAACAGGATTAGATACCCTGGTAGTCCACGCTGTAAACGTTGGGCACTAGGTGTGGGGGCCACCCGTGGTTTCTGCGCCGTAGCTAACGCTTTAAGTGCCCCGCCTGGGGAGTACGGCC
```

FASTQ and FASTA files are often used in combination to map reads to a genome or transcription. You've seen three ways to read large files. Next, we will learn how to use and manipulate the files.  


## Wildcards

Sometimes you know a file or directory exists, but you can't find it. Sometimes you want to find many files with similar properties. This is where the wildcard (`*`) comes in handy.  What do the following commands do?

 
1. `ls *` 
1. `ls MiSeq/F3D*`
1. `ls MiSeq/*fasta`
 

!!! spoiler

     1. `ls *` lists files in the working directory and 1 level down. 
1. `ls MiSeq/F3D*` lists files in the data/MiSeq directory that start with "F3D".
1. `ls MiSeq/*fasta` lists files in the data/MiSeq directory that end with "fasta".



### `grep`

A lot of the time we want to know if a file contains what we expect. A useful thing to do is to be able to **search the contents of files** for a particular string of characters you would like to find.  We can use the file pattern searcher `grep` to find things.

The `MiSeq/` directory contains many of the sequence files ending in`.fastq`. We expect these files to contain information in a particular format throughout the file with four lines of information for each sequence string. Looking through a million line file using less will take a long time. Rather than manually looking at the whole file, we can print only a portion of the file's contents to standard output. 

Let's say you'd like to find the sequence `CATTAG` in your MiSeq files. We can use function `grep` to search for  `CATTAG` in one or all all of the fastq files located in our current working directory.

```
cd ../data/MiSeq/
grep CATTAG F3D0_S188_L001_R2_001.fastq
grep CATTAG *.fastq
```


=== "Challenge"

What line does `CATTAG` occur on in `F3D141_S207_L001_R1_001.fastq`? 

=== "Hint"

Use `grep --help` to search for `grep` options related to line number.
`grep -n [filename]`` will print the line number.


### `wc`

How many lines are in the file total? The word count (`wc`) command will give us that answer.

```bash
wc F3D0_S188_L001_R2_001.fastq
```

By default, `wc` prints the characters, words, and lines in a file. To extract just the line numbers, we use `wc -l`. We will use this more later. 


```bash
wc -l F3D0_S188_L001_R2_001.fastq
```


### Key points

|Command|Description|
|-|-
|`find [filename]` | finds files with specific properties that match patterns|
| `grep [option] [filename]`  | selects lines in files that match patterns|
| `wc` [filename] | will print the total characters, words, and lines in a file.


