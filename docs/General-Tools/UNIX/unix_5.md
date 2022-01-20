# Finding things

A big part of data science is making sure what you _expect_ in a particular file is what you _have_ in that file. This is fairly easy when your files are small but is challenging when the files are much larger than your screen. 




#### Books 

Let's navigate to the books directory to explore this topic. 


```
cd ~/books
pwd
```


Earlier, when we used `head` and `tail` to view the books, we read extraneous information about Project Gutenberg. Let's say you want to extract some information about each book such as the author or title? We can use the command `grep` to search for a pattern

=== "Input"

	```
	grep Title *txt
	```

=== "Expected Output"

	```
	Alice_in_wonderland.txt:Title: Alice’s Adventures in Wonderland
	A-tale-of-two-cities.txt:Title: A Tale of Two Cities
	book.txt:Title: A Tale of Two Cities
	PeterPan.txt:Title: Peter Pan
	WizardOfOz.txt:Title: The Wonderful Wizard of Oz
	```

What are the Chapter titles of each book?

=== "Input"

	```
	grep Chapter *txt
	```

=== "Expected Output"

	You should see some text including the following:

	```
	PeterPan.txt: Chapter I. PETER BREAKS THROUGH
	PeterPan.txt: Chapter II. THE SHADOW
	...
	WizardOfOz.txt: Chapter I. The Cyclone
	WizardOfOz.txt: Chapter II. The Council with the Munchkins

	```

Why do we only see Chapter titles for 2 books? Because the other Chapter titles are written in all caps. Let's modify grep to match Chapter 1, chapter 1, and CHAPTER 1. We use the `-i` option to "ignore case" and the `-w` option to match the word. 

=== "Input"

	```
	grep -i -w "chapter i" *txt
	```

=== "Expected Output"


	This returns the first chapter for each book. Shown below is the  first chapter for each book. This pattern occurs once in the table of contents and once in the main text. A Tale of Two Cities has 3 Chapter 1s. 


	```
	Alice_in_wonderland.txt: CHAPTER I.     Down the Rabbit-Hole
	Alice_in_wonderland.txt:CHAPTER I.
	A-tale-of-two-cities.txt:     CHAPTER I      The Period
	A-tale-of-two-cities.txt:     CHAPTER I      Five Years Later
	A-tale-of-two-cities.txt:     CHAPTER I      In Secret
	A-tale-of-two-cities.txt:CHAPTER I.
	A-tale-of-two-cities.txt:CHAPTER I.
	A-tale-of-two-cities.txt:CHAPTER I.
	book.txt:     CHAPTER I      The Period
	book.txt:     CHAPTER I      Five Years Later
	book.txt:     CHAPTER I      In Secret
	book.txt:CHAPTER I.
	book.txt:CHAPTER I.
	book.txt:CHAPTER I.
	PeterPan.txt: Chapter I. PETER BREAKS THROUGH
	PeterPan.txt:Chapter I.
	WizardOfOz.txt: Chapter I. The Cyclone
	WizardOfOz.txt:Chapter I
```




#### MiSeq

To explore this topic in more detail and in a biological context, navigate to the `data/MiSeq/` directory.


```
cd ~/MiSeq
ls
```

This directory contains multiple [**FASTQ** files](https://en.wikipedia.org/wiki/FASTQ_format)). A FASTQ file normally uses four lines per sequence.

* Line 1 begins with a '@' character and is followed by a sequence identifier and an optional description (like a FASTA title line).
* Line 2 is the raw sequence letters.
* Line 3 begins with a '+' character and is optionally followed by the same sequence identifier (and any description) again.
* Line 4 encodes the quality values for the sequence in Line 2 and must contain the same number of symbols as letters in the sequence.

A FASTQ file containing a single sequence might look like this:


!!! info

An example FASTQ file
```
@SEQ_ID
GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65
```


We can use the `cat` command to print .fastq files to the screen, but thousands of lines of text would crowd your screen. Instead, let's use the `head` and `less` command to view a portion of a file.  By default, they print the first 10 lines. We can use the `-n` flag to specify how many lines to print. Because each entry of a .fastq file consists of 4 lines, printing the first 4 and last 4 lines of the file will confirm that the file is properly formatted.

You can copy the file name and paste it into the console or you can type and use tab complete to pick a particular file. 


=== "Input"
	
	```
	head -n 4 F3D0_S188_L001_R1_001.fastq
	```

=== "Expected Output"

	```
	@M00967:43:000000000-A3JHG:1:1101:18327:1699 1:N:0:188
	NACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGCCTGCCAAGTCAGCGGTAAAATTGCGGGGCTCAACCCCGTACAGCCGTTGAAACTGCCGGGCTCGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACCCCGATTGCGAAGGCAGCATACCGGCGCCCTACTGACGCTGAGGCACGAAAGTGCGGGGATCAAACAG
	+
	AABABBFFFGGGGGGGGGGGGGGGGHHHHHHHGGGHHHHHGHGGGGGGGHGGGGGGHHHHHHHHHHGGGGGHHHHGHGGGGGGHHBGHGDGGGGGHHHGGGGHHHHHHHHGGGGGHG@DHHGHEGGGGGGBFGGEGGGGGGGG.DFEFFFFFFFDCFFFFFFFFFFFFFFFFFFFFFFFFFFDFDFFFEFFCFF?FDFFFFFFFFAFFFFFFFFFFFBDDFFFFFEFADFFFFFBAFFFA?EFFFBFF
	```

=== "Input"

	```
	tail -n 4 F3D0_S188_L001_R1_001.fastq
	```

=== "Expected Output"


	```
	@M00967:43:000000000-A3JHG:1:2114:11799:28499 1:N:0:188
	TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGATGCCAAGTCAGCGGTAAAAAAGCGGTGCTCAACGCCGTCGAGCCGTTGAAACTGGCGTTCTTGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCCCTACTGACGCTGAGGCACGAAAGCGTGGGTATCGAACAG
	+
	3AAA?AADAFFFCGCGGGFEGCHA?EG?FHHGHGHGGEFHGFHGHF?EFA?EBFGC?EGEFHHHHHH3EEGEEGHFH@E0BCA/CGFHHHDGGGFFF/@DGGDGFHHHHBGH.<<AGGHHHHGHEGE?-ABGF;FFGGDGGGGGGG.CCFEFFF/9;9BFFFFFFFFFFFFFFFFFFFFFFFFFFBDFFFFFFFFCBAF9.AFF/FFAAFFADAFFEFFFFFBDDFFFF.DFFFFFFDDFA;BFFDEFFFF
```




**FASTQ** files should not be confused with **FASTA** files. FASTQ files contain information about the quality of the sequence, but FASTA files only contain the sequence and an identifier.

::: info
	
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


Let's look at a synthetic FASTA file. Because each entry of a .fasta file consists of 2 lines, let's modify head and tail to look at the first and last two lines of HMP_MOCK.v35.fasta.

=== "Input"

	```
	head -n 2 HMP_MOCK.v35.fasta 
	```
=== "Expected Output"

	```
	> A.baumannii.1 
	TGGGGAATATTGGACAATGGGGGGAACCCTGATCCAGCCATGCCGCGTGTGTGAAGAAGGCCTTATGGTTGTAAAGCACTTTAAGCGAGGAGGAGGCTACTTTAGTTAATACCTAGAGATAGTGGACGTTACTCGCAGAATAAGCACCGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGGTGCGAGCGTTAATCGGATTTACTGGGCGTAAAGCGTGCGTAGGCGGCTTATTAAGTCGGATGTGAAATCCCCGAGCTTAACTTGGGAATTGCATTCGATACTGGTGAGCTAGAGTATGGGAGAGGATGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGATGGCGAAGGCAGCCATCTGGCCTAATACTGACGCTGAGGTACGAAAGCATGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCATGCCGTAAACGATGTCTACTAGCCGTTGGGGCCTTTGAGGCTTTAGTGGCGCAGCTAACGCGATAAGTAGACCGCCTGGGGAGTACGGTC
	```

=== "Input"

	```
	tail -n 2 HMP_MOCK.v35.fasta 
	```
	
=== "Expected Output"

	```
	>S.pneumoniae.1
	TAGGGAATCTTCGGCAATGGACGGAAGTCTGACCGAGCAACGCCGCGTGAGTGAAGAAGGTTTTCGGATCGTAAAGCTCTGTTGTAAGAGAAGAACGAGTGTGAGAGTGGAAAGTTCACACTGTGACGGTATCTTACCAGAAAGGGACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTCCCGAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGCAGGCGGTTAGATAAGTCTGAAGTTAAAGGCTGTGGCTTAACCATAGTAGGCTTTGGAAACTGTTTAACTTGAGTGCAAGAGGGGAGAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAGGAACACCGGTGGCGAAAGCGGCTCTCTGGCTTGTAACTGACGCTGAGGCTCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCTGTAAACGATGAGTGCTAGGTGTTAGACCCTTTCCGGGGTTTAGTGCCGTAGCTAACGCATTAAGCACTCCGCCTGGGGAGTACGACC
	```

=== "Exercise"

	Sometimes you know a file or directory exists, but you can't find it. Sometimes you want to find many files with similar properties. This is where the wildcard (`*`) comes in handy.  What do the following commands do?
	 
	1. `ls *` 
	1. `ls F3D*`
	1. `ls *fasta`
 

=== "Answer"

	1. `ls *` lists files in the working directory and 1 level down. 
	1. `ls MiSeq/F3D*` lists files in the data/MiSeq directory that start with "F3D".
	1. `ls MiSeq/*fasta` lists files in the data/MiSeq directory that end with "fasta".


A lot of the time we want to know if a file contains what we expect. A useful thing to do is to be able to **search the contents of files** for a particular string of characters you would like to find.  We can use the file pattern searcher `grep` to find things.

The `MiSeq/` directory contains many of the sequence files ending in`.fastq`. We expect these files to contain information in a particular format throughout the file with four lines of information for each sequence string. Looking through a million line file using less will take a long time. Rather than manually looking at the whole file, we can print only a portion of the file's contents to standard output. 

Let's imagine you would like to find the sequence `CATTAG` in your MiSeq files. We can also use the **wildcard** regular expression to search `CATTAG` in all of the .fastq files located in our current working directory:

```
grep CATTAG F3D0_S188_L001_R2_001.fastq
grep CATTAG *.fastq
```

=== "Exercise"

	What line does `CATTAG` occur on in `F3D141_S207_L001_R1_001.fastq`? 

=== "Answer"

	Use `grep --help` to search for `grep` options related to line number. 
	`grep -n [filename]` will print the line number.



In addition to searching for nucleotide sequences, you may want to search for information in the first line of a .fastq or .fasta file. The `^` (shift + 6) can be used to specify "the beginning of the line".


=== "Input"

	```
	grep "^>" *fasta
	```

=== "Expected Output"

	This will print the name associated with a given sequence in the searched files. In this case, there is only one fasta file, so the name is not printed. 

	```
	>A.baumannii.1
	>A.odontolyticus.1
	>B.cereus.1
	...
	>S.agalactiae.1
	>S.mutans.1
	>S.pneumoniae.1
```


We can also print the line before or after the line that matches a pattern with `-B 1` `-A 1`, respectively.

=== "Input"

	```
	grep -A 1 "^>" *fasta
	
	```


=== "Expected Output"

	This will print the name and sequence for every entry. The first is shown here.

	```
	>A.baumannii.1
	TGGGGAATATTGGACAATGGGGGGAACCCTGATCCAGCCATGCCGCGTGTGTGAAGAAGGCCTTATGGTTGTAAAGCACTTTAAGCGAGGAGGAGGCTACTTTAGTTAATACCTAGAGATAGTGGACGTTACTCGCAGAATAAGCACCGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGGTGCGAGCGTTAATCGGATTTACTGGGCGTAAAGCGTGCGTAGGCGGCTTATTAAGTCGGATGTGAAATCCCCGAGCTTAACTTGGGAATTGCATTCGATACTGGTGAGCTAGAGTATGGGAGAGGATGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGATGGCGAAGGCAGCCATCTGGCCTAATACTGACGCTGAGGTACGAAAGCATGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCATGCCGTAAACGATGTCTACTAGCCGTTGGGGCCTTTGAGGCTTTAGTGGCGCAGCTAACGCGATAAGTAGACCGCCTGGGGAGTACGGTC
	```



As you have seen, `grep` is very useful for finding things within files, and the `*` or wildcard is useful for listings files that match a partial pattern. But, how do we find files when we don't know their location? The `find` command.

Let's navigate back to our home directory and use `find` to command to look for .fasta files. Use the`-name` flag to specify that you are looking for a file with the name listed in double quotes. Use the `*` wildcard to only search for files with a specific extension.

=== "Input"

	```
	find ~ -name "*.fasta"
	```
	
=== "Expected Output"	

	This reveals that .fasta file was found in both the books and the MiSeq directory.  
	
	```
	/home/jovyan/books/yeast.fasta
	/home/jovyan/MiSeq/HMP_MOCK.v35.fasta
	```


=== "Exercise"

	1. Which directories contain a `README.md` file? 
	2. Which directories contain images?

=== "Input"

	Use the commands:
	
	```
	find . -name "README.md"
	find . -name "*.png"
	```
	
=== "Expected Output"	

	to show the following README.md files 

	```
	./seattle/README.md
	./books/README.md
	./MiSeq/README.md
	./README.md
	./southpark/README.md
	```
	
	and the following images.
	
	```
	./images/rstudio-binder-setup.png
	./images/MiSeq-readcount-Mothur.png
	./rstudio-terminal.png
	./CFDE-logo.png
	```

!!! note "Key Points"

	|Command|Description|
	|-|-
	| `grep [option] [filename]`  | selects lines in files that match patterns|
	| `wc` [filename] | prints the total characters, words, and lines in a file
	|`find [path] [conditions]` | finds files with specific properties that match patterns|
