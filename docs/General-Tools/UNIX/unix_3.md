# Reading large files

Now that we know what files exist on our computer, it is time to look at the contents of the file. There are multiple ways to look at the contents of a file. 

`head` prints, by default the first 10 lines of a file, and `tail` prints the last 10 lines of a file by default. 


The following commands use a similar syntax:

```
head [filename]
tail [filename]
cat [filename]
less [filename]
```


Let's navigate to the `books` directory. 


```
cd ~/books/
```


Let's look at the first 10 lines with `head`. 

=== "Input"

	```
	head book.txt
	```

=== "Expected Output"

	The Project Gutenberg eBook of A Tale of Two Cities, by Charles Dickens

	This eBook is for the use of anyone anywhere in the United States and
	most other parts of the world at no cost and with almost no restrictions
	whatsoever. You may copy it, give it away or re-use it under the terms
	of the Project Gutenberg License included with this eBook or online at
	www.gutenberg.org. If you are not located in the United States, you
	will have to check the laws of the country where you are located before
	using this eBook.

Let's look at the last 10 lines with `tail`. 


=== "Input"

	```
	tail book.txt
	```

=== "Expected Output"

	Most people start at our website which has the main PG search
	facility: www.gutenberg.org
	
	This website includes information about Project Gutenberg-tm,
	including how to make donations to the Project Gutenberg Literary
	Archive Foundation, how to help produce our new eBooks, and how to
	subscribe to our email newsletter to hear about new eBooks.
	
		
The `cat` command prints the entirety of a file to the stdout of our 	computer, which can be overwhelming. The `less` command provides a safe way of looking at the contents of a file without the ability to change it. Use the up and down arrows to scroll through the file. Type `q` to exit the lesson program
	

```
cat book.txt
less book.txt
```

!!! info

	You can use TAB to do filename completion, so if you type `cat R` and then press your Tab key once, it will autocomplete if there is a unique match. If there is more than one match, the first Tab will do nothing, and the second will show all the possible matches.


We can see there are a lot more books, and we can look at the first line of all the .txt files with the *. 


=== "Input"

	```
	head -n 1 *.txt
	```


=== "Expected Output"


	==> Alice_in_wonderland.txt <==
	﻿The Project Gutenberg eBook of Alice’s Adventures in Wonderland, by Lewis 	Carroll
	
	==> A-tale-of-two-cities.txt <==
	﻿The Project Gutenberg eBook of A Tale of Two Cities, by Charles Dickens
	
	==> book.txt <==
	﻿The Project Gutenberg eBook of A Tale of Two Cities, by Charles Dickens
	
	==> PeterPan.txt <==
	﻿The Project Gutenberg eBook of Peter Pan, by James M. Barrie
	
	
	
	
Notice that there is one book that is compressed. We can uncompress it with the command `gunzip`. The `-k` option will keep the original file. 

```
gunzip -k WizardOfOz.txt.gz
```

Now we can return the previous command and also see the first lines of The Wizard of Oz. 


=== "Input"

	```
	head -n1 *.txt
	```

=== "Expected Output"

	==> Alice_in_wonderland.txt <==
	﻿The Project Gutenberg eBook of Alice’s Adventures in Wonderland, by Lewis 	Carroll
	
	==> A-tale-of-two-cities.txt <==
	﻿The Project Gutenberg eBook of A Tale of Two Cities, by Charles Dickens
	
	==> book.txt <==
	﻿The Project Gutenberg eBook of A Tale of Two Cities, by Charles Dickens
	
	==> PeterPan.txt <==
	﻿The Project Gutenberg eBook of Peter Pan, by James M. Barrie
	
	==> WizardOfOz.txt <==
﻿	The Project Gutenberg eBook of The Wonderful Wizard of Oz, by L. Frank Baum






!!! note "Key Points"

	| Command [OPTION] | Description |
	| -------- | -------- | 
	|`head [filename]` | print first 10 lines of  `FILENAME` | 
	|`cat [filename]`| print `FILENAME`'s contents to stdout|
	|`less [filename]`| view `FILENAME` without printing  to stdout |
	| `gunzip [filename]` | uncompress filename
	| `gzip -k [filename]` | compress a file and keep the original |
