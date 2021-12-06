# Building pipelines

By default, the output of many UNIX commands is sent to the
standard output, or "stdout". (Some commands also produce a standard error, or "stderr", which is where errors are printed; and standard input, or "stdin", which is where input comes from.)

Much of the power of the UNIX command line comes from working with
stdout output, and if you work with UNIX a lot, you'll see characters
like the `|` (pipe), `>` (redirect), `>>` (append)  thrown around. These
are redirection commands that say, respectively, "send stdout to a new
file", "append stdout to an existing file", and "send stdout from one
program to another program's stdin."

If you know you want to save an output file, you can use the redirect symbol `>`. Note, if you want to save a file in a different directory, that directory must exist.

Let's practice building pipelines using different datasets

=== "books"

```
cd ~/books
wc -l *txt
wc -l *txt | sort
wc -l *txt | sort -nr
wc -l *txt | sort -nr
mkdir results
wc -l *txt | sort -nr > results/booklengths.txt

```



=== "MiSeq"


```
cd ~/MiSesq
for file in *fastq
do
echo $file
grep CATTAG $file | wc -l
done
```


=== "southpark"


````
cd ~/southpartk
for character in Kenny Cartman Chef Kyle
do
echo $character
cut -d, -f3 All-seasons.csv | grep $character | wc -l
done
````
