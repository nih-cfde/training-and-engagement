---
title: Appendix
---

All the code from this lesson in one chunk.

```bash
clear
PS1='$ '
clear

ls
ls --help
ls -l books
ls -l MiSeq
ls -lh MiSeq
ls -lh *

cd books
pwd
ls
cd ../Miseq
pwd
ls
ls *fastq
ls F3*

cd ~books

head README.md

head book.txt
tail book.txt
cat book.txt
less book.txt

head *txt

wc -l *txt
wc -l *txt | sort
wc -l *txt | sort -nr

gunzip TheWonderfulWizardofOz.txt

rm book.txt
mkdir results

wc -l *txt | sort -nr > results/book_lengths.txt
cat results/book_lengths.txt

grep "Chapter" *txt
grep -i "Chapter" *txt
grep "CHAPTER" *txt > results/chapter_titles.txt
cat results/chapter_titles.txt

grep "The" *txt
grep "^The" *txt
grep -w "^The" *txt
grep -w -A 1 "^The" *txt
grep -w  "^The" *txt > results/The.txt
cat results/The.txt

head README.md
head -4 F3D0_S188_L001_R1_001.fastq
tail -4 F3D0_S188_L001_R1_001.fastq
less F3D0_S188_L001_R1_001.fastq

wc README.md
wc F3D0_S188_L001_R1_001.fastq
wc -l F3D0_S188_L001_R1_001.fastq	
wc -l *fastq
wc -l *R1*fastq
wc -l *R1*fastq | sort -nr

head -4 F3D0_S188_L001_R1_001.fastq
grep "^@M" F3D0_S188_L001_R1_001.fastq
grep "^@M" F3D0_S188_L001_R1_001.fastq | wc -l
grep "^@M" *R1*.fastq | wc -l

for file in *R1*.fastq
do
echo $file
grep "^@M" $file | wc -l
done

mkdir results


for file in *R1*.fastq
do
echo $file >> results/read_count.txt
grep "^@M" $file | wc -l >> results/read_count.txt
done

head results/read_count.txt

for file in *R1*.fastq
do
echo $file >> results/samples.csv
grep "^@M" $file | wc -l >> results/count.csv
paste -d , results/samples.csv results/count.csv > results/read_count.csv
done


```

