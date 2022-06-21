---
title: Appendix
---

All the code from this lesson in one chunk.

```bash
PS1='$ '
clear

pwd

ls

ls -l

ls --help

cd books/
pwd
ls -lh

cd ~/books/
head README.md

head *.txt

gunzip WizardOfOz.txt.gz

cp book.txt book-copy.txt
ls

mv book-copy.txt book-2cities.txt

mkdir data results images/
mkdir -p data/results/images

rmdir data/results/images
rmdir data/results

cd ~/MiSeq
ls

head -n 4 F3D0_S188_L001_R1_001.fastq

tail -4 F3D0_S188_L001_R1_001.fastq

head -4 HMP_MOCK.v35.fasta 

cd ../data/MiSeq/
grep CATTAG F3D0_S188_L001_R2_001.fastq
grep CATTAG *.fastq

wc F3D0_S188_L001_R2_001.fastq

wc -l F3D0_S188_L001_R2_001.fastq

mkdir results
grep CATTAG *.fastq > results/files-with-CATTAG.txt

grep CATTAG *.fastq 
grep CATTAG *.fastq | wc -l

for file in *fastq
do
grep CATTAG $file | wc -l
done

for file in *fastq
do
echo $file
grep CATTAG $file | wc -l
done

history

history > history.txt
```

