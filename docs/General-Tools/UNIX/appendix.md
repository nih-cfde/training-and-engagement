---
title: Appendix
---

All the code from this lesson in one page.

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

grep CATTAG F3D0_S188_L001_R2_001.fastq
grep CATTAG *.fastq
grep "^>" *.fasta
grep -A 1 "^>" *.fasta

cd ~
find . -name "*.fasta"

cd ~/MiSeq
wc -l *.fastq

head -n 1 *.fatsq
grep "^@M00967" *R1*.fastq | wc -l
grep "^@M00967" *R1*.fastq  | head

grep "^@M00967" F3D0_S188_L001_R1_001.fastq | wc -l
grep "^@M00967" F3D0_S188_L001_R1_001.fastq | wc -l
grep "^@M00967" F3D142_S208_L001_R1_001.fastq | wc -l

for file in *R1*.fastq
do
echo $file
grep "^@M00967" $file | wc -l
done


history

history > ~/history.txt
```

