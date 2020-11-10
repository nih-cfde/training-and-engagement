# Simulating Illumina Reads

InSilicoSeq is a software that simulates Illumina reads from genomes, and was especially built to simulate reads from metagenomes. InSilicoSeq is written in Python, is fairly well-documented and easily installed via Python’s package manager `pip(3)`.

## Why use simulated data?

There are a myriad of bioinformatics tools currently available to users. Before use, every tool must be tested for general usability and suitability for a given task. Simulated data can be ideal for the purposes of testing bioinformatics tools because the user can model the data set based on their specified parameters. With simulated data, the user knows the final answer, and so they know when the tool doesn't work. Simulated data sets are also great for use in the classroom, especially when protected (human) data can't be distributed.

## Selecting the right AWS instance to run InSilicoSeq

To run InSilicoSeq, you will need to launch a 64 bit Ubuntu Server 20.04 LTS (HVM), SSD Volume Type instance as described in [the AWS tutorial](Introduction_to_Amazon_Web_Services/introtoaws2.md) with a **few modifications**:

 1) You must select the the t2.xlarge instance (instead of the default t2.micro instance selection described in the tutorial):

![t2 xlarge](../images/Simulated_Data_t2xlarge.png)

2) Add additional storage by selecting the "4.Add Storage" tab on the instance launch page and then changing the number on the "Size (GB)" tab to read "16", as shown in the image below:

![Add storage to t2 xlarge](../images/Simulated_Data_t2xlarge_storage.png)

3) Click "Review and launch"

4) Then go back to the [AWS tutorial](Introduction_to_Amazon_Web_Services/introtoaws2.md) and follow instructions on how to access the instance via the MacOS terminal window

## Installing InSilicoSeq

!!! Important

    You must run all these commands and installations on your AWS Ubuntu instance and **not** your local terminal


InSilicoSeq can be installed using Python's package manager `pip3`. However, AWS's Ubuntu 20.04 instance does not come pre-installed with `pip3`. Let's begin by updating the instance and installing `pip3`.

=== "AWS Instance Code"

    ```
    sudo apt update
    sudo apt install python3-pip
    ```

Now you can install InSilicoSeq:

=== "AWS Instance Code"

    ```
    pip3 install InSilicoSeq
    ```

## Download the human reference genome

InSilicoSeq simulates reads based on one or more input reference genomes. You will use the human reference genome GRCh38.

First, make a directory called "reference_genome" using the command `mkdir`. Then download the compressed (".gz" extension) human reference genome inside the new folder and unzip the ".fna" file.

=== "AWS Instance Code"

    ```
    mkdir reference_genome
    cd reference_genome
    curl -LO ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz
    gunzip GRCh38_latest_genomic.fna.gz
    ```

## Set the path variable

For the Ubuntu machine to recognize the InSilicoSeq executable, we need to add the bin (where InSilicoSeq was downloaded) to Ubuntu's search path.

=== "AWS Instance Code"

    ```
    cd
    export PATH=$PATH:~/.local/bin
    ```

Now your Ubuntu machine should be able to find InSilicoSeq.

## Running InSilicoSeq

Make a directory called "fastq" in which the simulated fastq files can be saved.

=== "AWS Instance Code"

    ```
    mkdir fastq
    cd fastq
    ```

Finally, run the code to make your simulated fastq file.

=== "AWS Instance Code"

    ```
    iss generate --genomes ../reference_genome/GRCh38_latest_genomic.fna --model hiseq --n_reads 5M --cpus 4 --output my_sim --compress true
    ```

=== "Expected Output"

    ```
    ubuntu@ip-172-31-26-57:~/fastq$ iss generate --genomes ../reference_genome/GRCh38_latest_genomic.fna --model hiseq --n_reads 5M --output 1_sample
    INFO:iss.app:Starting iss generate
    INFO:iss.app:Using kde ErrorModel
    INFO:iss.util:Stitching input files together
    INFO:iss.app:Using lognormal abundance distribution
    INFO:iss.app:Using 2 cpus for read generation
    INFO:iss.app:Generating 5000000 reads
    INFO:iss.app:Generating reads for record: NC_000001.11
    INFO:iss.app:Generating reads for record: NT_187361.1
    INFO:iss.app:Generating reads for record: NT_187362.1
    INFO:iss.app:Generating reads for record: NT_187363.1
    INFO:iss.app:Generating reads for record: NT_187364.1
    INFO:iss.app:Generating reads for record: NT_187365.1
    ....
    ```

!!! note "What are all these tags?"

    --model lets you input an error model file. The default is None. You're using hiseq for a pre-computed error model suitable for simulating data generated from a HiSeq sequencing run. Other available options are NovaSeq and MiSeq. If you do not wish to use a model, use –mode basic. The name of the built-in models is case insensitive.

    --n_reads specifies the number of reads to be simulated. The default is 1 million (1M).

    --cpus specifies the number of CPUs that should be used. Default is 2.

    --output specifies the name of the output files. You're calling it my_sim, but you can change it to whatever suits your needs.

    --compress true compresses output files to ".gz" format

!!! Important
    It takes approximately 40 mins to run a single simulation on the AWS instance

If your run is successful, you will see two ".fastq" files, "my_sim_1.fastq" and "my_sim_2.fastq", and an abundance file called "my_sim_abundance.txt".

If you want to repeat the simulation multiple times, modify this code block by changing the `n` to an integer corresponding to the number of time you wish to run the simulation. Your output files will be names 1_my_sim, 2_my_sim, ..., n_my_sim.

=== "AWS Instance Code"

    ```
    for i in {1..n};
    do iss generate --genomes ../reference_genome/GRCh38_latest_genomic.fna --model hiseq --n_reads 5M --cpus 4 --compress true --output ${i}_my_sim;
    done
    ```

!!! Warning
    Depending on how many times you run the simulation and/or how many reads you simulate, the code may take a really long time to complete!




## Quality control

With real data, one must always perform some simple quality control to ensure that the raw data looks good. With simulated data, you can do some quality control to make sure nothing went terrible wrong with the simulation itself and that you have the desired read length.

Let's use FastQC to generate a quality control report.

### Install FastQC

To install FastQC run this code:

=== "AWS Instance Code"

    ```
    sudo apt install fastqc
    ```

### Run FastQC
Now run FastQC on your fastq files:

=== "AWS Instance Code"

    ```
    cd fastqc
    fastqc *.gz
    ```
You are using the `*` wildcard to help specify all of the .fastq.gz files here.

### Transfer to Local Computer

The easiest way to visualize the output is to transfer it to your local computer. First, let's list all the files and pull out the report summaries (.html) using ls:

=== "AWS Instance Code"

    ```
    ls *fastqc.zip
    ls *.html
    ```

Now you can move the ".html" files to your local computer to visualize:

    ```
    scp -i ~/Desktop/amazon.pem ubuntu@ec2-??-???-???-??.us-east-2.compute.amazonaws.com:/home/ubuntu/fastq/\*.html ~/Desktop/fastqc/.
    ```

!!! Important

    Remember to replace:

    - `~/Desktop/amazon.pem` with the [path to your amazon.pem](GWAS-in-the-cloud/aws_instance_setup.md)

    - `??-???-???-??.us-east-2.compute.amazonaws.com` with your specific instance. For more details on how to connect to an instance, visit the [AWS set up page of the "GWAS in the Cloud" tutorial](GWAS-in-the-cloud/download_accessAWS.md).


Here is an image of the top of an example FastQC report.

![FastQC report](../images/Simulated_Data_Fastqc.png)



### Interpretation

**Watch this [YouTube video](https://www.youtube.com/watch?v=bz93ReOv87Y) on how to interpret FastQC results**:

[![FastQC video tutorial](http://img.youtube.com/vi/bz93ReOv87Y/0.jpg)](https://www.youtube.com/watch?v=bz93ReOv87Y "FastQC Interpretation")

All simulated reads are 125 bp long :)

## Additional Resources

[Gourlé, H., Karlsson-Lindsjö, O., Hayer, J., & Bongcam-Rudloff, E. (2019). Simulating Illumina metagenomic data with InSilicoSeq. Bioinformatics, 35(3), 521-522](https://academic.oup.com/bioinformatics/article/35/3/521/5055123)

[InSilicoSeq Documentation](https://insilicoseq.readthedocs.io/en/latest/index.html)

[InSilicoSeq GitHub Repo](https://github.com/HadrienG/InSilicoSeq)

[Human reference genome download](https://www.ncbi.nlm.nih.gov/genome/guide/human/)

[Recommended read dept and coverage for NGS applications](https://genohub.com/recommended-sequencing-coverage-by-application/)

[What is FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

[FastQC video tutorial](https://www.youtube.com/watch?v=bz93ReOv87Y)
