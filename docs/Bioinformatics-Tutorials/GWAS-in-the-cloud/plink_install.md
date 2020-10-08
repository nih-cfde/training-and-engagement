# Install PLINK

[PLINK](http://zzz.bwh.harvard.edu/plink/index.shtml) is a free, open-source whole-genome association analysis toolset. The software is designed flexibly to perform a wide range of basic, large-scale genetic analyses including GWAS.


## Step 1: Download PLINK with wget
Type these commands into your [Ubuntu AWS terminal](./download_accessAWS.md):

=== "AWS Instance"
    ```
    cd /usr/local/bin/
    sudo wget https://zzz.bwh.harvard.edu/plink/dist/plink-1.07-x86_64.zip
    ```
The first line changes your directory to the "bin". The second line of code downloads the software in the "bin" directory.

## Step 2: Unzip

Now, uncompress the newly downloaded file:

=== "AWS Instance"
    ```
    sudo unzip -o plink-1.07-x86_64.zip
    sudo rm plink-1.07-x86_64.zip
    cd plink-1.07-x86_64
    ```
The first line of code does the uncompressing, the second line of code gets rid of the ".zip" version (you don't need it anymore) and the last line of code takes you into the newly created "plink-1.07-x86_64" directory.

## Step 3: Add PLINK to the .bashrc

After cd-ing into the "plink-1.07-x86_64" folder, add plink to the `.bashrc`. By doing so, you are telling the computer that PLINK can be accessed from any directory. Type:

=== "AWS Instance"
    ```
    echo export PATH=$PATH:$(pwd) >> ~/.bashrc
    source ~/.bashrc
    ```

Check to make sure it's installed:

=== "AWS instance"
    ```
    plink -h
    ```
=== "Expected Output"
    ```
    The top portion of the help page:
    @----------------------------------------------------------@
    |        PLINK!       |     v1.07      |   10/Aug/2009     |
    |----------------------------------------------------------|
    |  (C) 2009 Shaun Purcell, GNU General Public License, v2  |
    |----------------------------------------------------------|
    |  For documentation, citation & bug-report instructions:  |
    |        http://pngu.mgh.harvard.edu/purcell/plink/        |
    @----------------------------------------------------------@


    Please visit the PLINK website for a complete list of options

    A few common options are listed here:

    plink --file {fileroot}     Specify .ped and .map files
      --bfile {fileroot}        Specify .bed, .fam and .map

      --out {fileroot}          Specify output root filename  

    ```

PLINK has now been successfully installed!
