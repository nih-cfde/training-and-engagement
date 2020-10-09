---
layout: page
title: Data Download via Cavatica
---

Download KF data from the Cavatica portal
================================================

The easiest way to download a large number of KF files is through [Cavatica](https://cavatica.sbgenomics.com/).

## Step 1: Download the list of links

* On your local machine, create a folder on the Desktop called `KF_Data`.

* Open up the [KF portal](https://portal.kidsfirstdrc.org/dashboard) on a web browser and [choose your files of interest](../KF_7_PushToCavatica.md).

* Then [push those files to Cavatica](../KF_7_PushToCavatica.md). If you have multiple Cavatica projects, be sure to remember which project you chose to transfer your files into.

* Switch to the [Cavatica interface](https://cavatica.sbgenomics.com/) and click on the project name.

* Click the "Files" tab on the top left to view the list of files you just pushed from the KF portal.

* Select all files by checking the box next to "Name". A tiny pop-up window tells you how many files you have selected.

![](../../images/Data_Download_selecting_files_cavatica.png "Check boxes to select files")

* Click on the **drop down arrow** *next to* the "Download" button (shown in the figure below)

!!! Warning
    Do not click on the "Download" button

    Clicking on the "Download" button instead of the arrow next to it will result in the download of individual files (one at a time). While this works for a small number of file, it is unfeasible for a large number of files. If you selected 50 files, 50 pop-up windows will appears for manual download of individual files. If this happens on accident, you can close one pop-up window at a time.


![](../../images/Data_Download_Cavatica_Download_Links.png "Click Download arrow for download links")


* A pop-up window will appear after a few seconds.

![](../../images/Data_Download_Cavatica_Download_Links_2.png "Download links window")


* Click "Download link" and save the text file in the "KF_Data" folder on the "Desktop" of your local machine.

* A text file called "download-links.txt" will appear on your computer. It contains a list of links that you can use for downloading vcf files.

## Step 2: Download data with wget

!!! Important
    The download links expire after 48 hours. If you need to re-download files after 48 hours of first download, you must go back to Cavatica and download the "dowload-link.txt" file once again.

* Open up your Mac or Linux terminal window and navigate to the `KF_Data` folder by typing:

```
cd ~/Desktop/KF_Data
```

* Run the following `wget` command to download all the files:

```
wget -i download-links.txt
```
The `-i` flag tells `wget` to look for a list of URLs from a file.

* Depending on the number of files you selected, the download command may take a few seconds to minutes. When it is complete, all your files will appear in your "KF_Data" folder.

### Terminal vidlet: wget to download data

<script id="asciicast-355576" src="https://asciinema.org/a/355576.js" async></script>

## Step 3: Clean up file names

Notice that the file names look odd? There's a lot of text after the `.vcf.gz` part:

`cb1676e6-0f56-4972-92b4-f343a83dfd3e.strelka.PASS.vep.vcf.gz?X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=AKIAIZAHWHQWDB7XV4RA%2F20200825%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Date=20200825T150049Z&X-Amz-Expires=3600&X-Amz-SignedHe`

Not only does this extra text clutter your filenames, it also prevents your computer from recognizing the file format.

One way to clean up file names is to change each file name manually. However, manual editing takes too much time if you are working with multiple files. Luckily, you can run one line of code to clean up all file names at once:

```
for i in `ls`; do mv $i `ls $i | grep 'vcf.gz' | cut -d'?' -f1`; done
```

This command loops through the list of file names in "KF_Data" and renames the files by deleting all characters after the `vcf.gz` extension starting from the `?` in the file name (`-d` specifies a text delimiter, `-f1` keeps the first element which is the file name up to the `vcf.gz` part).

!!! Important
    This clean up command was written to look through files with the pattern `vcf.gz` in them. To apply the same renaming command to files with other extensions, such as a `.tsv` file, modify the code as follows:

    ```
    for i in `ls`; do mv $i `ls $i | grep 'tsv.gz' | cut -d'?' -f1`; done
    ```

Check your folder - you have now successfully downloaded KF files to your local machine!

### Terminal vidlet: clean up file names

<script id="asciicast-355593" src="https://asciinema.org/a/355593.js" async></script>
