# Upload Data to Cavatica

There are several ways in which users can upload data from their local computers or academic clusters to Cavatica. Cavatica provides many [tutorials](https://docs.cavatica.org/docs/upload-your-data-to-cavatica) on how to do so.

In this tutorial, we expand on using the Cavatica's Command Line Uploader to move fastq files from your AWS Instance on to Cavatica via the command line interface.


!!! note "Learning Objectives"

    - Learn to upload files to Cavatica

    - Learn to edit metadata of files on Cavatica

=== "Est. Time"

    ~ 30 min

=== "Est. Cost"

      < $1.00

=== "Prerequisites"

    - AWS account

    - Cavatica account (check out all requirements in the [Register for Cavatica](Portal-Setup-And-Permissions/KF_4_Cavatica_Registration.md) page)

    - Basic command line


## Step 1: Update Instance

Visit the [AWS tutorial webpage](../Introduction_to_Amazon_Web_Services/introtoaws2.md) to launch a 64 bit `Ubuntu Server 20.04 LTS (HVM), SSD Volume Type` instance. `LTS 20.04` is frozen at `version 20.04`, and so you always needs to be updated the app installer before you start:

=== "AWS Instance Code"

    ```
    sudo apt update
    ```

## Step 2: Download example data

Next, make a directory called "fastq" using the command mkdir, and then download some example fastq data that you can move from AWS to Cavatica:

=== "AWS Instance Code"

    ```
    mkdir fastq
    cd fastq
    curl -LO https://osf.io/5daup/download ERR458493.fastq.gz
    curl -LO https://osf.io/8rvh5/download ERR458494.fastq.gz
    ```

[Curl](https://curl.se/) is an open source software that transfers data. The -L flag redirects the user to the right URL if the server reports that the requested page has been moved. The -O flag names the local file the same as its remote counterpart. These example files are from Schurch et al, 2016 yeast RNAseq study. The exact nature of the files does not matter for this tutorial. Any file type may be used instead of the fastq.

## Step 3: Install Java and Download Command Line Uploader

Cavatica Command Line Uploader needs `java version "1.8.0_20"`. Ubuntu Server 20.04 LTS (HVM) does not come with java pre-installed. You will need to install it with this command:

=== "AWS Instance Code"

    ```
    cd
    sudo apt install openjdk-8-jre-headless
    ```
The command `cd` takes you to the home directory.

Next, download the Cavatica Uploader by running this code:

=== "AWS Instance Code"

    ```
    curl -LO https://cavatica.sbgenomics.com/downloads/cli/cavatica-uploader.tgz
    ```

Now uncompressed the Cavatica Uploader by running:

=== "AWS Instance Code"

    ```
    tar zxvf cavatica-uploader.tgz -C ~
    ```

## Step 4: Test the Command Line Uploader

Check if the Uploader works by running this code:

=== "AWS Instance Code"

    ```
    ~/cavatica-uploader/bin/cavatica-uploader.sh -h
    ```

=== "Expected Output"

    ```
    ubuntu@ip-172-31-26-145:~$ ~/cavatica-uploader/bin/cavatica-uploader.sh -h
    Upload files to Cavatica
    usage: cavatica-uploader.sh [-h] [-l] [-p id] [-t token] [-x url] file ...
    -a,--automation                 Start automation from manifest file.
                                 This option must be used together with
                                 --manifest-file.
    --dry-run                    Dry run the upload (manifest) and/or
                                 metadata setting process.
    -f,--folder <arg>            Specify optional folder, inside of a
                                 specified project, to upload the files
                                 into.
                                 You can specify nested folder structure
                                 separated by the path separator `/`.
                                 If any of the specified folders is
                                 missing it will be created.
    -h,--help                    Print a short usage summary.
    -l,--list-projects           Print a list of projects available as
                                 upload targets. The output is a
                                 tab-separated list of project IDs and
                                 names.
    --list-tags                  Print a list of tags present in a project
                                 and exit.
                                 This option must be used together with
                                 --project.
    -mf,--manifest-file <arg>       Specify manifest tabular file to set
                                 metadata.
                                 This option must be used together with
                                 --project.
    -mm,--manifest-metadata <arg>   Parse metadata from manifest file.
                                 You can list metadata field names as
                                 argument to this option.
                                 This option must be used together with
                                 --manifest-file.
    -p,--project <arg>              Specify the ID of the project to upload
                                 files to.
    -pf,--preserve-folders          Should the folder structure for specified
                                 input folders be preserved while
                                 uploading recursively.
                                 By default, files encountered in the
                                 nested folders are `flattened`, and
                                 uploaded into root target folder.
                                 With this flag, inner folders are created
                                 along the way, and files are uploaded
                                 into them.
    --skip-partial               Do not attempt to resume incomplete
                                 uploads.
                                 If omitted, the uploader will resume an
                                 upload when the local file matches in
                                 name and size.
    -t,--token <arg>                Specify an authorization token.
    --tag <arg>                  Apply tag <arg> to all the files in this
                                 upload.
                                 This option may appear multiple times.
    -u,--username <arg>             Specify username.
                                 If omitted and not using the -t option,
                                 user will be prompted for a username.
    -x,--proxy <arg>                Specify a proxy server through which the
                                 uploader should connect.
                                 The URL to the proxy server in the form
                                 proto://[user:pass]@host[:port].
                                 Proto can be `http' or `socks'. Supports
                                 SOCKS4 and SOCKS5.
    The program outputs a tab-separated list of newly created remote file IDs
    and local file names.
    Complete documentation is available at:
    http://docs.sevenbridges.com/docs/upload-via-the-command-line

    ```

## Step 5: Move files

### Find your Cavatica Authentication Token and Username

The Authentication Token is a personalized code in Cavatica that allows other programs to get access to your Cavatica account. You can [find the Cavatica Authentication](Portal-Setup-And-Permissions/KF_5_ConnectingAccounts.md) token on you Cavatica account under the "Developer" tab.

Copy the Authentication token. You will replace `a???????????????????????????????` in the code block below with your own token.

Next, find and remember your username at the top right of your Cavatica account page. Replace `username` in the code block below with your personal username.


### Choose a Cavatica Project

You can either create a new project or choose an existing project.

#### New project
Create a new Cavatica project by clicking on the "Projects" tab on the Cavatica homepage and selecting the " + Create a project" option. You can name your new project whatever you like.  Remember your project name and replace "project-name" with the name of your project.


#### Existing project
Alternatively, you may choose to select an existing project. To get a list of all existing projects to choose from, run this code:

=== "AWS Instance Code"
  ```
  ~/cavatica-uploader/bin/cavatica-uploader.sh -t a??????????????????????????????? --list-projects
  ```
  Please replace `a???????????????????????????????` with your own Authentication token.

Remember your project name and replace "project-name" with the name of your project.

!!! Important
    If you have `_`s in your project name, replace them with `-` in the uploader code.


### Moving Files

Finally, you can transfer files by running this code:

=== "AWS Instance Code"
    ```
    cd ~/fastq
    ~/cavatica-uploader/bin/cavatica-uploader.sh -t a??????????????????????????????? -p username/project-name *.fastq.gz
    ```
=== "Expected Output"

    ```
    Initializing upload...
    Starting upload of 2 item(s) to 'project-name'
    Uploading file '/home/ubuntu/fastq/ERR458493.fastq.gz'
    5fadb1cee4b05495de67d7ea	/home/ubuntu/fastq/ERR458493.fastq.gz    100%
    Uploading file '/home/ubuntu/fastq/ERR458494.fastq.gz'
    5fadb1dce4b05495de67d7f1	/home/ubuntu/fastq/ERR458494.fastq.gz   45.52%
    ```

By running this code, you are moving files from your AWS instance to the "project-name" project in Cavatica. The `*` wildcard copies all the files ending in ".fastq.gz". The `-t` flag tells AWS to look for an Authentication token, the `-p` flag tells AWS which Cavatica project to put the files into. If you wish, you could use the `-f` flag to create a subfolder inside the project (not shown here).

You're all done! Log in to Cavatica and look for your files in the the "Files" tab in your project.
