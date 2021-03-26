# Demo: Custom Workflows

In this section, we will demonstrate how to build a custom workflow for analysis on Terra. For this demo, we will show you how we ran the same workflow as in the [Genome-wide Association Study (GWAS) lesson](../GWAS-in-the-cloud/index.md), but on Terra instead of AWS.

To follow along, please watch the vidlets, which show the components of this process:

- Part 1: Create docker containers with specific software environments for steps of the workflow
- Part 2: Write a workflow script to specify the analysis steps
- Part 3: Test the workflow before uploading it to Terra
- Part 4: Run workflow on Terra

!!! note

    It's beyond the scope of this lesson to teach docker and Workflow Description Language (WDL) syntax, however, check out the resources linked throughout the lesson. We hope this demo helps to show the general steps and get you started with your own workflows!

## Part 1: Setting up Docker

- what is docker
- why containers are great for creating software environments to share

<add vidlet here of all the steps below for running docker>


### Step 1: Install docker

- for Mac: https://docs.docker.com/docker-for-mac/install/
    - is there a command line installation option?

### Step 2: Make Dockerhub account
- need to make an account on Docker hub if you want to push docker there
- there's also a google storage option

- https://hub.docker.com/


### Step 3: Write Dockerfile

- make Dockerfile (it must be called "Dockerfile")



### Step 4: Build container


```
# usage: docker build -t <user name>/<docker name>:<version number> <file location of Dockerfile>
docker build -t mlim13/gwas_test:tag0 .
```


### Step 5: Run container

```
docker run -it mlim13/gwas_test:tag0
```

### Step 6: Test container

- add some test commands to demo how commands are run


### Step 7: Push container to Dockerhub

```
# push to dockerhub (username: mlim13, repo: gwas_test, version: tag0)
docker push mlim13/gwas_test:tag0
```
it will be here:
- https://hub.docker.com/repository/docker/mlim13/gwas_test

<!-- ![](./custom-workflow-images/push_to_dockerhub.png) --> <- make this


## Part 2: Write a WDL

Workflow Description Language (WDL, pronounced "widdle") is one of many computer programming languages used to define analysis workflows. It is maintained by an [open online community](https://openwdl.org/) that you can contribute to!

If your WDL is using a docker container, add it to the `runtime` block of a given WDL task. The syntax will look like this:
```
runtime {
  docker: 'mlim13/gwas_test:tag0'
}
```

To test the workflow on local computer, we need 2 files:

- the WDL script, saved with ".wdl" extension: [GWAS_WDL.wdl](./GWAS_WDL.md)
- a JSON file that specifies inputs, saved with ".json" extension: [GWAS_WDL.json](./GWAS_json.md)




## Part 3: Testing a WDL workflow

<add vidlet here of all the steps below for running womtool/cromwell>

### Step 1: Install tools

Install Cromwell and WOMtool:

```
wget https://github.com/broadinstitute/cromwell/releases/download/47/cromwell-47.jar
wget https://github.com/broadinstitute/cromwell/releases/download/47/womtool-47.jar
```

Also need to have Java installed: https://www.oracle.com/java/technologies/javase-jdk15-downloads.html
```
# check installed
java --version
```

!!! warning

    Need to have docker installed and app open otherwise get docker daemon errors.

<alternative is miniwdl, have not tested yet, may add it just as a note for now.>


### Step 2: Test script syntax with Womtool

```
# check WDL
java -jar womtool-47.jar validate hello.wdl

# womtool will make the basic input json file!
# but you will need to customize it
java -jar womtool-47.jar inputs hello.wdl > hello.json
```


### Step 3: Run WDL script


```
# run WDL
# -i = input file
java -jar cromwell-47.jar run hello.wdl -i hello.json
```



#### outputs
- each task gets a call-<task name> directory that has `execution` and `inputs` directories
   - the outputs, logs, stderr files are in `execution`
   - when defining the inputs in WDL script, check how the inputs are partitioned in the `inputs` directory. this was the issue for plink commands and the reason i had to define each input file flag for plink.



## Run custom workflow on Terra

<add vidlet here of all the steps below for running on terra>

We followed these [steps](https://support.terra.bio/hc/en-us/articles/360031366091-Create-edit-and-share-a-new-workflow#h_bc0175df-adb7-422f-b2fe-efab18fd598b) from the Terra documentation. There is also a [video lesson](https://www.youtube.com/watch?v=VtKlYqWBW6A&ab_channel=Terra) from Terra about workflows.

For the demo in the vidlet, these are the specific settings we used:

### Step 1: upload workflow

To add WDL as a new workflow on Firecloud (aka Broad Method Repository):

- namespace = "cfde-workflows"
- name = "test-gwas"
- note that you can redact a snapshot (version) if you don't want people using outdated versions. Select the snapshot version you want to export to the workspace, or update in the workflows tab on Terra.
- set whether it can be shared publicly or not - since this is a test, we have made it private
- note that we do not need to upload the JSON file, the input files are specified in the Terra interface (although there is an option to upload a JSON file if you would like to)

### Step 2: export workflow

To export the new workflow from Firecloud to workspace:

- we used a blank configuration
- set root to participant (this is how inputs are set/selected in the Terra interface)

### Step 3: import data to workspace

The demo data files (".pheno" and ".vcf") are small enough to manually upload (click <span class="highlight_txt">+</span> sign).

!!! tip

    For larger datasets, you can upload files to the Terra workspace Google bucket location with the command line. An example of this process is detailed in the [GCP lesson example](../Introduction-to-GCP/gcp3.md).

    If you want to move local (on laptop) data to the cloud, you need a local installation of `gsutil` ([instructions](https://cloud.google.com/storage/docs/gsutil_install) for installation from GCP). Be sure to add `gcloud` to your system's PATH when asked during set up of the `gsutil` tool!

    Use the storage bucket location already associated with a Terra workspace when you upload the files. For example, we have a workspace with this bucket name, "gs://fc-a2cdb170-5198-4a53-a631-86d22564d199", which we renamed with the alias, "mybucket".

    ```
    # add alias
    export mybucket="gs://fc-a2cdb170-5198-4a53-a631-86d22564d199"

    # then we use gsutil to copy (cp) a file (file.txt) to the bucket
    gsutil cp file.txt $mybucket

    # alternatively, the -r flag recursively copies files from a directory to the bucket
    gsutil cp -r lots_of_files $mybucket
    ```

    After the files are in the bucket, they are accessible on Terra! Your files should be in the Data tab's Files page. To use the files in an analysis on Terra, you will then need to format the sample table section (next step below).

### Step 4: Set up data table

We've uploaded data, but we still have to tell Terra how to reference the input files. This will make it easier to select inputs when we set up the workflow.

- copy paste this table to the add table (alternatively, upload a tab-delimited file (TSV) with this information):

```
entity:sample_id	sample	pheno	vcf
coatColor	coatColor	coatColor.pheno	pruned_coatColor_maf_geno.vcf
```

- then edit the file names so they include the google bucket path (otherwise, it's just a string. You'll know it's correct when the file name has a hyperlink and opens the download window). The bucket path is available from the File tab by clicking on the uploaded file name. Copy the part that looks like: "gs://fc-a2cdb170-5198-4a53-a631-86d22564d199/coatColor.pheno".

### Step 5: Run the workflow!

- click row of data table > open with WDL
- make sure input files are correct ("this.vcf", etc.) > save
- click run!
- check the job manager to see if job is running
- Terra currently does not notify you when a job successfully runs or fails, so check the job manager for status updates


You can view information about the virtual machine Terra spun up for your job from the Job Manager. The compute information will look like this:

<!-- ![](./custom-workflow-images/terra-compute-info.png "Terra VM info") -->


!!! note "Key Points"

    - *ADD*
