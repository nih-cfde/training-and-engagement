# Demo: Custom Workflows

In this section, we will demonstrate how to build a custom workflow for analysis on Terra. For this demo, we will show you how we ran the same workflow as in the [Genome-wide Association Study (GWAS) lesson](../GWAS-in-the-cloud/index.md), but on Terra instead of AWS.

To follow along, please watch the vidlets, which show the components of this process:

- Part 1: Create docker containers with specific software environments for steps of the workflow
- Part 2: Write a workflow script to specify the analysis steps
- Part 3: Test the workflow before uploading it to Terra
- Part 4: Run workflow on Terra

!!! important

    It's beyond the scope of this lesson to teach general docker and Workflow Description Language (WDL) syntax, however, check out the resources linked throughout the lesson.

    We hope this demo helps to show the general steps and get you started with your own workflows on Terra!

<!-- undecided whether to expect ppl to be able to follow all WDL testing steps because it is kind of OS dependent; could run it in a GCP VM, but that requires additional set up. i am considering making this entirely a demo, except perhaps the parts directly on Terra.. -->

## Part 1: Setting up Docker

- what is docker
- why containers are great for creating software environments to share


The core steps of building, testing, and sharing the container are shown in the vidlet below:

<iframe id="kaltura_player" src="https://cdnapisec.kaltura.com/p/1770401/sp/177040100/embedIframeJs/uiconf_id/29032722/partner_id/1770401?iframeembed=true&playerId=kaltura_player&entry_id=1_5pvutxvq&flashvars[localizationCode]=en&amp;flashvars[leadWithHTML5]=true&amp;flashvars[sideBarContainer.plugin]=true&amp;flashvars[sideBarContainer.position]=left&amp;flashvars[sideBarContainer.clickToClose]=true&amp;flashvars[chapters.plugin]=true&amp;flashvars[chapters.layout]=vertical&amp;flashvars[chapters.thumbnailRotator]=false&amp;flashvars[streamSelector.plugin]=true&amp;flashvars[EmbedPlayer.SpinnerTarget]=videoHolder&amp;flashvars[dualScreen.plugin]=true&amp;flashvars[Kaltura.addCrossoriginToIframe]=true&amp;&wid=1_cww3jbnf" width="608" height="402" allowfullscreen webkitallowfullscreen mozAllowFullScreen allow="autoplay *; fullscreen *; encrypted-media *" sandbox="allow-forms allow-same-origin allow-scripts allow-top-navigation allow-pointer-lock allow-popups allow-modals allow-orientation-lock allow-popups-to-escape-sandbox allow-presentation allow-top-navigation-by-user-activation" frameborder="0" title="Kaltura Player"></iframe>

More information about each step:

=== "1. Install docker"

    - for Mac: https://docs.docker.com/docker-for-mac/install/


=== "2. Dockerhub account"

    In order to share the docker, we are pushing the image to Dockerhub. There are other options (i.e., Dockstore). For the workflows used in Terra, either of these platforms will work!

    Go to <https://hub.docker.com/> to create an account.

=== "3. Dockerfile"

    We tell docker what we want installed in the container with a file called ["Dockerfile"](./Dockerfile.md). We installed dependency software tools, vcftools, and the R package "qqman". There's an existing docker container for plink (`gelog/plink:latest`), so we use that in the workflow for plink steps.

    The base computer image of our container is ubuntu 20.04, and the software installation commands are similar to those used in the [GWAS lesson](../GWAS-in-the-cloud/index.md), with some modifications for docker syntax.

    - `WORKDIR /usr` specifies the working directory - since dockers are closed containers, it helps to know where commands will be run so you can locate output files.
    - we use `ENV DEBIAN_FRONTEND noninteractive` because docker doesn't allow interactive installation.
    - `RUN` precedes all software install commands, and strings of commands are connected by `&&`. The `\` allows commands to be written over multiple lines in the Dockerfile file. After the first `RUN` line, there is an indent for the others in the installation block.

=== "4. Build"

    We built the docker with a specified tag list (`-t`). In this case, the tag list contains the user name (`mlim13`), the docker name (`demo_gwas`), and the version (`tag0`). The last element of the command specifies the location of our "Dockerfile"; since we're in the directory that contains it, we can use the `.` notation.

    ```
    docker build -t mlim13/demo_gwas:tag0 .
    ```

=== "5. Test"

    We set up an interactive (`-it`) session of our docker container with `docker run`:

    ```
    docker run -it mlim13/demo_gwas:tag0
    ```

    And ran a test command to make sure the software `vcftools` was installed:

    ```
    vcftools -h
    ```

=== "6. Share"

    Pushing to Dockerhub takes a few minutes to complete:

    ```
    docker push mlim13/demo_gwas:tag0
    ```

    After pushing the container, it is publicly available on [Dockerhub](https://hub.docker.com/repository/docker/mlim13/demo_gwas).

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

=== "1. Install tools"

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


=== "2. Womtool"

    We used Womtool to test the WDL script syntax. Note that the WDL may still fail to run, but Womtool helps to catch syntax errors.

    ```
    # check WDL
    java -jar womtool-47.jar validate hello.wdl

    # womtool will make the basic input json file!
    # but you will need to customize it
    java -jar womtool-47.jar inputs hello.wdl > hello.json
    ```


=== "3. Run WDL"


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
