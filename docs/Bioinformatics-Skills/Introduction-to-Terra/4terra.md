# Existing workflows in Terra

For this demo, we will show you how to run the same workflow as in the [Genome-wide Association Study (GWAS) lesson](../GWAS-in-the-cloud/index.md), but on Terra instead of AWS.

<work in cfde-teaching-2021/GWAS_demo workspace. check cost on this billing project after run GWAS>


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
