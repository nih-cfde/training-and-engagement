# Try all of this on Terra!

- before you can do this:
    - need GCP billing set up to make billing project on Terra
    (cannot do a demo of Terra without this unfortunately. but could make a screencast/vidlet of the process...)


1. add WDL as a new workflow on Firecloud (aka Broad Method Repository)
    - follow these steps: https://support.terra.bio/hc/en-us/articles/360031366091-Create-edit-and-share-a-new-workflow#h_bc0175df-adb7-422f-b2fe-efab18fd598b
    - i called the namespace: 'cfde-workflows' and the name: 'test-gwas'
    - you can redact a snapshot (version) if you don't want ppl going to wrong version. make sure you're on the correct snapshot when the workflow is exported to Terra

2. export workflow from Firecloud to workspace
    - follow steps from: https://support.terra.bio/hc/en-us/articles/360031366091-Create-edit-and-share-a-new-workflow#h_bc0175df-adb7-422f-b2fe-efab18fd598b
    - i used a blank configuration
    - set root to participant (don't really understand what this is yet)

3. import data to workspace
    - for this data (.pheno and .vcf) it's small enough to do manual upload (click '+' sign).
    - i have steps on GCP VM hackmd to use gsutil to programmatically upload data to an existing google bucket via command line

4. make data table to tell Terra what data to use
    - copy paste this table to the add table
```
entity:sample_id	sample	pheno	vcf
coatColor	coatColor	coatColor.pheno	pruned_coatColor_maf_geno.vcf
```
    - then edit the file names so they include the google bucket path (otherwise, it's just a string. you'll know it's right when the file name has a hyperlink and opens the download window). you can get the bucket paths from the File tab by clicking on the uploaded file name. just copy the part that looks like: 'gs://fc-a2cdb170-5198-4a53-a631-86d22564d199/coatColor.pheno'

5. try running workflow!
    - click row of data table > open with WDL
    - make sure input files are correct ('this.vcf', etc.) > save
    - click run!
    - check the job manager to see if job is running
    - as far as i can tell, terra doesn't notify (email) when jobs finish/fail. so you just have to go back and check.
