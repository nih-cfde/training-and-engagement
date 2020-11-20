# Uploading and downloading files

4. Set up gcloud authorization so you can move files to/from your instance.

- `gcloud init`
- type "2"
- click the link that appears on the terminal. A new web browser page will open, log in with your GCP google account
- copy/paste the verification code back to the terminal
- enter the number that corresponds to your project
- you can configure a default Compute Region and Zone or not



5. This is how to move files to instance.

- for example, this command downloads the book demo to instance (there are a lot of files, you can cancel the download with Ctrl+z).
```
# -r is a recursive flag to download all files in the v1 directory
# -m parallel multi-threaded/multi-processing copy, use when copying many large files
gsutil -m cp -r gs://genomics-in-the-cloud/v1/* ~/book/
```


**We need data to work with - follow these steps to upload files to a Google storage bucket, either to analyze on a custom VM on the GCP interface or in Terra.**

*This might be useful for users bringing their own data (e.g., broader use of Terra/GCP besides CFDE purposes) or loading test data.*

1. If you want to move local (on laptop) data to the cloud, you need a local installation of `gsutil`. Instructions for installing on a Mac: https://cloud.google.com/storage/docs/gsutil_install

FYI: To upload or download data with the command line, users need to install `gsutil`. After installation, users need to initiate `gcloud` with `gcloud init`, however this command will not work unless users choose to put gcloud in their system’s PATH. Otherwise, the `gcloud` installation lives in the directory `/bin/gcloud` where google-cloud-sdk was installed. After running the gcloud init command, users need to log in to their GCP account and authorize access. According to the instructions, the next step is to run `gsutil` commands, however I get an error message: “ServiceException: 401 Anonymous caller does not have `storage.objects.get` access to the Google Cloud Storage object.” Again, the problem is that I didn’t set `gcloud` in my PATH, so I need to run `./google-cloud-sdk/bin/gsutil cp` to move files. Now it works!

2. Either create a new storage bucket or
```
# create google bucket.
gsutil mb gs://my-bucket

# create bucket alias
export BUCKET="gs://my-bucket"

# copy local files to the bucket, can use -m for multi-thread and -r for recursive copying from a directory
gsutil cp [files] $BUCKET
```

3. Use a storage bucket already associated with a Terra workspace
```
# if you want to access these files on Terra
# create a Terra workspace, copy the bucket name that is
# auto-assigned to the workspace
# for example, i have a workspace with this bucket name
export mybucket="gs://fc-a2cdb170-5198-4a53-a631-86d22564d199"
gsutil cp test_moving_this.txt $mybucket

# if you didn't set the PATH to find gcloud, the command will look like this
./google-cloud-sdk/bin/gsutil cp test_moving_this.txt $mybucket
```

After the files are in the bucket, they are accessible on Terra! Check the Data tab > Files page. You should see your files! To use the files in an analysis on Terra, you will then need to format the sample table section.

Globus can also be used to access files on Google Drive or Google Storage; have to set up the end points to point to Google - read more here: https://docs.globus.org/how-to/access-google-storage/
