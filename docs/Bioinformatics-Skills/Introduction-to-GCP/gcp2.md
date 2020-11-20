# Setting up a GCP instance

## Step 1: Create a project

- Go to the "Manage Resources" page (it's in the IAM & Admin section) and click "Create Project".

**TO DO: add image**
**TO DO: some images have "My First Project" and one has "test-GCP-VM" as project name...let's make those all the same**

- Enter a project name. The name must be unique. Click "Create".

## Step 2: Configure custom VM

- On the GCP console left-hand panel, scroll down the "Compute Engine" section and select "VM instances". It may take a few minutes to load.

![](../../images/gcp_images/gcp_vm.png "VM instances")

- Click "Create". There are several configuration options to set up:

![](../../images/gcp_images/gcp_vmconfig.png "VM configuration options")

### 1. Name your VM

Names must be in lowercase letters, numbers. Use hyphens "-" instead of spaces.

### 2. Choose a Region

Several regions are available from the dropdown menu. It is less important which zone you choose, and the interface automatically selects a zone based on the region you choose. We are using "us-west1 (Oregon)".

!!! note "Machine regions"

    **TO DO: add what region and zones mean for GCP**

    This tool can help to choose the region closest to you (http://www.gcping.com/). You may need to refresh several times.

### 3\. Choose machine type and configuration

The estimated monthly cost for each machine type is shown on the right side panel of the console page. For this tutorial, select Series "N2" and Machine type `n2-standard-2`.

!!! note "Machine types"

    **TO DO: add info about what these machine types mean, how to select, etc.**

    **TO DETERMINE: is this a good machine to use?**

### 4\. Customize boot disk

Click on "Change". The default operating system is Debian, change it to "Ubuntu" and select version "Ubuntu 20.04 LTS". Set the [persistent disk storage](https://cloud.google.com/persistent-disk); the default is 10Gb. Depending on the tasks you will use the VM for, you may need to increase the storage amount (e.g., to 100Gb).

!!! note "Persistent disk storage"

    **TO DO: add info about what this is for and how to choose**

When you're done configuring the VM, click "Create".

## Step 3: Connect to your VM

- Click on the "SSH" dropdown option for your new instance and select "Open in browser window".

![](../../images/gcp_images/gcp_vmssh.png "VM SSH dropdown")

A new window will open with the instance terminal. The gear icon at the top right can be used to customize terminal window settings (e.g., appearance, text size).

!!! note "Google Cloud Shell"

    Alternatively, you can also access the instance from the Google Cloud Shell or a local laptop terminal (will require installing and setting up [`gcloud`](https://cloud.google.com/sdk/docs/install)).

    The GCP console provides a free Google Cloud Shell. This shell environment is useful for small tasks that do *not* require a lot of CPU or memory (as most bioinformatic analyses do). For example, it is a good place to learn how to use the Google shell environment without incurring cost or to access Google Cloud products (e.g., a Google Storage bucket or GCP virtual machine).

    ![](../../images/gcp_images/gcp_gcshell.png "Activate Cloud Shell button")

    Click on the "Activate Cloud Shell" icon. The first time you start the shell, you'll need to agree to the Google Cloud terms of service and privacy policy. After you start the shell, it may take a few minutes to connect. Check the Google support [documentation](https://cloud.google.com/shell/docs/using-cloud-shell) for more information.

- Access the instance, by entering the `gcloud compute ssh` command in the instance terminal window. You'll need your project ID, zone, and instance name.

**TO DO: test this command again, don't remember having to set the zone**

```
gcloud compute ssh --project <PROJECT_ID> --zone <ZONE> <INSTANCE_NAME>
```

## Step 4: Use the VM

You now have a custom configured VM! Things we could do at this point include:

- [Upload or download files](./gcp3.md)
- [Pull in a docker container](./gcp4.md) to get software without a bunch of software installations
- Run analysis. For example, we could run the [command-line Blast tutorial commands](../Command-Line-BLAST/BLAST1.md) in this GCP instance (instead of the AWS instance used in the tutorial).

**TO DO: test if blast steps work in GCP instance - probably can choose a smaller VM**

## Step 5: Terminating the instance

When you're finished using the virtual machine, be sure to stop or terminate it, otherwise it will continue to incur costs.

There are two options:

- You can "Stop" the instance (on the console page, click the 3 vertical dots). This will pause the instance, so it's not running, but it will still incur storage costs. This is a good option if you want to come back to this instance (click "Start") without having to reconfigure and download files every time.

- If you're completely done with the instance, you can "Delete" it. This will delete all files though, so [download](./gcp3.md) the files you want to keep!

**TO DO: add image of these options**
