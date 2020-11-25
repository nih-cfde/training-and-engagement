# Setting up a GCP instance

## Step 1: Create a project

- Under the "IAM & Admin" section, select the "Manage Resources" page

![](../../images/gcp_images/gcp_project1.png "Manage resources tab")

- Click on "Create Project", enter a unique project name, and click "Create".

![](../../images/gcp_images/gcp_project2.png "Create Project button")

- The new project is now listed in the table, along with a project ID. You'll need the project ID later to connect to a virtual machine (VM).

![](../../images/gcp_images/gcp_projectid.png "Project ID")

## Step 2: Configure custom VM

- On the GCP console left-hand panel, scroll down the "Compute Engine" section and select "VM instances". It may take a few minutes to load.

![](../../images/gcp_images/gcp_vm.png "VM instances")

- Click "Create". There are several configuration options to set up:

### a. Name your VM

![](../../images/gcp_images/gcp_vmconfig1.png "VM configuration name and region")

Names must be in lowercase letters, numbers. Use hyphens "-" instead of spaces. You'll need the VM name to connect to it.

### b. Choose a Region

Several regions are available from the dropdown menu. It is less important which zone you choose, and the interface automatically selects a zone based on the region you choose. You'll need the zone name to connect to the VM. We are using the "us-west1 (Oregon)" region.

!!! note "Machine regions"

    **TO DO: add what region and zones mean for GCP**

    This tool can help to choose the region closest to you (http://www.gcping.com/). You may need to refresh several times.

### c. Choose machine type and configuration

![](../../images/gcp_images/gcp_vmconfig2.png "VM configuration machine type")

For this tutorial, select Series "E2" and Machine type `e2-micro`. The estimated monthly cost for each machine type is shown on the top right side panel of the console page next to where you entered the VM name.

!!! note "Machine types"

    **TO DO: add info about what these machine types mean, how to select, etc.**

    **TO DETERMINE: is this a good machine to use?**

### d. Customize boot disk

![](../../images/gcp_images/gcp_vmconfig3.png "VM configuration boot disk")

Click on "Change". The default operating system is Debian, change it to "Ubuntu" and select version "Ubuntu 20.04 LTS". For this tutorial, we'll leave the [persistent disk storage](https://cloud.google.com/persistent-disk) as the default 10Gb. Depending on the tasks you will use the VM for, you may need to increase the storage amount (e.g., to 100Gb).

!!! note "Persistent disk storage"

    **TO DO: add info about what this is for and how to choose**

### e. Firewall setting

![](../../images/gcp_images/gcp_vmconfig4.png "VM configure firewall")

Check the box by "Allow HTTP traffic" under the Firewall configuration, which opens port 80 (HTTP) and allows you to access the virtual machine. If this box is not checked, you will get an error message when trying to connect to the VM ("Insufficient Permission: Request had insufficient authentication scopes.").

### f. Complete configuration

When you're done configuring the VM, click "Create".

## Step 3: Connect to your VM

![](../../images/gcp_images/gcp_vmGCS.png "VM connect with Google Cloud Shell")

### a. Check box by VM

### b. Open Google Cloud Shell

The GCP console provides a free Google Cloud Shell. This shell environment is useful for small tasks that do *not* require a lot of CPU or memory (as most bioinformatic analyses do). For example, it is a good place to learn how to use the Google shell environment without incurring cost or to access Google Cloud products (e.g., a Google Storage bucket or GCP virtual machine).

Click on the "Activate Cloud Shell" icon. A new panel will open on the bottom half of your screen. The first time you start the shell, you'll need to agree to the Google Cloud terms of service and privacy policy. After starting the shell, it may take a few minutes to connect. Go to the Google support [documentation](https://cloud.google.com/shell/docs/using-cloud-shell) for more information.

The Google Cloud Shell command prompt format will show: "<username@cloudshell:~ (<project id>)$".

### c. Connect to VM

Use the `gcloud compute` command to connect to your virtual machine (`gcloud` is a tool from the Google Cloud SDK toolkit). You'll need your project ID, zone, and instance name. In the command below, the project ID is `"my-first-project-296523"`, the zone is `us-west1-b`, and the instance name is `test-vm`. Replace these values to run the command for your virtual machine.

```
gcloud compute --project "my-first-project-296523" ssh --zone us-west1-b test-vm
```

### d. Set up SSH keys

Next, set up the SSH public/private keys. Follow the prompts in the terminal:

- "Do you want to continue (Y/n)?": type ++Y++
- "Enter passphrase (empty for no passphrase)": the passphrase can be left empty, type ++enter++
- "Enter same passphrase again": type ++enter++

When this process is complete, your command prompt in the terminal should switch to your "<user name>@<VM name>:~$". You're now in the VM space!
