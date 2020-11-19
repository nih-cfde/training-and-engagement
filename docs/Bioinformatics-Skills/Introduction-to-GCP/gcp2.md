# Configure custom VM

- need to make new demo files for the bucket steps

**The Google Cloud Shell is free, but doesn't have enough CPU or memory to do bioinformatics analyses. Here's how to set up a custom VM with more CPU/memory and other configurations.**

1. On the GCP console, go to "VM instances" in the "Compute Engine" section. It may take a few minutes to load up.
![image](https://user-images.githubusercontent.com/5659802/96937594-57b6d700-147d-11eb-9220-321151fc5b98.png)
2. Click "Create". There are several configuration options to set up
![image](https://user-images.githubusercontent.com/5659802/96938921-72d71600-1480-11eb-84c4-cf4134eca73b.png)

- name your VM
- choose Region (less important which zone you choose). This tool can help to choose the region closest to you (http://www.gcping.com/). You may need to refresh several times.
- choose machine type. The cost for each machine type is shown on the right side panel of the console page. From the book example, they set used Series N1 and machine type `n1-standard-2`. There are more options on the current interface. Series N1 isn't an option anymore, it's now N2. I chose `n2-standard-2`.
- customize boot disk. Click on "Change". The default OS is Debian, change it to Ubuntu. For the version, to match AWS tutorials, I chose Ubuntu 20.04 LTS. Set the amount of persistent disk storage you want (the book suggests 100Gb for its tutorials; I set it at 50Gb for testing).

When you're done configuring, click "Create".

3. Click on the "SSH" dropdown option for your new instance and select "Open in browser window"
![image](https://user-images.githubusercontent.com/5659802/96941450-4d013f80-1487-11eb-9951-2ccc6950a899.png)

A new window will open with the instance terminal. The gear icon at the top right can be used to customize terminal window settings (e.g., appearance, text size). You can also access the instance via command line, e.g., from the Google Cloud Shell or local laptop terminal (but you then need to install/set up/authorize `gcloud`):
```
gcloud compute ssh --project [PROJECT_ID] --zone [ZONE] [INSTANCE_NAME]
```
4. Set up gcloud authorization so you can move files to/from your instance.

- `gcloud init`
- type "2"
- click the link that appears on the terminal. A new web browser page will open, log in with your GCP google account
- copy/paste the verification code back to the terminal
- enter the number that corresponds to your project
- you can configure a default Compute Region and Zone or not

*After this step, you now have a custom configured VM. Things we could do at this point are:*
- [upload or download files](./gcp3.md)
- [pull in a docker container](./gcp4.md) to get software without a bunch of software installations
- do some analysis





## Terminating the instance

When you're finished using the virtual machine, be sure to stop or terminate it, otherwise it will continue to incur compute engine and/or storage costs (?).

There are two options:

- You can "Stop" the instance (on the console page, click the 3 vertical dots). This will pause the instance, so it's not running, but it will still incur storage costs. This is a good option if you want to be able to come back to the instance (click "Start") without having to reconfigure and download files every time.

- If you're completely done with the instance, you can "Delete" it. This will delete all files though, so download whatever you want to keep!

ADD IMAGE HERE - can i show both stop and delete and add numbers to them?
