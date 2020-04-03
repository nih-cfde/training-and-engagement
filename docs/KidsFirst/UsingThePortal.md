# Kids First Data Portal Training

[toc]

### Let's Join and Get Logged In 
We will need to make two accounts: 
(1) Kids First Data Portal (KFDP)
(2) Cavatica 
This is for our ultimate goal of analyzing KF data in Cavatica. 
KFDP is home to the *data*. Cavatica is home to the *analysis tools*. 
You will also see Cavatica called the "cloud computing resource."

#### Kids First Data Portal
kidsfirstdrc.org 
Click 'Join now'
Fill out your info and log in
Logon choices:
* ORCID
* Google
* Facebook

#### Cavatica
cavatica.squarespace.com
Click "Access Data"
Fill out your info to create an account
Cavatica will send you an email to verify your account, please logon to your email and accept.
Logon choices:
* eRa commons
* email


All set? Great! Now that we have accounts, let's start familiarizing ourselves with KFDP and Cavatica. We're going to start with KFDP. 

### Kids First Data Portal Orientation
#### Portal Home Page
1. #### Site Overview

Let's take a look at the KFDP Homepage. At the top we see this bar: 

Does your page look like this? F5 (Windows) or Cmd+Shift+R (MacOS) to refresh, or the refresh button.

![](https://i.imgur.com/mAnoGT1.png)

Where do we start? Remember: the data we want lives on KFDP but our tools live in Cavatica. So the first thing we're going to do is make KFDP and Cavatica talk to each other. This is what KFDP calls "Application Integration."

* Click on your name (top right) --> Select Settings --> Scrool down to Application Integration --> Click "Connect"

Need authentication token for API endpoint from Cavatica
Go back to Cavatica --> Developer --> Authentication token --> Generate --> copy and paste this code into KFDP 

All set? Not quite. Even though KFDP and Cavatica are now linked, we need to give *permission* for the KFDP data to get pushed over to Cavatica, our cloud computing resource. This is what KFDP calls "Data Repository Integration."

Need eRa commons, NIH login. 
***This is a holdup for trainees. 
***Can we skip this part? It's for controlled access data. 

* From the same page as the Application Integration, scroll up slightly and choose "Connect" for both the KF Framework Services and NCI CRDC. 
* You will need your NIH logon credentials. 

2. #### Dashboard Overview

Click on the KF logo (top left) to get back to the KFDP Homepage. 
You'll notice that besides the navigation bar at the top of the page, we also have a Dashboard. Let's take a closer look at what the Dashboard can tell us. 

![](https://i.imgur.com/2I0WsES.png)

Since this is our first time on the portal, let's start with the bottom row. The categories on this row are generalized to the KFDP as a whole. You'll notices each topic has an associated graphic. Mouse over each of the graphics below to get an idea of what is studied using KFDP: 
* Studies/Participants
* Member Research Interests
* Most Frequent Diagnoses

Now let's look at the top row. Here we see information specific to your own portal experience. 
* Saved Queries
* Authorized Studies
* Cavatica Projects

#### Browsing the Portal

Now that we're all connected, let's explore KFDP. Our goal is twofold: 
(1) Find data that is of interest to us, and 
(2) Get that data onto Cavatica

* #### Explore Data

    * This is where we can get a bird's eye view of the data available. It is important to get a sense of what KFDP houses prior to designing an experiment. 
    * The dashboard on this page gives us **Cohort Results** based on the dataset features we select.
    * The most expedient way to do this is to click the **Quick Filters** button. Let's try that now.
    * Choose some different options from the **Quick Filters** category and watch how your **Cohort Results** change. 
    * To further refine your **Cohort Results** select more categories along the navigation bar.
    * Clicking on the graphics also adds filters to the data
    * Need to cover how to remove the filter

![](https://i.imgur.com/HIuuO67.png)



* #### File Repository
    * Now that we have a sense of what data is on KFDP, let's choose some data to work with. 
    * **The most important feature of data is whether or not you have access to it** 
        * Filter (top left) --> Browse All --> Access --> Click "Open" --> View Results
        * Now we are searching exclusively Open data
    *  Let's select some data from the **Pediatric Brain Tumor Atlas: Childrens' Brain Tumor Tissue Consortium (CBTTC)**
    *  You will see that we can choose Clinical Filters and/or File Filters
    *  First, lets choose *Cancer* from "Diagnosis Category" and *High Grade Glioma* from "Diagnosis" from Clinical Filters.
    *  That's enough Clinical Filters for now. Let's see what this returns for us in terms of available files. *It's important to not get too specific until we know what's available.* 
    *  Under File Filters, we see two types of "Experiment Strategies", WGS and RNASeq. Click each one at a time and see how the File Filter categories below change. 

![](https://i.imgur.com/UtjJeYp.png)



#### File Format Sidebar: 
* I want to go over what data/file types are associated with what type of analysis
* I also want to go over what file types can actually be used for further downstream analysis
* Example, the RNASeq data returns .tsv, .pdf and "Other"
    * ie: all files returned are not going to be useful in Cavatica, but may be useful in other ways

![](https://i.imgur.com/MIlY0DJ.png)


#### Selecting Files and Pushing Files to Cavatica
* We are currently looking at data under the CBTTC. The data has the following features: 
    * Cancer
    * Patients diagnosed with high grade glioma
    * RNASeq 
* Under the File Filter "File Format" choose .tsv
* On the main screen, scroll to the bottom and select the number of files you wish to view at one time
* On the main screen, click the small box to the left of "File ID" to select all of the .tsv files
* Now the exciting part! Click the large purple **"ANALYZE IN CAVATICA"** button

![](https://i.imgur.com/So5jl7k.png)


* We see a dialogue box pop up. In this case, we are authorized to copy all 100 files we selected, but we need a Cavatica project to assign them to. 
    * Click "+ Create a project" and name your project "{Last name} practice lesson"

![](https://i.imgur.com/Df1aDGx.png)

*Note: can currently only copy the files to Cavatica that you have selected. So, if you have 400 files, but can only select the first 100, you will have to copy 100 files at a time.

*Note: There is a message button that pops up that says the files were transferred successfully. But I haven't been successful yet. 

#### Confirm KF data can get onto Cavatica
* Log onto Cavatica, navigate to your project and look for your files. Do you seem them? 
* This is currently a problem
* Last night my 100 RNASeq .tsv files didn't transfer (again)
* FAQ guidance: In Cavatica --> click name (upper right) --> account settings --> dataset access
    * Look for green check mark next to KF
    * I currently have a red X. Not sure who to contact to change this to a green check mark. 
* There appear to be CBTTC files in Cavatica but are not the same as on KFDP. For example, the RNASeq files are now bam files.

### Random Final Thoughts

#### "About the Portal"
* Why is this such a tiny button in the middle of nowhere? There's some good stuff in here, but users are unlikely to click on it. 

#### Here's how to apply for data access
* I haven't made this section yet 
