---
layout: page
title: Adding Data to CAVATICA with DRS IDs
---

# Adding data to CAVATICA with DRS IDs

[CAVATICA](https://cavatica.sbgenomics.com/) is a storage, sharing, and analysis platform designed by Seven Bridges. This tutorial shows how to add data to CAVATICA using DRS IDs obtained from [the previous lesson](find-export/). Please consult the [CAVATICA Documentation](https://docs.cavatica.org/docs) for a full description of its features.

## Login

Login to CAVATICA https://cavatica.sbgenomics.com/home

![drs1](https://i.imgur.com/KlGNa2Y.png)

Click "Yes, I authorize" to allow CAVATICA to view your Researcher Auth Service (RAS) information

![drs2](https://i.imgur.com/E3LjLsr.png)

## Create a Project

Create a new project by clicking "+ Create Project". Provide a name and billing group then click "Create".  

## Add files

Click the "Files" tab then click "+ Add files" drop down menu select "GA4GH Data Repository Service (DRS)".  

![drs-3](https://i.imgur.com/zWX7W9M.png)

Enter the DRS identifiers (separated by a comma, enter, or tab) and click "Import"
* _(Follow the steps in [the pentultimate lesson](./find-export) to get 2 DRS identifiers)_
* drs://drs.hmpdacc.org/mZBm6TYnDQoS
* drs://drs.hmpdacc.org/1CzGt5DyGwcOe

Wait for the import to finish then click "Download" to download a raw file locally confirm that it is the correct format.

![drs-4](https://i.imgur.com/mLcE2wZ.png)
