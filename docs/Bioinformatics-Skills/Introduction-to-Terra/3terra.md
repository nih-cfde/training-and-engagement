# Let's talk about costs

The three elements of cloud computing the incur cost on Terra are:

1. running an analysis (workflow, jupyter notebook, Rstudio)
2. downloading files from Terra (i.e., to a local computer); also known as egress
3. storing files on the Google Cloud Storage bucket associated with Terra workspaces

## Running a workflow

For small data sets, analysis costs are quite reasonable. In our tests, running through two tutorials cost less than $1. However, Terra does not display any cost tracking information. Billing is only available in the Google Cloud Platform (GCP) console, and only “Billing Account Administrator” roles have access to billing reports. We could not find a way for users to track their costs unless they are added as “Billing Account Viewer”. This means that for multi-user institutions like universities or hospitals, most users will not know how much they’ve spent unless they ask the credit card holder, or are specifically assigned the extra permissions to do so.

Conversely, Terra does provide several ways for users to assess their upcoming costs. Data analysis price rates are set by Google, and rates differ slightly depending on the cloud service region’s zone and virtual machine type. To help navigate those differences, Terra offers several support articles about “Controlling costs (cloud storage, compute, and egress)”, and the “Cloud environment” configuration panel shows the cost per hour for the various cloud compute configurations. Similarly, the quickstart and featured workspaces on Terra often state the estimated costs of running through a given workflow. For example, this demo exercise will cost $0.29 to execute:


## Egress

Users can preview the cost to download individual data files before downloading, and in our tests these costs were low, typically less than a penny per file. However, while this per file cost is low, there are two important caveats. First, Terra does not have native display capabilities for all file types, and download may be the only way to preview/view some files. So egress charges will be unavoidable for some users. Second, our test files were created by running through tutorials, which are typically made with much smaller and far fewer files than a real analysis would use. Egress expenses associated with downloading local copies of data or moving data between storage cloud regions, typically cost $0.01/Gb in the United States, and it will be important to assess the volume of data that must be egressed in real life analysis pipelines before promoting this platform for all users.

## Storage
Data storage costs are set by Google Cloud Platform price rates, and are currently low but variable. Data storage prices are more expensive for data that need to be accessed quickly and frequently - $0.026/Gb/month (multi-regional) to $0.020/Gb/month (regional). Data that is accessed infrequently (e.g., archived processed sequence files) can be stored on Nearline or Coldline storage systems, where the prices range from $0.01/Gb (Nearline) to $0.007/Gb (Coldline) for storage and $0.01/Gb (Nearline) to $0.05/Gb (Coldline) for data retrieval.



In the next lesson, we'll run an existing workflow on Terra.
