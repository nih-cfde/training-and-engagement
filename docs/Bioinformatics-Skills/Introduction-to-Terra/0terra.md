# How does Terra work?

In a nutshell, Terra provides you a consolidated, shareable space for storage of data, analysis workflows, and analysis output results. When you submit a workflow job or use interactive coding tools (Jupyter notebooks, Rstudio), Terra does the heavy-lifting of managing your Google Cloud Platform (GCP) virtual machine for the analysis.

The underlying cloud computing services are provided by the Google Cloud (primarily the Compute Engine and Storage Bucket services). For an introduction to the GCP, see our [tutorial](../Introduction-to-GCP/index.md) on setting up a GCP virtual machine.


## Why should I use Terra?

### Scenarios

- don't have access (or long-term access) to a local high performance computer (HPC)
- have HPC access, but takes too long to wait in queue to process data
- have very large datasets that you need to process in the same way and want to scale up the process
- want to use a point & click interface and/or already built data processing workflows (i.e. GATK workflows for variant calling)
- need a secure place to share data/workflows with collaborators


### Interface features
- Interactive web-based GUI & data viz
- Provenance tracking and versioning of workflows, data
- Secure place to share datasets (designed to meet the security standards required for using and storing human data), workflows, and analysis workspaces with collaborators
- Long-term persistence of data/workflow availability to share with the community (e.g., custom workflows can be uploaded and shared on the platform)

### Cloud and scaling
- Scale up sample processing
- Cost optimization of cloud resources

### Data access
- Apply for access to use controlled-access datasets, which in some cases, are only available via cloud platforms like Terra (e.g., Gene-Tissue Expression (GTEx) project [data](https://anvilproject.org/data?query=consortium%3DGTEx%2B%2528v8%2529)
- open datasets available on workspaces in Terra making it easier to pull into workflows (no need to upload or worry about multiple copies of data)






## Finding help on Terra

The [support center](https://support.terra.bio/hc/en-us) and [youtube channel](https://www.youtube.com/channel/UCkXAqpR5Hk1ZmNd2-1K2l5Q/videos) of short tutorial videos provide helpful documentation for all users looking to get started with Terra’s basic functions. Importantly, both are actively updated, so users can be sure they are getting the most up to date information. Users can post questions to the [community help forum](https://support.terra.bio/hc/en-us/community/topics) or directly to Terra by [submitting a request](https://support.terra.bio/hc/en-us/requests/new). There is a [service notification page](https://support.terra.bio/hc/en-us/sections/360003692231-Service-Notifications) for scheduled maintenance or service interruptions/errors related to the Terra platform.
