# Demo: Custom Workflows

In this section, we will demonstrate how we built the workflow used in the last lesson to run an example GWAS analysis on Terra. The docker and workflow script files were tested locally on a MacOS computer.

- Part 1: Create docker containers with specific software environments for steps of the workflow
- Part 2: Write a workflow script to specify the analysis steps
- Part 3: Test the workflow before uploading it to Terra
- Part 4: Set up workflow on Terra

To see how this process works, watch the vidlets below! Each vidlet is followed by written descriptions about the steps.

!!! important

    It's beyond the scope of this lesson to teach general docker and Workflow Description Language (WDL) syntax, however, check out the resources linked throughout the lesson.

    We hope this demo helps to show the general steps and get you started with your own workflows on Terra!

## Part 1: Setting up Docker

We are using docker containers to define isolated software environments for specific steps of the GWAS workflow. Read more about using [docker](https://support.terra.bio/hc/en-us/articles/360037340472-Docker-container-overview) from Terra's support docs.

The core steps of building, testing, and sharing the docker container are shown in the vidlet below:

<iframe id="kaltura_player" src="https://cdnapisec.kaltura.com/p/1770401/sp/177040100/embedIframeJs/uiconf_id/29032722/partner_id/1770401?iframeembed=true&playerId=kaltura_player&entry_id=1_t7exj3sj&flashvars[localizationCode]=en&amp;flashvars[leadWithHTML5]=true&amp;flashvars[sideBarContainer.plugin]=true&amp;flashvars[sideBarContainer.position]=left&amp;flashvars[sideBarContainer.clickToClose]=true&amp;flashvars[chapters.plugin]=true&amp;flashvars[chapters.layout]=vertical&amp;flashvars[chapters.thumbnailRotator]=false&amp;flashvars[streamSelector.plugin]=true&amp;flashvars[EmbedPlayer.SpinnerTarget]=videoHolder&amp;flashvars[dualScreen.plugin]=true&amp;flashvars[Kaltura.addCrossoriginToIframe]=true&amp;&wid=1_3nzotpfs" width="608" height="402" allowfullscreen webkitallowfullscreen mozAllowFullScreen allow="autoplay *; fullscreen *; encrypted-media *" sandbox="allow-forms allow-same-origin allow-scripts allow-top-navigation allow-pointer-lock allow-popups allow-modals allow-orientation-lock allow-popups-to-escape-sandbox allow-presentation allow-top-navigation-by-user-activation" frameborder="0" title="Kaltura Player"></iframe>

More information about each step:

=== "1. Install docker"

    We followed installation instructions for Mac: <https://docs.docker.com/docker-for-mac/install/>. See Terra docs for more about [installing docker](https://support.terra.bio/hc/en-us/articles/360036000631-Install-Docker-and-test-that-it-works).

=== "2. Dockerhub account"

    In order to share the docker, we are pushing the image to Dockerhub. We created an account at <https://hub.docker.com/>. We need the account username to run docker commands.

    There are other options (i.e., Dockstore). For the workflows used in Terra, either of these platforms will work!

=== "3. Dockerfile"

    We told docker what we wanted installed in the container using a file called ["Dockerfile"](./Dockerfile.md). Read Docker's docs for [best practices for writing Dockerfiles](https://docs.docker.com/develop/develop-images/dockerfile_best-practices/) and Terra's lesson on [writing a Dockerfile](https://support.terra.bio/hc/en-us/articles/360024737591-Make-a-Docker-container-image-the-easy-way-using-a-base-image).

    We installed dependency software tools, vcftools, and the R package "qqman". There's an existing docker container for plink (`gelog/plink:latest`), so we use that in the workflow for plink steps.

    The base computer image of our container is ubuntu 20.04, and the software installation commands are similar to those used in the [GWAS lesson](../GWAS-in-the-cloud/index.md), with some modifications for docker syntax.

    - `WORKDIR /usr` specifies the working directory - since dockers are closed containers, it helps to know where commands will be run so you can locate output files.
    - we use `ENV DEBIAN_FRONTEND noninteractive` because docker doesn't allow interactive installation.
    - `RUN` precedes all software install commands, and strings of commands are connected by `&&`. The `\` allows commands to be written over multiple lines in the Dockerfile file. After the first `RUN` line, there is an indent for the others in the installation block.

=== "4. Build"

    We built the docker with a specified tag list (`-t`). In this case, the tag list contains a user name (`demo`), the docker name (`demo_gwas`), and the version (`tag0`). The last element of the command specifies the location of our "Dockerfile"; since we're in the directory that contains it, we can use the `.` or `./` notation.

    ```
    docker build -t demo/demo_gwas:tag0 ./
    ```

=== "5. Test"

    !!! warning

        The docker app must be running to use docker. In the vidlet, you can see that a docker container (`kind_elion`; it's a randomly generated name!) pops up in the docker window (top left).

    We set up an interactive (`-it`) session of our docker container with `docker run`:

    ```
    docker run -it demo/demo_gwas:tag0
    ```

    And ran a test command to make sure the software `vcftools` was properly installed:

    ```
    vcftools -h
    ```

=== "6. Share"

    Pushing to Dockerhub takes a few minutes to complete:

    ```
    docker push demo/demo_gwas:tag0
    ```

    After pushing the container, it is publicly available on Dockerhub.

## Part 2: Write a WDL

Workflow Description Language (WDL, pronounced "widdle") is one of many computer programming languages used to define analysis workflows. It is maintained by an [open online community](https://openwdl.org/) that you can contribute to!

There are many resources from Terra to help you get started with the Workflow Description Language:

- [Step-by-step tutorials for writing WDL pipelines](https://support.terra.bio/hc/en-us/sections/360007347652-WDL-Tutorials)
- [Youtube playlist with WDL lessons](https://www.youtube.com/watch?v=RtcW2Zdn_28&list=PL4Q4HssKcxYv5syJKUKRrD8Fbd-_CnxTM)

To test the workflow on local computer, we need 2 files:

=== "1. WDL script"

    The WDL script we wrote: [GWAS_WDL.wdl](./GWAS_WDL.md)

    WDL workflow structures are nested. This script has a main `workflow` block that specifies the steps (`call`) of the workflow to run. Each step is specified by tasks (`task`) that include input, command, output, and runtime information. We created a task for each step of the GWAS workflow.

=== "2. Inputs JSON"

    The JSON file we wrote to specify inputs for each workflow step: [GWAS_WDL.json](./GWAS_json.md)

    Note that this file is only needed for local testing of the WDL. Terra has a graphical user interface for specifying inputs/outputs.

## Part 3: Testing a WDL workflow

<iframe id="kaltura_player" src="https://cdnapisec.kaltura.com/p/1770401/sp/177040100/embedIframeJs/uiconf_id/29032722/partner_id/1770401?iframeembed=true&playerId=kaltura_player&entry_id=1_vzjg3it3&flashvars[localizationCode]=en&amp;flashvars[leadWithHTML5]=true&amp;flashvars[sideBarContainer.plugin]=true&amp;flashvars[sideBarContainer.position]=left&amp;flashvars[sideBarContainer.clickToClose]=true&amp;flashvars[chapters.plugin]=true&amp;flashvars[chapters.layout]=vertical&amp;flashvars[chapters.thumbnailRotator]=false&amp;flashvars[streamSelector.plugin]=true&amp;flashvars[EmbedPlayer.SpinnerTarget]=videoHolder&amp;flashvars[dualScreen.plugin]=true&amp;flashvars[Kaltura.addCrossoriginToIframe]=true&amp;&wid=1_449aslzx" width="608" height="402" allowfullscreen webkitallowfullscreen mozAllowFullScreen allow="autoplay *; fullscreen *; encrypted-media *" sandbox="allow-forms allow-same-origin allow-scripts allow-top-navigation allow-pointer-lock allow-popups allow-modals allow-orientation-lock allow-popups-to-escape-sandbox allow-presentation allow-top-navigation-by-user-activation" frameborder="0" title="Kaltura Player"></iframe>

=== "1. Install tools"

    There are a few tools for testing WDL workflow scripts. One option is to use [WOMtool](https://cromwell.readthedocs.io/en/stable/WOMtool/) and [Cromwell](https://github.com/broadinstitute/cromwell/), which we'll use in this demo, and another is to use [miniwdl](https://github.com/chanzuckerberg/miniwdl).

    The commands for installing Cromwell and WOMtool apps on MacOS are:

    ```
    wget https://github.com/broadinstitute/cromwell/releases/download/47/cromwell-47.jar
    wget https://github.com/broadinstitute/cromwell/releases/download/47/womtool-47.jar
    ```

    These tools require Java ([installation instructions](https://www.oracle.com/java/technologies/javase-jdk15-downloads.html)). We checked the java installation with this command:

    ```
    java --version
    ```

=== "2. Inputs"

    The 2 input files for this workflow are:

    - "coatColor.pheno"

    - "pruned_coatColor_maf_geno.vcf.gz"

    We added the workflow input files to the directory that contains:

    - the WOMtool app ("womtool-47.jar")

    - the Cromwell app ("cromwell-47.jar")

    - the WDL script ("GWAS_WDL.wdl")

    - the input JSON file ("GWAS_WDL.json")

=== "3. WOMtool"

    We used the WOMtool app to test the WDL script syntax. Note that the workflow may still fail to run, but WOMtool helps to catch WDL syntax errors. We have tested the script to make sure it runs, so the output of the command below is "Success!".

    ```
    java -jar womtool-47.jar validate GWAS_WDL.wdl
    ```

=== "4. Run WDL"

    !!! warning

        The docker app must be running to use docker. In the vidlet, you can see that a docker container (`bold_wilbur`) pops up in the docker window (top left).

    We used the Cromwell app (".jar" file) to run the WDL script. We have set up the input (`-i`) files for each workflow task in the ".json" file.

    ```
    java -jar cromwell-47.jar run GWAS_WDL.wdl -i GWAS_WDL.json
    ```


### WDL Outputs

<iframe id="kaltura_player" src="https://cdnapisec.kaltura.com/p/1770401/sp/177040100/embedIframeJs/uiconf_id/29032722/partner_id/1770401?iframeembed=true&playerId=kaltura_player&entry_id=1_7ggwnimg&flashvars[localizationCode]=en&amp;flashvars[leadWithHTML5]=true&amp;flashvars[sideBarContainer.plugin]=true&amp;flashvars[sideBarContainer.position]=left&amp;flashvars[sideBarContainer.clickToClose]=true&amp;flashvars[chapters.plugin]=true&amp;flashvars[chapters.layout]=vertical&amp;flashvars[chapters.thumbnailRotator]=false&amp;flashvars[streamSelector.plugin]=true&amp;flashvars[EmbedPlayer.SpinnerTarget]=videoHolder&amp;flashvars[dualScreen.plugin]=true&amp;flashvars[Kaltura.addCrossoriginToIframe]=true&amp;&wid=1_r2c3qxy8" width="608" height="402" allowfullscreen webkitallowfullscreen mozAllowFullScreen allow="autoplay *; fullscreen *; encrypted-media *" sandbox="allow-forms allow-same-origin allow-scripts allow-top-navigation allow-pointer-lock allow-popups allow-modals allow-orientation-lock allow-popups-to-escape-sandbox allow-presentation allow-top-navigation-by-user-activation" frameborder="0" title="Kaltura Player"></iframe>

If the workflow completes successfully, it will display the "Final Outputs" paths and show the message "workflow finished with status 'Succeeded'" in the terminal. Each `task` gets a `call-<task name>` directory:

```
call-run_vcftools
call-plink_missing_rates
call-plink_binary
call-plink_association
call-run_R
```

For each task, the inputs are in an `inputs` directory and the outputs are in an `execution` directory, which also includes helpful files for troubleshooting errors when testing the WDL script (log and stderr files).

The final output of this particular GWAS workflow is a Manhattan plot ("coatColor_man.bmp").


## Part 4: Set up workflow on Terra

After ensuring the WDL works, we followed these [steps](https://support.terra.bio/hc/en-us/articles/360031366091-Create-edit-and-share-a-new-workflow#h_bc0175df-adb7-422f-b2fe-efab18fd598b) from the Terra documentation to set up our workflow on Terra. There is also a [video lesson](https://www.youtube.com/watch?v=VtKlYqWBW6A&ab_channel=Terra) from Terra about workflows.

<iframe id="kaltura_player" src="https://cdnapisec.kaltura.com/p/1770401/sp/177040100/embedIframeJs/uiconf_id/29032722/partner_id/1770401?iframeembed=true&playerId=kaltura_player&entry_id=1_4k2vtulg&flashvars[localizationCode]=en&amp;flashvars[leadWithHTML5]=true&amp;flashvars[sideBarContainer.plugin]=true&amp;flashvars[sideBarContainer.position]=left&amp;flashvars[sideBarContainer.clickToClose]=true&amp;flashvars[chapters.plugin]=true&amp;flashvars[chapters.layout]=vertical&amp;flashvars[chapters.thumbnailRotator]=false&amp;flashvars[streamSelector.plugin]=true&amp;flashvars[EmbedPlayer.SpinnerTarget]=videoHolder&amp;flashvars[dualScreen.plugin]=true&amp;flashvars[Kaltura.addCrossoriginToIframe]=true&amp;&wid=1_2jjqquzc" width="608" height="402" allowfullscreen webkitallowfullscreen mozAllowFullScreen allow="autoplay *; fullscreen *; encrypted-media *" sandbox="allow-forms allow-same-origin allow-scripts allow-top-navigation allow-pointer-lock allow-popups allow-modals allow-orientation-lock allow-popups-to-escape-sandbox allow-presentation allow-top-navigation-by-user-activation" frameborder="0" title="Kaltura Player"></iframe>

=== "1. Upload workflow"

    On Firecloud (aka Broad Method Repository):

    - namespace = "cfde-workflows"
    - name = "GWAS-demo"
    - we copy/pasted the WDL script, but note that we do not need to upload the JSON file because the input files will be specified with Terra's interface
    - in "Permissions", we set the workflow to be "Publicly Readable" so you can use it!

=== "2. Export workflow"

    To export the new workflow from Firecloud to workspace, we used a blank configuration, kept root as "participant" (this is how inputs are identified in the Terra interface), and selected the workspace to export to.

After this point, the workflow is ready to be used, as shown in the [previous section](./4terra.md)!

In the next section, we'll talk about cloud costs.
