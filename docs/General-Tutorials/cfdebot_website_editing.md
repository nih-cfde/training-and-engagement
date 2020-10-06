# Editing MkDocs websites with cfde-bot

This is a tutorial on how to edit MkDocs websites from Github repos that use the cfde-bot robot. The cfde-bot handles changes in each website GitHub repo to check  edits and most importantly, prevent changes that could break the website.

!!! note "Learning Objectives"

    - learn how to edit Mkdocs website docs on Github using the cfde-bot

=== "Est. Time"
    
    ~20 mins (depending on the website, it make take longer for all the bot's checks to complete)

=== "Prerequisites"
    
    The cfde-bot was constructed for managing CFDE websites. To edit any of these websites, you must be onboarded to the CFDE to gain access to the nih-cfde Github repo files. 

## General steps for editing CFDE websites with cfde-bots

For the CFDE, the key websites we maintain are: `welcome-to-cfde`, `training-and-engagement`, `usecases`, and `published-documentation`. *See section below for specific instructions on editing the `published-documentation` website.*
    
The `welcome-to-cfde`, `usecases`, and `published-documentation` CFDE websites currently have a single PR and review stage to merge changes to the public-facing website. For these websites, the public-facing website is rendered from the `stable` branch.
    
The `training-and-engagement` repo has a 2-stage PR and review stage. Changes are initially pushed to a development branch `dev` and then merged to `stable` for new releases of the public-facing website. Consequently, the steps of the PR process are slightly different compared to the other CFDE websites. The details are available in the website repo's style guide and release plan. 
    
In this tutorial, we lay out the general steps for the single stage PR process.

### Step 1: Go to website's GitHub repository.

### Step 2: Create a new branch

If working from Github, type a new branch name (e.g., `my_branch`) in the `Branch:stable` dropdown button and hit ++enter++. The new branch should be created with the `stable` branch as its base.

![](../images/github-branch-stable.png "create new github branch")

### Step 3: Make edits

You should now be on your new branch, otherwise click `Branch:stable` and switch to your new branch. 

Find the file(s) you want to change and make changes directly in the web interface (while not necessary, you can also make changes via a local copy of the repo with `git` commands). When you’re done editing, scroll to the bottom of the page to commit the changes; Github knows to commit the changes to your branch.

### Step 4: Create pull request

Push your changes from `my_branch` to the `stable` branch by clicking on "Compare & pull request". Add a description about the pull request (PR) and submit it. Depending on the website, follow style guidelines for the PR and review process.

### Step 5: Wait for cfde-bot checks

The cfde-bot will now make sure the website can build with these changes. The PR page will update the progress as the bot goes through its checks, starting with a yellow circle and "Some checks haven't completed yet". If the checks complete successfully, the circle turns green with a check mark and the message switches to "All checks have passed". If any checks failed, there will be a red circle with an "X" mark and the message "Some checks were not successful". In the latter case, check the build logs to track the error.

### Step 6: Check rendered changes

If the checks pass, changes should be automatically merged into a preview version of the readthedocs website. Click on the "Details" button of the "docs/readthedocs.com:<Github repo> — Read the Docs build succeeded!" check to view the rendered preview (the website is rendered for the specific PR).

### Step 7: Request reviews

If you are satisfied with the edits, request reviewers to check, request changes, and approve the changes. 

### Step 8: Merge changes

When the PR has been approved, click "Merge pull request" and confirm merge. Be sure to delete your branch. Your edits should now be viewable on the public-facing website!

## Specific steps for editing the CFDE documentation website

The website created by the `published-documentation` repo pulls some docs that are in its repo AND some from a sub-module (`the-fair-cookbook`). There are two ways to make changes to this website.

### A) **To edit documents that are *in* the `published-documentation` repo**

Follow the general steps above. Note that this repo currently labels branches as `master` instead of `stable`. This will be updated in the near future.

### B) **To edit documents that are in the sub-module `the-fair-cookbook` repo**

The cfde-bot's process for checking changes to this repo is slightly different:

- The `published-documentation` cfde-bot monitors the `the-fair-cookbook` repo `master` branch every 6 hours. Note that `the-fair-cookbook` website is not one of the websites that the Brown lab maintains, though we may contribute material. Also, this repo currently labels branches as `master` instead of `stable`.

- It pulls any changes into the `preview` branch of the `published-documentation` repo as a PR.

- If the website build checks pass, the bot auto-merges changes into the `preview` branch and renders the preview website.

- The cfde-bot simultaneously creates a PR of the changes to the `published-documentation` `master` branch so you can request a reviewer to merge to the `master` branch if you decide the preview looks good.

#### Step 1: Go to `the-fair-cookbook` repo: https://github.com/nih-cfde/the-fair-cookbook

#### Step 2: Make changes *directly* on the `master` branch.

#### Step 3: Check progress from the cfde-bot.

#### Step 4: Check preview website

The website link is https://cfde-published-documentation.readthedocs-hosted.com/en/preview. You must have admin permissions to acces the readthedocs website.

#### Step 5: Request review 

To keep these changes, you need to request a reviewer for the PR to the `master` branch the cfde-bot created in the `published-documentation` repo. 

#### Step 6: Merge changes

After approved, merge changes and check that the changes are on the public version of the website!
