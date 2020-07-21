Editing MkDocs websites with cfde-bot
---------------------------------------

This is a tutorial on how to edit MkDocs websites from Github repos that use robots. The cfde-bots handle changes in each website GitHub repo so you never have to render the website locally to check the edits and most importantly, prevents changes that could break the website.

General steps for editing CFDE websites with cfde-bots
-------------------------------------------------------

1. Go to website's GitHub repository.

    For the CFDE, these are the key web pages we work with: [welcome](https://github.com/nih-cfde/welcome-to-cfde/), [training](https://github.com/nih-cfde/training-and-engagement), [use cases](https://github.com/nih-cfde/usecases/), and [documentation](https://github.com/nih-cfde/published-documentation/). *See section below for specific instructions on editing the documentation website.*

2. Create a new branch, "my_branch" (type new branch name, hit enter)

    <div><img src="https://i.imgur.com/5pv1CQL.png" alt="drawing" width="300" align="top"/></div>

3. You should now be on your new branch, otherwise click `Branch:master` and switch to your new branch. Find the file(s) you want to change and make changes directly in the web interface (it's not necessary to make changes locally). When you’re done editing, scroll to the bottom of the page to merge; it will force you to commit changes to your branch.

4. Then push those changes from "my_branch" to the permanent "preview" branch by clicking on "Compare & pull request". Change the base branch to "preview" instead of "master". Add a description about the pull request (PR) and submit it.

5. The robot will now make sure the website can build with these changes. The PR page will show the following as the bot is going through its checks:

    <div><img src="https://i.imgur.com/Jqdofvs.png" alt="drawing" style="border:1px dotted grey;" width="500"  align="top"/></div>
    then,
    <div><img src="https://i.imgur.com/GN4unNE.png" alt="drawing" style="border:1px dotted grey;" width="500"  align="top"/></div>
    and finally!
    <div><img src="https://i.imgur.com/XljjGj3.png" alt="drawing" style="border:1px dotted grey;" width="500"  align="top"/></div>

    This is what happens if the bot does not approve:
    <div><img src="https://i.imgur.com/nJBrhPl.png" alt="drawing" style="border:1px dotted grey;" width="500"  align="top"/></div>

6. If the checks pass, changes should be automatically merged into the preview branch so you can look at the preview render of the website. Click "Merge pull request" and confirm merge.

7. Merging will re-render on the preview version on the [readthedocs](www.readthedocs.com) website. Login to [readthedocs](www.readthedocs.com), go to website Overview page, and click on "preview" under Versions (*requires admin privileges*). You might need to refresh and wait a few minutes for the changes to show up.

8. To keep those changes, go back to the GitHub repo and click "Compare & pull request" to push "my_branch" to master. This time, request a reviewer for this PR. After approved, merge changes and check that the changes are on the "latest" version of the website!

Specific steps for editing the CFDE documentation website
----------------------------------------------------------

The website created by the `published-documentation` repo pulls some docs that are in its repo AND some from a sub-module (`the-fair-cookbook`). There are two ways to make changes to this website.

### A) **To edit documents that are *in* the `published-documentation` repo**

Follow the general steps above for this repo: https://github.com/nih-cfde/published-documentation/.

### B) **To edit documents that are in the sub-module `the-fair-cookbook` repo**

The cfde-bot's process for checking changes to this repo is slightly different:

- The `published-documentation` cfde-bot monitors the `the-fair-cookbook` repo "master" branch every 6 hours.

- It pulls any changes into the "preview" branch of the `published-documentation` repo as a PR.

- If the website build checks pass, the bot auto-merges changes into the "preview" branch and renders the preview website.

- The cfde-bot simultaneously creates a PR of the changes to the `published-documentation` "master" branch so you can request a reviewer to merge to the "master" branch if you decide the preview looks good.

Steps:

1. Go to `the-fair-cookbook` repo: https://github.com/nih-cfde/the-fair-cookbook

1. Make changes *directly* on the "master" branch.

1. Check progress from the cfde-bot.

1. Check preview website: `https://cfde-published-documentation.readthedocs-hosted.com/en/preview`

1. To keep these changes, you need to request a reviewer for the PR to the "master" branch the cfde-bot created in the `published-documentation` repo. After approved, merge changes and check that the changes are on the "latest" version of the website!
