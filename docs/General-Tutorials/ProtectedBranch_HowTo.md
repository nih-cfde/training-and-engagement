---
title: "Working with Protected Branches"
layout: page
---

Working with Protected Branches
============================

Why use protected branches
--------------------------

Protected branches ensure that rules are enforced on any changes made to
that branch in a repo. A common branch protection rule is for pull requests (PRs) to be reviewed by at least one other person before they get merged. Protected branches stop you
from making unauthorized changes to that branch. However, you can make multiple changes locally or to other non-protected branches. For the local changes to reflect in the protected branch all the set criteria should be met and the changes should be purposefully merged.
This type of set up is useful for repositories (repos) that render websites. The branch that hosts the live site can be protected, so that no one can incorporate changes to the
site until someone else checks them. That's how we make sure that no one
accidentally breaks the websites.

How to work with protected branches
-----------------------------------

For the purposes of this tutorial, a practice repo will be used to showcase the git commands. First, copy the repo(s) to your local computer using command line like so:

```
git clone https://github.com/nih-cfde/play-with-github.git
```

!!! note "Repo name"
    The example repo is called "play-with-github". You can replace the GitHub URL with any other repo of your choice. The `.git` extension designates a [bare repo](http://www.saintsjd.com/2011/01/what-is-a-bare-git-repository/)

Once you have cloned your repos (now called directories on your local
computer) and want to get started, you need to navigate into the newly created directory and work on a [branch](https://github.com/nih-cfde/organization/blob/master/GitHubUsage.md#definitions).

The command to create a new branch is as follows (you can name your
branch whatever you like):

```
git branch <name_of_new_branch>
git checkout <name_of_new_branch>
```

!!! note
    It is critical to navigate to the directory in which you wish
    to make these changes.

The first line of code creates the branch and the second line of code switches your current directory to the new branch. Now if you make changes to the files, the changes will appear as this new branch. Please note: when you change to a new branch, nothing
will physically change on your computer. But GitHub will recognize your directory as a new branch.

Once you have made the necessary changes (or any changes at all), you
can push changes to GitHub. Here are the basic commands for pushing changes made to a  file:

```
git add <filename>
git commit -m <your_message>
git push --set-upstream origin <name_of_new_branch>
```

If changes have been made to multiple file, you can first check the status using:

```
git status
```
You can use `--all` flag to stage all the changes or list multiple filenames:

```
git add --all
git add <filename1> <filename2> <filename3>
```

!!! note "git commit -am"
    Adding the `-a` flag saves the changes made on tracked files ONLY.

The first line of code adds your changes and the second line of code
saves your changes. The last line pushes your changes to the main GitHub
repo as a new branch.

Until your [pull request](https://github.com/nih-cfde/organization/blob/master/GitHubUsage.md#definitions) is [merged](https://github.com/nih-cfde/organization/blob/master/GitHubUsage.md#merging-pull-requests), you can continue to work on the same branch and push multiple changes. It is not required to create a new branch for ever change. However, after your pull request has been merged and the branch is deleted, any new changes will require to be tracked on a new branch, which incidentally can have the same branch name as previously used.

!!! note "Important housekeeping notes"
    -   It is important to keep the workspace clean by deleting abandoned
        branches and branches that have already been merged. In addition to
        cluttering the workshop, abandoned branches can cause collisions
        with new work.
    -   As you work on your branch, we encourage you to continuously push
        your changes to GitHub. This enables other people with access to the
        repo to see your active branch. If you have a problem (e.g. a link won't
        work), a team member(s) could access your branch and help fix the problem
        prior to the changes being merged.

How to work with a previously cloned repo
-----------------------------------------

The first step is to update your local version by typing:

```
git pull https://github.com/nih-cfde/play-with-github.git
```

Then switch to your working branch:

```
git checkout --track origin/<branch_name>
```

If you don't have an active branch, you can create a new branch (see
above).

Now you can make you changes locally, add, commit and then push those
changes.

How to make changes to branch created by someone else
------------------------------------------------------

Let's assume that Bobby (a random CFDE employee) wants to make some
changes to the theme of *play-with-github.git* repo. She creates a branch called
*newtheme*. Then she adds some files, edits some existing files and
creates a PR. You review Bobby's PR and love her work,
but wish to make a couple of small edits.

First, clone the repo to your local computer using:

```
git clone https://github.com/nih-cfde/play-with-github.git
```

Or pull her latest changes to your local repo (if previously cloned):

```
git pull https://github.com/nih-cfde/play-with-github.git
```

Then switch to the branch that Bobby has been working on. Remember that
her specific branch was called *newtheme* :

```
git checkout --track origin/newtheme
```

Now you can make changes locally, add and commit those changes. When you push your changes, they will be added to Bobby's PR.

The most basic work flow in GitHub will look something like this:
-----------------------------------------------------------------

-   Copy the GitHub repo to your local computer `git clone`
-   Edit, create and/or make changes in the newly created directory
-   Stage your changes `git add`
-   Save your changes `git commit`
-   Continue editing same file(s) after initial tracking
-   Save the recent changes enabling git to track all versions `git add` `git commit`
-   Add more changes/edits to the file but the previous version was better
-   Revert back to the last working version of the file `git log` `git revert <commit hash>`
-   Publish your local changes to GitHub `git push`

Preview website on GitHub branch
--------------------------------

!!! note "Important"
    You will require admin privileges on [readthedocs.com](www.readthedocs.com) for previewing website changes from a GitHub repo.

This tutorial applies to GitHub repos that render as websites. As
described above, your changes to the website repo must be pushed to a
new branch. Before merge, the master branch has NOT been updated yet,
and so you cannot view the changes on the actual website. 

This is stepwise guide to previewing changes to a website on a GitHub
pull request through readthedocs.

(1) Create a new branch, 'preview' on that GitHub repo

(2) Go configure readthedocs to publish that branch as well as master.
    You'll need maintainer status on the readthedocs site to do that.

(3) That will create a link that shows the latest preview branch. Now,
    whenever you push changes to that branch, it will update the
    readthedocs preview branch link. You will need to use:

```
git push origin <my branch>:preview
```

!!! note "git --force option"
    If you are pushing a really large change (like revamping a website), you may want to use `preview -f`. The `-f | --force` option tells git to make the change regardless of the history of changes.
    Since this gives the branch a new starting point for its version tracking, please *ONLY* use it on the preview branch. For most edits, you want to preserve the version history so you do NOT need this option.

(4) Click on the preview branch in readthedocs. This should take you to
    the website. Copy and paste the link on your GitHub PR.
