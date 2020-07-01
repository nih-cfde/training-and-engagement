---
layout: page
title: Working with Protected Branches
---

Working with Protected Branches
============================

Why use protected branches
--------------------------

Protected branches ensure that rules are enforced on any changes made to
that branch in a repo. A common branch protection rule is for PRs to be reviewed
by at least one other person before they get merged. Protected branches stop you
from making unauthorized changes to that branch, however you can make as
many changes as you want locally or to other branches. They just won't show up
in the protected branch until all the criteria are met and they are purposefully merged.
This type of set up is useful for repos that render websites. The branch that hosts
the live site can be protected, so that no one can incorporate changes to the
site until someone else checks them. That's how we make sure that no one
accidentally breaks the websites.

How to work with protected branches
-----------------------------------

First, copy the repo(s) to your local computer like so:

    git clone https://github.com/nih-cfde/example.git

Once you've cloned your repos (now called directories on your local
computer) and want to start work, you need to go into this directory and
work on a
[branch](https://github.com/nih-cfde/organization/blob/master/GitHubUsage.md#definitions)

The command to create a new branch is as follows (you can name your
branch whatever you like):

    git branch <name_of_new_branch>
    git checkout <name_of_new_branch>

*Type this command once you are inside the directory in which you wish
to make these changes.* The first line of code creates the branch and
the second line of code switches your current directory to the new
branch. Now if you make changes to the files, the changes will appear as
this new branch. Please note: when you change to a new branch, nothing
will physically change on your computer. But GitHub now knows that your
directory is a new branch.

Once you have made the necessary changes (or any changes at all), you
can push changes to git hub:

    git add .
    git commit -m <your_message>
    git push --set-upstream origin <name_of_new_branch>

The first line of code adds your changes and the second line of code
saves your changes. The last line pushes your changes to the main GitHub
repo as a new branch.

Until your [pull request](https://github.com/nih-cfde/organization/blob/master/GitHubUsage.md#definitions)
is [merged](https://github.com/nih-cfde/organization/blob/master/GitHubUsage.md#merging-pull-requests),
you can keep working on the same branch, and do as many pushes as you
want. You don't need to make a new branch for every change. However,
once your pull request is merged, a branch is closed, and you will have
to start a new branch (using the same branch name is ok).

### Important housekeeping notes:

-   It is important to keep the workspace clean by deleting abandoned
    branches and branches that have already been merged. In addition to
    cluttering the workshop, abandoned branches can cause collisions
    with new work.
-   As you work on your branch, we encourage you to push your changes to
    GitHub a lot. This is because other people can see your branch, so
    if you have a problem (e.g. a link won't work), someone else could
    go onto your branch and help fix the problem without it having to be
    merged in broken.

How to work with a previously cloned repo
-----------------------------------------

The first step is to update your local version by typing:

    git pull https://github.com/nih-cfde/example.git

Then switch to your working branch:

    git checkout --track origin/<branch_name>

If you don't have an active branch, you can create a new branch (see
above)

Now you can make you changes locally, add, commit and then push those
changes

How to make changes to someone else\'s branch
---------------------------------------------

Let's assume that Bobby (a random CFDE employee) wants to make some
changes to the theme of *examples.git* repo. She creates a branch called
*newtheme*. Then she adds some files, edits some existing files and
creates a pull request (PR). You review Bobby's PR and love her work,
but wish to make a couple of small edits.

First, clone the repo to your local computer using:

    git clone https://github.com/nih-cfde/example.git

Or pull her latest changes to your local repo (if previously cloned):

    git pull https://github.com/nih-cfde/example.git

Then switch to the branch that Bobby has been working on. Remember that
her specific branch was called *newtheme* :

    git checkout --track origin/newtheme

Now you can make changes locally, add and commit those changes

When you push your changes, they will be added to Bobby's PR

The most basic work flow in GitHub will look something like this:
-----------------------------------------------------------------

-   clone the GitHub repo to your local computer
-   do some work on your local git directories
-   git add your changes
-   git commit your changes
-   make some edits to the file.... you like it better now
-   git add and commit that file (now git is saving both versions)
-   make more edits to a file...decide you don't like them
-   just tell git to go back to the last version instead of clicking
    undo 50 times
-   git push all your changes

Preview website on GitHub branch
--------------------------------

You will need admin privileges on readthedocs for this!

This tutorial applies to GitHub repos that render as websites. As
described above, your changes to the website repo must be pushed to a
new branch. Before merge, the master branch has NOT been updated yet,
and so you cannot view the changes on the actual website. 

Here is a quick tutorial to previewing changes to a website on a GitHub
pull request through readthedocs.

(1) Create a new branch, 'preview' on that GitHub repo

(2) Go configure readthedocs to publish that branch as well as master.
    You'll need maintainer status on the readthedocs site to do that.

(3) That will create a link that shows the latest preview branch. Now,
    whenever you push changes to that branch, it will update the
    readthedocs preview branch link. You will need to use:

<!-- -->

    git push origin <my branch>:preview

If you are pushing a really large change (like revamping a website), you may want to use `preview -f`. The `-f | --force` option
tells git to make the change regardless of the history of changes. Since this gives the branch a new starting point
for its version tracking, please *ONLY* use it on the preview branch. For most edits, you want to preserve the version history
so you do NOT need this option.

(4) Click on the preview branch in readthedocs. This should take you to
    the website. Copy and paste the link on your GitHub PR.
