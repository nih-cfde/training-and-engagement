---
layout: page
title: Rendering a GitHub website locally with Jekyll
---

Rendering a GitHub website locally with Jekyll
===============================================

This tutorial includes instructions for:

- [MacOS](#global-install-of-jekyll-for-macos-using-command-line)
- [Windows OS](#install-and-website-build-for-windows-os)
- [Windows Subsystem for Linux](#rendering-a-jekyll-website-on-a-windows-10-machine-using-windows-subsystem-for-linux-wsl)

Global install of Jekyll for MacOS using command line
-----------------------------------------------------

### We will need command line tools and ruby

#### Command line tools

First, we need to install the command-line tools to be able to compile
native extensions, open a terminal and run:

    xcode-select --install

#### Ruby

Jekyll requires Ruby \> 2.5.0. macOS Catalina 10.15 comes with ruby
2.6.3. If we\'re running a previous macOS system, we\'ll have to
[install a newer version of
Ruby](https://jekyllrb.com/docs/installation/macos/#brew).

### Now install Jekyll

First, we run:

    gem install --user-install bundler jekyll

Second, get the ruby version on our computer:

    ruby -v

Then append our path file with the following, replacing the X.X with the
first two digits of our Ruby version. E.g. my `ruby -v` was 2.6.3, so I
replaced X.X.0 with 2.6.0.

    echo 'export PATH="$HOME/.gem/ruby/X.X.0/bin:$PATH"' >> ~/.bash_profile

To check that our gem paths point to our home directory, we run:

    gem env

Gem paths should look something like this:

    - GEM PATHS:
        - /Library/Ruby/Gems/2.6.0
        - /Users/abbysmith/.gem/ruby/2.6.0
        - /System/Library/Frameworks/Ruby.framework/Versions/2.6/usr/lib/ruby/gems/2.6.0

### Building the website locally:

Open Terminal.

Navigate to the publishing source for our site. *The publishing source
is the folder where the source files for your site live*

Run Jekyll site locally:

    bundle exec jekyll serve

The output looks like this:

    (base) abbys-MacBook-Pro:welcome-to-cfde abbysmith$ bundle exec jekyll serve
    Configuration file: /Users/abbysmith/Desktop/GitHub/welcome-to-cfde/_config.yml
                Source: /Users/abbysmith/Desktop/GitHub/welcome-to-cfde
           Destination: /Users/abbysmith/Desktop/GitHub/welcome-to-cfde/_site
     Incremental build: disabled. Enable with --incremental
          Generating...
            Pagination: Pagination is enabled, but I couldn't find an index.html page to use as the pagination template. Skipping pagination.
                        done in 0.952 seconds.
     Auto-regeneration: enabled for '/Users/abbysmith/Desktop/GitHub/welcome-to-cfde'
        Server address: http://127.0.0.1:4000
      Server running... press ctrl-c to stop.

Copy and paste the server address into a web browser.

Done.

Install and Website Build for Windows OS
----------------------------------------

These instructions are incomplete and further input is requested. This
version was performed on a Windows 10 machine using VS Code.

#### Setup - Basic requirements - in terminal:

    install jekyll
    install ruby
    install bundle

#### Special Windows Requirements:

    gem install tzinfo
    gem install tzinfo-data

In the Gemfile edit to include:

    gem 'tzinfo-data'

#### Build the site locally:

    bundle exec jekyll serve

 Rendering a jekyll website on a Windows 10 machine using Windows Subsystem for Linux (WSL)
-------------------------------------------------------------------------------------------

I followed some of the steps from this tutorial
(<https://connelhooley.uk/blog/2018/03/11/installing-jekyll>) but with
modifications for updated software. I used the Ubuntu terminal for this
tutorial.

1.  Install bash. Search for \'Turn Windows features on or off\' with
    the Windows start search bar. Check the \'Windows Subsystem for
    Linux\' box. Follow prompts, I had to restart my computer.
2.  Enable developer mode. Go to Settings \> Update & Security \> For
    developers. Select the \'Developer mode\' option. This may take a
    few minutes to finish.
3.  Get Ubuntu. Open the Microsoft Store and search for \'Ubuntu\'.
    Select the blue \'Get\' box. After installation is complete, select
    \'Launch\'. A Ubuntu terminal window will appear and take a few
    minutes to finish installing. Set up a username and password (just
    for using this bash window, doesn\'t have to match your computer
    logins).
4.  Set up local package database for Ubuntu.

<!-- -->

    sudo apt-get update -y

5.  Install ruby and Jekyll!

<!-- -->

    # edit based on most up-to-date ruby version
    sudo apt-get install ruby2.7 ruby2.7-dev build-essential dh-authoreconf -y
    # update ruby gems (packages)
    sudo gem update
    # install jekyll and bundler gems
    sudo gem install jekyll bundler

6.  Now, navigate to the website directory. In this case,
    <https://github.com/nih-cfde/welcome-to-cfde>. The file system is
    structured such that you would do this:

<!-- -->

    cd /mnt/c/Users/[your path to]/welcome-to-cfde/

7.  On Windows, there are some time/zone issues. They are fixed by doing
    a few things:

<!-- -->

    # Windows has timezone issues
    sudo gem install tzinfo
    sudo gem tzinfo-data
    # gem 'tzinfo-data' should be added to the Gemfile
    # this was required to get rake
    bundle install --path vendor/bundle
    # but this now creates a directory called vendor in your repo which you do NOT want to render. Thus, it must be included in the exclude list of the config file. Use a text editor to add 'vendor/bundle' to the exlude list near the end of the config file. Note: until I got this formatted correctly, I kept getting invalid time errors.
    nano _config.yml

8.  Render the website with:

<!-- -->

    bundle exec jekyll serve

9.  Copy/paste the server address to a web browser. Voila! There\'s your
    website. When you are done checking the local version, ctrl-c to
    close the server.
