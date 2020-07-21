---
layout: page
title: Rendering a GitHub website locally with Jekyll
---

Rendering a GitHub website locally with Jekyll
===============================================

Jekyll is a static site generator written in Ruby. It can be easily integrated into GitHub pages to host project documentation, blogs or other relevant content. The first part of the tutorial includes instructions for downloading Jekyll on:

- [MacOS](#install-jekyll-on-macos-using-command-line)
- [Windows OS](#install-jekyll-on-windows-os)
- [Windows Subsystem for Linux](#install-jekyll-using-windows-subsystem-for-linux-wsl)

Install Jekyll on MacOS using command line
-------------------------------------------

The instructions follow the official Jekyll tutorial posted [here](https://jekyllrb.com/docs/installation/macos/#brew).

#### Command line tools

First, we install the command line utility tool for OS platform, by
running the following code in the terminal:

```
xcode-select --install
```

#### Ruby

Jekyll requires Ruby \> 2.5.0. macOS Catalina 10.15 comes with ruby
2.6.3. For older versions of macOS system, instructions for [installing newer Ruby version are available](https://jekyllrb.com/docs/installation/macos/#brew).

#### Jekyll

A software package with packaged Ruby library is called a `gem`. In addition to Jekyll, we will also install `bundler` which ensures consistency of environment for Ruby projects. To install a local version, we run:

```
gem install --user-install bundler jekyll
```

Next, obtain the Ruby version on your computer:

```
ruby -v
```

Append your path file with the following, replacing the X.X with the
first two digits of our Ruby version. E.g. my `ruby -v` was 2.6.3, so I
replaced X.X.0 with 2.6.0.

```
echo 'export PATH="$HOME/.gem/ruby/X.X.0/bin:$PATH"' >> ~/.bash_profile
```

To check that our gem paths point to our home directory, we run:

```
gem env
```

Gem paths should look something like this:

```
- GEM PATHS:
   - /Library/Ruby/Gems/2.6.0
   - /Users/abbysmith/.gem/ruby/2.6.0
   - /System/Library/Frameworks/Ruby.framework/Versions/2.6/usr/lib/ruby/gems/2.6.0
```
To check if the installation worked, run this:

```
jekyll -v
```

This should return a string similar to this : `jekyll 4.x.x`

Install Jekyll on Windows OS
------------------------------

This version was performed on a Windows 10 machine using VS Code on the terminal.

#### Install Jekyll and Ruby

```
install jekyll
install ruby
install bundle
```

#### Additional Windows requirements: Fix TimeZone issues

```
gem install tzinfo
gem install tzinfo-data
```

In the Gemfile edit to include:

```
gem 'tzinfo-data'
```

Install Jekyll using Windows Subsystem for Linux (WSL)
---------------------------------------------------------

The installation follows some of the steps from this [tutorial](https://connelhooley.uk/blog/2018/03/11/installing-jekyll) but with
modifications for updated software. A Ubuntu terminal on Windows 10 OS was used for this
tutorial.

#### Setup Ubuntu terminal

Search for `Turn Windows features on or off` on the
the Windows start search bar. Check the `Windows Subsystem for
Linux` box. Follow prompts and restart the computer.

To enable developer mode, go to Settings \> Update & Security \> For
developers. Select the `Developer mode` option. This may take a
few minutes to finish.

To get Ubuntu, open the Microsoft Store and search for `Ubuntu` (you will need to login to your Microsoft account).
Select the blue `Get` box. After installation is complete, select
`Launch`.

An Ubuntu terminal window will appear with this message: `Installing, this may take a few minutes...`. 
Set up a username and password.

!!! note
    The credentials setup are for using the bash window and need not match the computer logins.

Finally, set up local package database for Ubuntu:

```
sudo apt-get update -y
```

#### Install Ruby and Jekyll

The following commands are edited based on most up-to-date ruby version:

```
sudo apt-get install ruby2.7 ruby2.7-dev build-essential dh-autoreconf -y
```

Next, update RubyGems, the Ruby package manager

```
sudo gem update
```

Finally, install Jekyll and bundler (this step may take several minutes to complete):

```
sudo gem install jekyll bundler
```

Build Jekyll site using template from GitHub
---------------------------------------------

For the second half of the tutorial, we will use a template website hosted on GitHub that was generated using Jekyll to modify and add content. Navigate (`cd <directory name>`) to the location (e.g., Desktop) where you want to save the template directory on your computer and create a local copy of the repo:

```
git clone https://github.com/nih-cfde/Jekyll-demo.git
```

!!! note "Building Jekyll site in Windows OS using WSL"
    
      In Windows OS using the WSL, you would type this to navigate to the repo directory if you cloned it to your Desktop:

      ```
      cd /mnt/c/Users/<user name>/Desktop/Jekyll-demo/
      ```

Navigate to the newly created directory with the name of the repo which contains default files and folders that are created when building a new site with Jekyll. This is the directory structure:

```
├── 404.html
├── Gemfile
├── Gemfile.lock
├── README.md
├── _config.yml
├── _posts
│   └── 2020-07-10-welcome-to-jekyll.markdown
├── about.markdown
└── index.markdown
```

Jekyll template comes populated with generic content in form of markdown files i.e. `.md` extension.
The `_posts` folder will host all the content for the website and contains a markdown file.
The `_config.yml` file is in `yaml` format and stores attributes about the site as key value pairs.
The `Gemfile` is a Ruby file which stores all the dependencies for the Jekyll site.

#### Additional Windows requirements: Fix TimeZone issues

Ruby runs into timezone issues on Windows due to absence of native zoneinfo data. Sites generated with the Jekyll \> v3.4 will have instructions for handling the missing data added to the `Gemfile`. For older versions, the workaround is as follows:

```
sudo gem install tzinfo
sudo gem tzinfo-data
```

This line should be added to the `Gemfile`. It is already in the template Gemfile:

```
gem 'tzinfo-data'
```

To get rake, a task runner in Ruby, set up correctly in Windows OS, you may have to install gem dependencies on a location apart from the system's default. To do so use this code:

```
bundle install --path vendor/bundle
```

This will create a folder called `vendor` in your repo which you do NOT want to render. Thus, it must be added to the exclude list at the end of the `_config.yml` file.
The `vendor/bundle` has already been added to the template's `_config.yml` file. Check the file using a text editor:

```
nano _config.yml
```

!!! warning
    Until exclude list in `_config.yml` is correctly formatted, you are bound to run into invalid time errors.

#### Building the website locally

We are ready to build the site locally and the following code sets it up:

```
bundle exec jekyll serve
```

The output on MacOS looks like this:

```
(base) scanchi@MacBook-Pro Jekyll-demo % bundle exec jekyll serve
Configuration file: /Users/scanchi/Desktop/Jekyll-demo/_config.yml
            Source: /Users/scanchi/Desktop/Jekyll-demo
       Destination: /Users/scanchi/Desktop/Jekyll-demo/_site
 Incremental build: disabled. Enable with --incremental
      Generating...
       Jekyll Feed: Generating feed for posts
                    done in 0.269 seconds.
 Auto-regeneration: enabled for '/Users/scanchi/Desktop/Jekyll-demo'
    Server address: http://127.0.0.1:4000/
  Server running... press ctrl-c to stop
```

The output on Windows OS WSL should look like this (with your username filled in):

```
$ bundle exec jekyll serve
Configuration file: /mnt/c/Users/<username>/Desktop/Jekyll-demo/_config.yml
            Source: /mnt/c/Users/<username>/Desktop/Jekyll-demo
       Destination: /mnt/c/Users/<username>/Desktop/Jekyll-demo/_site
 Incremental build: disabled. Enable with --incremental
      Generating...
       Jekyll Feed: Generating feed for posts
                    done in 1.085 seconds.
                    Auto-regeneration may not work on some Windows versions.
                    Please see: https://github.com/Microsoft/BashOnWindows/issue216
                    If it does not work, please upgrade Bash On Windows or run Jekyll with --no-watch.
 Auto-regeneration: enabled for '/mnt/c/Users/<username>/Desktop/Jekyll-demo'
    Server address: http://127.0.0.1:4000/
  Server running... press ctrl-c to stop
```

Copy and paste the server address to a web browser to render the site. When you are done checking the local version, `ctrl-c` to close the server.

Site set up generates a `_site` folder that will have all the content associated with the rendered site.

Add content on Jekyll site
----------------------------

A neat functionality of the static site generator is the fast rendering. With the Jekyll site running (`bundle exec jekyll serve` step), navigate to the `_posts` folder and open the `.markdown` file in a text editor (in a separate terminal tab or outside editor, like Atom).

The top content delineated by \- is termed front matter.

```
---
layout: post
title:  "Welcome to Jekyll!"
date:   2020-07-10 13:11:00 -0700
categories: jekyll update
---
```

This is in `yaml` format and details included here are used by Jekyll to display the content on the main page along with creating the specific URL for this post.

We can make changes to this post and refresh the website to instantly view the rendered changes. The naming for posts in Jekyll follow the `YEAR-MONTH-DAY-title.MARKUP` format.

Make a copy of the existing file, modify the filename as
`2020-07-12-first-post.markdown` and open in a text editor.

Modify the front matter of the file to reflect the file name and date. You can add any text of your choice in this post. As an example, we have added some sample text along with code snippet highlighting. Refresh the website to see the rendered changes.  

```
---
layout: post
title:  "First Post"
date:   2020-07-12 13:11:00 -0700
categories: git
---
We used git to download the Jekyll site template. Git is a distributed version control system and will keep track of all changes made. Here are some commonly used commands:

{% highlight ruby %}
git add --all
git commit -m <message>
git push origin master
git pull
git status
{% endhighlight %}
```

Build local Jekyll site
-------------------------

So far we used an existing Jekyll template to make changes. We can also build the site locally. This requires installation of Ruby and Jekyll.

Run this code on the command line:

```
jekyll new newblog
```

!!! note
    The name following Jekyll is the name of the site, in this example "newblog".

Change into the newly created directory:

```
cd newblog
```

We are ready to build the site and locally serve:

```
bundle exec jekyll serve
```

Additional features and integration into GitHub are available as part of the [official Jekyll documentation](https://jekyllrb.com/docs/).  
