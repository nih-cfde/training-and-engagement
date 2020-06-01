# Quick Jekyll tutorial for local rendering of websites

This tutorial includes instructions for:
- [MacOS](#global-install-of-jekyll-for-macos-using-command-line)
- [Windows OS](#install-and-website-build-for-windows-os)


## Global install of Jekyll for macOS using command line

### We will need command line tools and ruby

#### Command line tools
First, we need to install the command-line tools to be able to compile native extensions, open a terminal and run:

```
xcode-select --install
```
#### Ruby
Jekyll requires Ruby > 2.5.0. macOS Catalina 10.15 comes with ruby 2.6.3. If we're running a previous macOS system, we'll have to [install a newer version of Ruby](https://jekyllrb.com/docs/installation/macos/#brew).

### Now install Jekyll

First, we run:

```
gem install --user-install bundler jekyll
```

Second, get the ruby version on our computer:

```
ruby -v
```

Then append our path file with the following, replacing the X.X with the first two digits of our Ruby version. E.g. my `ruby -v` was 2.6.3, so I replaced X.X.0 with 2.6.0.

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
    - /Users/abhijnaparigi/.gem/ruby/2.6.0
    - /System/Library/Frameworks/Ruby.framework/Versions/2.6/usr/lib/ruby/gems/2.6.0
```

### Building the website locally:

Open Terminal.

Navigate to the publishing source for our site.
*The publishing source is the folder where the source files for your site live. Visit [this page](https://help.github.com/en/github/working-with-github-pages/about-github-pages#publishing-sources-for-github-pages-sites) for more information about publishing sources.*

Run Jekyll site locally:

```
bundle exec jekyll serve
```

The output looks like this:

```
(base) Abhijnas-MacBook-Pro:welcome-to-cfde abhijnaparigi$ bundle exec jekyll serve
Configuration file: /Users/abhijnaparigi/Desktop/GitHub/welcome-to-cfde/_config.yml
            Source: /Users/abhijnaparigi/Desktop/GitHub/welcome-to-cfde
       Destination: /Users/abhijnaparigi/Desktop/GitHub/welcome-to-cfde/_site
 Incremental build: disabled. Enable with --incremental
      Generating...
        Pagination: Pagination is enabled, but I couldn't find an index.html page to use as the pagination template. Skipping pagination.
                    done in 0.952 seconds.
 Auto-regeneration: enabled for '/Users/abhijnaparigi/Desktop/GitHub/welcome-to-cfde'
    Server address: http://127.0.0.1:4000
  Server running... press ctrl-c to stop.
```

Copy and paste the server address into a web browser.

Done.


## Install and Website Build for Windows OS

These instructions are incomplete and further input is requested. This version was performed on a Windows 10 machine using VS Code.

#### Setup - Basic requirements - in terminal:

```
install jekyll
install ruby
install bundle
```

#### Special Windows Requirements:

```
gem install tzinfo
gem install tzinfo-data
```

In the Gemfile edit to include:

```
gem 'tzinfo-data'
```

#### Build the site locally:

```
bundle exec jekyll serve
```
