# General tutorial text style guide

## Accessibility

Please ensure that our tutorials are [accessible to the blind or visually impaired](https://www.lehman.cuny.edu/academics/education/education-technology/online-accessibility/documents/LehmanCollege_BlindandVisuallyImpairedAccessibilityChecklist_002.pdf) by:

- writing descriptive text for hyperlinks. Avoid using "Click here".

- writing descriptive text for images that screenreaders can detect. If you want to check how images look, use a color blindness simulator, such as [toptal](https://www.toptal.com/designers/colorfilter/)

- when possible, avoiding color combinations or images that are not color-blind friendly. If you are creating figures in R, there is a [colorblind friendly palette](http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette).

More details on accessible formats below.

## Headers
All tutorials are written with Markdown syntax. Each page should have a title (`#`) and sub-sections (`##` or `###`). The sub-sections cannot use `#`, otherwise the right-side TOC will not render. The syntax must have a space between the `#` and the text to render as a header:
```
# Title

## Introduction

### Resources
```

Markdown syntax allows up to [6 levels of hierarchy](https://www.markdownguide.org/basic-syntax/#headings), however, if you need headings beyond level 3, this is a good indication that the tutorial page could benefit from being broken down into separate pages. Please try to only use level 2 (`##`) and 3 (`###`) headings.

## Lesson in Development banners

For website section landing pages use `>` quotes. For lesson landing pages + lesson pages use css banner syntax above the page title in the markdown document:

```
div class="banner"><span class="banner-text">Lesson in Development</span></div>
```

## Markdown tables

Tables columns are defined with `|` and rows with `-`. Markdown cell widths are set to be the same for each column when rendered by Mkdocs. The syntax looks like this:

```
col1 | col2 | col 3
--- | --- | ---
row 1 thing | row 1 thing | row 1 thing
row 2 thing | row 2 thing | row 2 thing
```

When rendered, the table looks like this:

col1 | col2 | col 3
--- | --- | ---
row 1 thing | row 1 thing | row 1 thing
row 2 thing | row 2 thing | row 2 thing

## Differentiating text

Use Markdown formatting to differenciate blocks of text, code, quotes, and hyperlinks. Several format guidelines are inspired from [DigitalOcean's technical writing guide](https://www.digitalocean.com/community/tutorials/digitalocean-s-technical-writing-guidelines#formatting).

### In-line text formatting with quotation marks

Use | Example
--- | ---
Text on webpages | look for the "Actions" column
File names | A text file called "download-links.txt" **&ast;**
Folder names | create a folder called "KF_Data"

&ast;If you refer to this specific file many times, you can put the file name in quotations for the first mention. The "Snakefile" in the Snakemake tutorial is an example where we state that the file is called "Snakefile" but thereafter simply refer to it as Snakefile.

### In-line code formatting with single backticks

Use | Example
--- | ---
Command names | `nano`
Package names | `plink`
Optional commands | `git branch -d`
Conda environment names | `base`
File names and paths | `./docs/images/` **&ast;**
File extensions | `.fastq.gz`
Ports | `:3000`
Github branches | `stable`

&ast;Use the backticks for folders and file names in relative or absolute paths. Use " " for the folder or file name itself

### Buttons to Click

When lessons ask users to click a button, please use the highlight text syntax.

The text color is set to purple website-wide. Use this syntax in markdown files, where the highlighted text here is "File Repository":
```
<span class="highlight_txt">File Repository</span>
```
and it will render as:

<span class="highlight_txt">File Repository</span>

The text highlighting colors are specified in the "extra.css" file with this block:

```
.highlight_txt {
  color: #9f0bde;
  background-color: #ededed;
  padding: .25em;
}
```

### Quotes

Use `> ` for quoted text.

```
> This is a quote!
```

renders as:
> This is a quote!

### Code blocks

Use triple backticks for lines of code that users need to execute to do the tutorial, example code syntax, and command-line terminal outputs.

The triple backtick format auto-generates a copy/paste feature in the top right corner of the gray text box. Do not include the command prompt symbol `$` in the command, as it will result in a syntax error if users copy/paste the command with the `$`.

Divide multiple lines of code into separate blocks (unless they can be executed all together if copy/pasted) or indicate that commands must be run one at a time. For multi-line code blocks, include a short description of what the command does using the `#` to comment within the code block, e.g.,:

```
# change directories
cd ~/Desktop/test

# make a new file
nano newfile.txt
```

For showing input and output code:

Code block | About | Example
--- | --- | ---
Usage | to demo code command syntax | `echo <text>`
Input | to demo actual code users should type |`echo 'hello world'`
Output | to demo the expected output results | `hello world!` in a javascript tabbed box format

Use javascript user-chooseable tabs to differentiate between input/expected output code or exercise/answer with the following syntax:

```
=== "Input"

    ```
    samtools --version
    ```

=== "Expected Output"

    ```
    samtools 1.10
    Using htslib 1.10.2
    Copyright (C) 2019 Genome Research Ltd.
    ```
```

This renders as:

=== "Input"

    ```
    samtools --version
    ```

=== "Expected Output"

    ```
    samtools 1.10
    Using htslib 1.10.2
    Copyright (C) 2019 Genome Research Ltd.
    ```
Use this structure for exercises:

```
=== "Exercise"

    Exercise question

=== "Answer"

    ```
    code answer
    ```
    You can also add regular text or images
```

For some tutorials, you may want to distinguish between code entered at the terminal vs. code entered in a text editor. The Snakemake rules are an example of this
distinction, where Snakemake rule code is added to a file and not entered at the terminal:

```
=== "Snakemake rule"

    Backticks and indentation important for rendering the entire rule in a code block

    ```
    rule rule_name:

    shell:

        # for single line commands

        # command must be enclosed in quotes

        "command"
    ```
```

We are also using this tabbed format for "Prerequisites" and "Tutorial Resources". This format helps to declutter the front matter so there are fewer admonition boxes:

```
=== "Prerequisites"

    Specify tutorial prerequisites (e.g., experience level, OS system, software installations)

=== "Tutorial Resources"

    e.g., cheat sheets, example scripts
```

### Keys syntax

Syntax for keyboard keys will follow the [Keys PyMdown extension format](https://facelessuser.github.io/pymdown-extensions/extensions/keys/). It is built around the `+` symbol. A single of combination of keys is surrounded by `++` while each key is separated by a single `+`.
List of all available key syntax are listed on the official Keys extension page.

Example syntax for `CTRL+C` would be ++ctrl+c++ which will render as keyboard keys in the site.

### Hyperlink syntax

Use hyperlinks to link images, other Markdown files, or websites. Provide short but descriptive text for the links instead of 'Click here' or 'Click this'.

Type | Syntax
--- | ---
images | `![](relative/path/to/docs/images/imagefile "short description of image")` **&ast;**
tutorial references | `[text about the reference](relative/path/to/docs/Resources/referencefile)`
other pages in the Github repo | `[text about the page](relative/path/to/page/in/Github/repo/filename)`
URLs with title | `[text about the website](URL)`
URL link itself | `<URL>`

&ast;By including a short description of the image in the quotations marks, [screen readers can provide information about the image without needing to see the image](https://blog.jwf.io/2019/06/markdown-accessible-images/). This is one way to make images accessible to the blind or vision-impaired.


Note that the image link format below results in a broken image link on Mkdocs, even though it renders on Github Markdown:

`<img src="path/to/img" title="short description of image" width="200" height="300">`

The `htmlproofer` plugin checks for broken page and URL links when a pull request (PR) is submitted to the Github repo. The results are shown in the `Checks` tab of the PR page. Under `Build documentation on PR` and under `build-and-deploy`, expand the `Build site` section to view any error messages. Common issues:

- `WARNING`: for links on pages that are broken, file paths that don't exist anymore because of file name changes
- `404 error`: page does not exist. This might be because it is a brand new file or the file was changed and neither are on the `latest` branch of the Github repo. You don't need to worry about this - it'll be solved when PR changes are merged.

### Admonition types

There are several [built-in admonition styles, and it is possible to add custom titles and styles](https://squidfunk.github.io/mkdocs-material/reference/admonitions/#types). Below is the syntax for commonly used built-in styles to highlight important information:

```
!!! info

    Important information
```

!!! info

    Important information

```
!!! tip

    Lesson tips/shortcuts/alternatives
```

!!! tip

    Lesson tips/shortcuts/alternatives

```
!!! note

    Lesson notes, learning objectives, prerequisites (see below)
```

!!! note

    Lesson notes, learning objectives, prerequisites (see below)

```
!!! warning

    Warnings (e.g., actions that result in data loss or incur costs to run)
```

!!! warning

    Warnings (e.g., actions that result in data loss or incur costs to run)

```
!!! error

    Errors users may encounter and what to do about them (e.g., software installation errors)
```

!!! error

    Errors users may encounter and what to do about them (e.g., software installation errors)

Customize admonition titles by specifying the admonition title in quotes. For example:

```
!!! note "Learning Objectives"

    Learning objectives for the lesson
```

!!! note "Learning Objectives"

    Learning objectives for the lesson

```
!!! note "Key Points"

    Lesson key points/take home messages
```

!!! note "Key Points"

    Lesson key points/take home messages

Admonition blocks do not render on Github. More examples are shown on the [supported types page](https://squidfunk.github.io/mkdocs-material/reference/admonitions/#supported-types).

### Emojis

Use this syntax in markdown docs to render emojis `:fontawesome-regular-trash-alt:` which renders as :fontawesome-regular-trash-alt:

Specified in mkdocs.yml by:
```
- pymdownx.emoji:
     emoji_index: !!python/name:materialx.emoji.twemoji
     emoji_generator: !!python/name:materialx.emoji.to_svg
```

Here is a gallery of available [fontawesome](https://fontawesome.com/icons?d=gallery&p=2&m=free) icons.
