# Website navigation

The website uses several forms of navigation:

Navigation type | Specified in | About
--- | --- | ---
Navigation tabs | name of tutorial folder, separated by hyphens (e.g., 'Bioinformatics-Skills/' is rendered as 'Bioinformatics Skills') **&ast;** | Currently organized by topic and user type
Left TOC | `.pages` | Table of contents for each tab and tutorial (if there are multiple pages)
Right TOC | `#`, `##`, etc. markdown headers | Table of contents for each tutorial page
Tab overview | `index.md` for tab section | Landing page describing the contents of the tab
Tutorial overview | `index.md` for tutorial section, must include a yaml header with `layout: page` and a `title:` | Landing page describing the contents of the tutorial
Next/previous page arrows | `index.md` yaml header `title:` for initial page of each section; page titles in markdown (`#`) used for subsequent pages | Navigate from one page to the next with the arrows at the bottom of each webpage

&ast;MkDocs will ignore separators like hyphens and underscores, but not periods. For consistency, only separate words with hyphens!

## Website style guide

### Formats
- the website is made with [MkDocs with Material theme](https://squidfunk.github.io/mkdocs-material/)
- the text follows [Markdown syntax](https://www.markdownguide.org/basic-syntax/)

### Website features
MkDocs extensions, plugins, and other features are specified in the `mkdocs.yml` file. Additional MkDocs features are specified in the `stylesheets/extra.css` file.

A useful resource is the [list of plugins for MkDocs categorized by functionality](https://github.com/mkdocs/mkdocs/wiki/MkDocs-Plugins).

These are the `markdown_extensions` the website uses:

Extension | About
--- | ---
`pymdownx.details` | [Details is an extension that creates collapsible elements that hide their content](https://facelessuser.github.io/pymdown-extensions/extensions/details/)
`pymdownx.emoji` | [The Emoji extension adds support for inserting emoji via simple short names enclosed within colons](https://facelessuser.github.io/pymdown-extensions/extensions/emoji/)
`admonition` | [admonition boxes](https://squidfunk.github.io/mkdocs-material/reference/admonitions/#types). Note that admonition boxes do not render in the markdown doc, you need to check the formatting in the rendered website
`codehilite` with `guess_lang: false` | [Add code/syntax highlighting to standard Python-Markdown code blocks](https://python-markdown.github.io/extensions/code_hilite/)
`toc` with `permalink: true` | [Table of contents formatting](https://www.mkdocs.org/user-guide/configuration/#formatting-options)
`pymdownx.superfences` | [Nesting of code/quote block fences and tabbed fences](https://facelessuser.github.io/pymdown-extensions/extensions/superfences/)
`pymdownx.extra` | [Improve compatibility of Python Markdown with PyMdown Extensions](https://facelessuser.github.io/pymdown-extensions/extensions/extra/)
`pymdownx.keys` | [Render keyboard keys](https://facelessuser.github.io/pymdown-extensions/extensions/keys/)

These are the plugins the website uses:

Plugin | About
--- | ---
`awesome-pages` | page organization
`git-revision-date-localized` | page last updated dates
`htmlproofer` | hyperlink checker

Check the MkDocs [changelog page](https://squidfunk.github.io/mkdocs-material/changelog/) for release, upgrade, and deprecation notes! When there are plugin updates, be sure to update your local version of the plugins. The Github repo always pulls the most up-to-date version of the plugins.

For example, these commands update `mkdocs-material` and the `awesome-pages` plugins:
```
pip install --upgrade mkdocs-material
pip install mkdocs-awesome-pages-plugin --upgrade
```

These are the Javascript features the website uses:

Feature | About
--- | ---
`javascripts/medium-zoom.min.js` and `javascripts/extra.js` files | specify image zoom functionality
markdown extensions ([`pymdownx.tabbed`](https://facelessuser.github.io/pymdown-extensions/extensions/tabbed/)) as well as custom javascript logic for the drop down boxes | specify user-chooseable boxes and dropdown boxes
`javascripts/asciinema-player.js` and `stylesheets/asciinema-player.css` files | enables integration of asciinema player for screencasts rendering

### Website host and Github branches
The website is hosted by readthedocs.com. The development branch `dev` renders to a private ReadtheDocs website (*Admin access required*): https://training.nih-cfde.org/en/dev/. The public-facing branch `stable` renders to a public ReadtheDocs website: https://training.nih-cfde.org

Guidelines for adding or fixing website content are detailed in the CFDE's training website release plan. To summarize, changes are merged from individual branches to the `dev` branch as pull requests (PR). For each public release of the website, `dev` is merged to `stable`. Each PR stage requires a form of review specified in the release plan. Release notes are generated for each new release and added to the "Release Notes" page of the website.
