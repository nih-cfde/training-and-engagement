## Pangeo binder style guide

We are using [pangeo binders](https://binder.pangeo.io/) to create pre-configured computing environments for tutorials where useful.

- Binders require their own Github repo. Create a new repo in the [`nih-cfde` organization](https://github.com/nih-cfde/). The repo name should begin with `training-<subject of binder>`. 
- The repo should have a README.md and the files you want included in the binder. 
- The README.md should state and/or link the tutorial the binder is for, the binder badge (generated on the pangeo website), and any useful information about how the binder was set up (e.g., the branch the binder was generated from). 
- After you build the binder, copy the text for the binder badge to the repo README.md and the tutorial markdown page.
