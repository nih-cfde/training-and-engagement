---
name: New Release
about: Create an issue with a new release checklist
title: Release Month Year
labels: 'Training-Release'
assignees: 
   - raynamharris
   - jeremywalter
   - ACharbonneau
   - ecjoye
   - jessicalumian
   - ctb
---

# New Release Checklist

The detailed release standard operating procedure (SOP) can be found [here](https://github.com/nih-cfde/training-and-engagement/blob/dev/docs/TrainingRepoReleasePlan/TrainingRepo-Release-Plan.md). Read it carefully, repeatedly. 

Very briefly

- [ ] Create a tag for the release (e.g month-year)
- [ ] Tag all (open and closed) PRs and Issues that relate to updates and features highlighed on the [training website](https://training.nih-cfde.org/en/latest/)
- [ ] Once all PRs are merged into `dev`, create a PR to merge `dev` into `stable`
- [ ] Fix any issues that may have arisen after the last step by creating and merging PRs to `dev`
- [ ] Merge `dev` into `stable`
- [ ] Create the GitHub auto-genearted Release (which will be visible [here](https://github.com/nih-cfde/training-and-engagement/releases)
- [ ] Write the human friendly version of the release notes for the website [here](https://training.nih-cfde.org/en/latest/Release-Notes/). Create and merge. PR. 
- [ ] Merge `stable` into `dev` to resync 
- [ ] Add a comment in the [open issue of the Announcements repo](https://github.com/nih-cfde/announcements/issues) for dissemation
