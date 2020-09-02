## Tutorial review and merge

The current workflow for the CFDE-CC team is to make contributions of tutorials, edits, vidlets, etc. on individual branches on the training GitHub repo. When tutorials/edits are ready for review, create a pull request (PR) to merge into `latest`. This merge requires at least two approving review. Once you have submitted the PR, the Github automation we have set up will check that your changes can be rendered without breaking the website, check for broken links, and auto-generate a preview of the rendered website.

When reviewing a PR, follow these guidelines:
- Check spelling, grammar, and overall flow of the content
- Check broken links and add suggestions for links if necessary
- Check for computing system requirements along with code syntax
- Document any errors as a result of code run through
- Make suggestions for code and content clarity 
- Check that text and images render on the website (not only the Markdown file)
- Check for functional navigation between new tutorial pages and rest of the website

If the requested changes have been made, approve the PR. After approval, the PR author should merge their changes into `latest` and delete their branch.

An issue will be posted in the training Github repo to keep track of new content for the next release of the website. After PRs are merged, add a summary of the updates in this issue. 

The [training website release plan](https://hackmd.io/O8k5wQvrQui_8_dqXOw2jA?both##release-checklist) will dictate the frequency of public releases. Releases are currently set to coincide with CFDE deliverable dates and occur in two to three month time periods. Release updates will be posted on the training website.

Merging changes into `stable` from `latest` requires group approval. Guidelines for reviewers are the same as above, and especially include:
- Check for any redundancies in new tutorial content with existing material - can the new tutorial link to an existing tutorial instead?
- Check for site wide feature integration and functionality if applicable
- Website rendering on multiple browsers in different modes
    
