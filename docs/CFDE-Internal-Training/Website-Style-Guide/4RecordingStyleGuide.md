# Tutorial recordings style guide

## Recording and video hosting tools
- Vidlets: for speaker appearances we prefer the 'circle' picture window   
    - [loom](https://www.loom.com/) (preferred)
    - [OBS Studio](https://obsproject.com/)
- Screencasts: [asciinema](https://asciinema.org/) (preferred)
- Video editing: [DaVinci Resolve](https://www.blackmagicdesign.com/products/davinciresolve/)
    - Computer Requirements: 8GB (16GB Recommended) RAM and a minimum 256GB memory to operate basic functions.
    Though we recommend this software, most generic video editors will work fine (Kaltura has a very basic editing suite).
- Hosting Vidlets: [Kaltura](https://video.ucdavis.edu/) (Currently via AggieVideo)

## Video format guide

Vidlets will have accompanying video, text translation, and audio. They will typically be a few hundred MB in size. We are currently using Kaltura via [AggieVideo Services](https://video.ucdavis.edu/) to host vidlets (requires membership). Vidlets with or without narrator view should be used to demonstrate a concept, tasks that may require multiple windows, or tasks that may require explanation of UI with/without clicks. Short vidlets can also be used for video walkthroughs for UI. The vidlet should be complete in itself with respect to content without relying on the associated tutorial text. Setup, installation, and introductory concepts are good choices for vidlets. Please use the [powerpoint slide template](./tutorial_template_docs/Template Powerpoint.pptx), if there is a slide presentation in the vidlet.

- pacing: total length should be ~15 minutes or less
- uploading and editing videos requires an account on AggieVideo
- add user uploaded video to the Common Fund Data Ecosystem (CFDE) Training Videos channel. Use "AggieChannels" tab to search for CFDE to get the channel. Alternatively, under your username on the right top corner, click on "My Channels". Click on Filters to set Channel Membership to "Member" and search for CFDE in the search bar.
- To embed the video into a tutorial, below the video, select the "Share" tab, then the "Embed" tab, set the video size (we are using 608x402), and add the `<iframe` link to the tutorial markdown file.
- captions: Kaltura provides automated transcription. To edit video captions, click the "Actions" button and select "+ Caption & Enrich". If a request for caption generation has already been made, click on the pencil icon to edit text, otherwise put in a request (it may take several hours to complete).
- adding text/image overlays to vidlets
    - vidlets should have titles overlaid on the video
    - used iMovies or DaVinci Resolve for editing

Screencasts have a `.cast` extension and are typically in the order of tens of KB. Screencasts can be used to showcase a code chunk for a tool/utility or some series of commands to accomplish a task. Screencasts should not be used for running long scripts or commands that have a long execute time. The screencasts should include the commands and the associated text comments for clarity. The screencasts are currently hosted in a subfolder within the associated tutorial directory.

- screencast pacing: total length should be ~5 minutes or less
     - default parameters are no autoplay, preload, or looping
     - set speed with `s` parameter which can be fractional value
     - customize the player width and length with `cols` and `rows`
     - add text overlay using `poster`
     - text apperance with `theme` and `font-size`
