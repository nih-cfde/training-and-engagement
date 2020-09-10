## Tutorial recordings style guide

### Recording and video hosting tools
- Vidlets: [loom](https://www.loom.com/) (preferred)
    - for future testing: [OBS Studio](https://obsproject.com/)
- Screencasts: [asciinema](https://asciinema.org/) (preferred)
- Video editing: [DaVinci Resolve](https://www.blackmagicdesign.com/products/davinciresolve/) 
    - Computer Requirements: 8GB (16GB Recommended) RAM and a minimum 256GB memory to operate basic functions
- Hosting Vidlets: [Kaltura](https://video.ucdavis.edu/) (Currently via AggieVideo)

### Video format guide

Vidlets will have accompanying video, text translation and audio and typically will be a few hundred MB in size. We are currently using Kaltura via [AggieVideo Services](https://video.ucdavis.edu/) to host vidlets. Vidlets with or without narrator view should be used to demonstrate a concept, tasks that may require multiple windows or tasks that may require explanation of UI with/without clicks. The vidlet should be complete in itself with respect to content without relying on the associated tutorial text. Setup, installation and introductory concepts are good choices for vidlets. 

- pacing: total length should be ~15 minutes or less
- Uploading and editing videos requires an account on AggieVideo. To embed the video into a tutorial, below the video, select the "Share" tab, then the "Embed" tab, set the video size (we are using 608x402), and add the `<iframe` link to the tutorial markdown file. 
- captions: Kaltura provides automated transcription. To edit video captions, click the "Actions" button and select "+ Caption & Enrich". If a request for caption generation has already been made, click on the pencil icon to edit text, otherwise put in a request (it may take several hours to complete). 
- adding text/image overlays to vidlets
    - vidlets should have titles overlaid on the video
    - used iMovies or DaVinci Resolve for editing
    
Screencasts have a `.cast` extension and are typically in the order of tens of KB. Screencasts can be used to showcase a code chunk for a tool/ utility or some series of commands to accomplish a task. Screencasts should not be used for running long scripts or commands that have a long execute time. The screencasts should include the commands and the associated text comments for clarity. The screencasts are currently hosted in a subfolder within the associated tutorial directory. 

- screencast pacing: total length should be ~5 minutes or less
     - default parameters are no autoplay, preload or looping
     - set speed with `s` parameter which can be fractional value
     - customize the player width and length with `cols` and `rows`
     - add text overlay using `poster`
     - text apperance with `theme` and `font-size`


     