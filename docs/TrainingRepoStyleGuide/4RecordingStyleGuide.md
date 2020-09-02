## Tutorial recordings style guide

### Recording and video hosting tools
- Vidlets: [loom](https://www.loom.com/), [OBS Studio](https://obsproject.com/)
- Screencasts: [asciinema](https://asciinema.org/)
- Video editing: [DaVinci Resolve](https://www.blackmagicdesign.com/products/davinciresolve/) 
    - Computer Requirements: 8GB (16GB Recommended) RAM and a minimum 256GB memory to operate basic functions
- Hosting Vidlets: [Kaltura](https://video.ucdavis.edu/) (Tentatively Aggie Video)

### Video format guide
- vidlet pacing: total length should be ~15 minutes or less
    - captions: Kaltura provides automated transcription 
    - adding text/image overlays to vidlets
        - vidlets should have titles overlaid on the video
        - used iMovies for editing
- screencast pacing: total length should be ~5 minutes or less
     - default parameters are no autoplay, preload or looping
     - set speed with `s` parameter which can be fractional value
     - customize the player width and length with `cols` and `rows`
     - add text overlay using `poster`
     - text apperance with `theme` and `font-size`

Screencasts have a `.cast` extension and are typically in the order of tens of KB. Screencasts can be used to showcase a code chunk for a tool/ utility or some series of commands to accomplish a task. Screencasts should not be used for running long scripts or commands that have a long execute time. The screencasts should include the commands and the associated text comments for clarity. The screencasts are currently hosted in a subfolder within the associated tutorial directory. 

Vidlets with or without narrator view should be used to demonstrate a concept, tasks that may require multiple windows or tasks that may require explanation of UI with/without clicks. The vidlet should be complete in itself with respect to content without relying on the associated tutorial text. Setup, installation and introductory concepts are good choices for vidlets. 
     
