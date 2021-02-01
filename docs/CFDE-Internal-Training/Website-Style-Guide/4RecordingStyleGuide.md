# Tutorial recordings style guide

## Recording and video hosting tools
- Vidlets: for speaker appearances we prefer the 'circle' picture window   
    - [loom](https://www.loom.com/) (preferred)
    - [OBS Studio](https://obsproject.com/)
- Screencasts: [asciinema](https://asciinema.org/) (preferred)
- Video editing: [DaVinci Resolve](https://www.blackmagicdesign.com/products/davinciresolve/)
    - Computer Requirements: 8GB (16GB Recommended) RAM and a minimum 256GB memory to operate basic functions. Though we recommend this software, most generic video editors will work fine (Kaltura has a very basic editing suite).
- Hosting Vidlets: [Kaltura](https://video.ucdavis.edu/) (Currently via AggieVideo)

## Vidlet format guide

Vidlets will have accompanying video, text translation, and audio. We are using Kaltura via [AggieVideo Services](https://video.ucdavis.edu/) to host vidlets (requires membership). They will typically be a few hundred MB in size. Vidlets with or without narrator view should be used to demonstrate a concept, tasks that may require multiple windows, or tasks that may require explanation of UI with/without clicks. Short vidlets can also be used for video walkthroughs for UI. The vidlet should be complete in itself with respect to content without relying on the associated tutorial text. Setup, installation, and introductory concepts are good choices for vidlets. Please use the "Template Powerpoint" [slide template](https://drive.google.com/drive/u/0/folders/14dOaf7-G4k7rCw5mL2Q5jdRWXrO0Y5i-), if there is a slide presentation in the vidlet (we provide a powerpoint and google slideshow version of the file).

**For Vidlets with narrator view,** be sure to include an adequate light source pointing towards the front of the narrator to make them visible to the audience. The camera should be at eye level with the narrator. It is also important the narrator have a simple solid background (e.g. wall, curtain) or a neat and tidy office space with minimal distractions.

- pacing: total length should be ~15 minutes or less
- adding text/image overlays to vidlets
    - vidlets should have
        - titles overlaid on the video ([use intro title slide](https://drive.google.com/drive/u/0/folders/14dOaf7-G4k7rCw5mL2Q5jdRWXrO0Y5i-)) at the start of the video
        - CFDE Logo watermark in the bottom left side of the screen

### Editing vidlets

- use iMovies or DaVinci Resolve for editing
- captions: Kaltura provides automated transcription. To edit video captions, click the "Actions" button and select "+ Caption & Enrich". If a request for caption generation has already been made, click on the pencil icon to edit text, otherwise put in a request (it may take several hours to complete).

### Saving vidlets

- uploading and editing videos requires an account on AggieVideo
- add user uploaded video to the "Common Fund Data Ecosystem (CFDE) Training Videos" channel. Use "AggieChannels" tab to search for "CFDE" to find the channel. Alternatively, under your username on the right top corner, click on "My Channels". Click on Filters to set Channel Membership to "Member" and search for "CFDE" in the search bar.
- to embed the video into a tutorial, select the "Share" tab below the video, then the "Embed" tab, set the video size (we are using 608x402), and add the `<iframe` link to the tutorial markdown file. The syntax should look like this:

```
<iframe id="kaltura_player" src="https://cdnapisec.kaltura.com/p/1770401/sp/177040100/embedIframeJs/uiconf_id/29032722/partner_id/1770401?iframeembed=true&playerId=kaltura_player&entry_id=1_yjudzmwr&flashvars[mediaProtocol]=rtmp&amp;flashvars[streamerType]=rtmp&amp;flashvars[streamerUrl]=rtmp://www.kaltura.com:1935&amp;flashvars[rtmpFlavors]=1&amp;flashvars[localizationCode]=en&amp;flashvars[leadWithHTML5]=true&amp;flashvars[sideBarContainer.plugin]=true&amp;flashvars[sideBarContainer.position]=left&amp;flashvars[sideBarContainer.clickToClose]=true&amp;flashvars[chapters.plugin]=true&amp;flashvars[chapters.layout]=vertical&amp;flashvars[chapters.thumbnailRotator]=false&amp;flashvars[streamSelector.plugin]=true&amp;flashvars[EmbedPlayer.SpinnerTarget]=videoHolder&amp;flashvars[dualScreen.plugin]=true&amp;flashvars[Kaltura.addCrossoriginToIframe]=true&amp;&wid=1_33bt0n0m" width="608" height="402" allowfullscreen webkitallowfullscreen mozAllowFullScreen allow="autoplay *; fullscreen *; encrypted-media *" sandbox="allow-forms allow-same-origin allow-scripts allow-top-navigation allow-pointer-lock allow-popups allow-modals allow-orientation-lock allow-popups-to-escape-sandbox allow-presentation allow-top-navigation-by-user-activation" frameborder="0" title="Kaltura Player"></iframe>
```

## Screencast format guide

Screencasts have a `.cast` extension and are typically in the order of tens of KB. Screencasts can be used to showcase a code chunk for a tool/utility or some series of commands to accomplish a task. Screencasts should not be used for running long scripts or commands that have a long execute time. The screencasts should include the commands and the associated text comments for clarity. The screencasts are currently hosted in a subfolder within the associated tutorial directory.

- screencast pacing: total length should be ~5 minutes or less
     - default parameters are no autoplay, preload, or looping
     - set speed with `s` parameter which can be fractional value
     - customize the player width and length with `cols` and `rows`
     - add text overlay using `poster`
     - text apperance with `theme` and `font-size`

### Editing screencasts

1. Download the screencast from your asciinema account (".cast" file)

2. Open in text editor (e.g. atom)

3. The file looks like this (each line has one typed character in the screencast):

     ```
     {"version": 2, "width": 75, "height": 18, "timestamp": 1598395763, "env": {"SHELL": "/bin/bash", "TERM": "xterm-256color"}}
     [0.01083, "o", "\u001b[?1034hbash-3.2$ "]
     [1.548746, "o", "\u001b[H\u001b[2Jbash-3.2$ "]
     [3.362292, "o", "#"]
     [3.849687, "o", " "]
     [4.105164, "o", "F"]
     [4.25718, "o", "i"]
     [4.328995, "o", "r"]
     [4.434286, "o", "s"]
     [4.61828, "o", "t"]
     [4.790068, "o", " "]
     [5.031393, "o", "l"]
     [5.125048, "o", "e"]
     [5.22469, "o", "t"]
     [5.396713, "o", "'"]
     [5.556953, "o", "s"]
     [5.755022, "o", " "]
     [5.932595, "o", "c"]
     [6.061395, "o", "h"]
     [6.14247, "o", "e"]
     [6.27824, "o", "c"]
     [6.402957, "o", "k"]
     [6.54971, "o", " "]
     [6.671492, "o", "t"]
     [6.794996, "o", "o"]
     [6.880366, "o", " "]
     [7.037239, "o", "m"]
     [7.132153, "o", "a"]
     [7.254111, "o", "k"]
     [7.377614, "o", "e"]
     [7.490531, "o", " "]
     [7.607784, "o", "s"]
     [7.762003, "o", "u"]
     [7.863834, "o", "r"]
     [7.926838, "o", "e"]
     ...
     ```

4. You can edit the words, capitalization, spacing, etc. that appear on the screencast. It is best to avoid edits to timestamps as that can be manipulated using the speed option of the player.

5. Save your changes and view those changes as you go by replaying on your terminal using this command:

     ```
     asciinema play /path/to/screencast/2_screencast.cast
     ```

### Saving screencasts

- for a given tutorial, save the ".cast" file(s) in a subfolder of the tutorial folder
- to embed the screencast into a tutorial:

     - change the `src` path tell the asciinema-player where to locate the ".cast" file. The relative path needs to start with `../`; `./` does not work.

     - change the `poster` screencast title

For example, if your screencast is located in `./docs/Bioinformatics-Skills/Kids-First/vidlets/example.cast`, the markdown page you're embedding the screencast in is located in `./docs/Bioinformatics-Skills/Kids-First/example.md`, and the title is "My screencast", the syntax would be:

```
<asciinema-player src="../vidlets/example.cast" speed="2" theme="tango" font-size="medium" cols="60" rows="15" poster="data:text/plain,\x1b[1;37mMy screencast"></asciinema-player>
```
