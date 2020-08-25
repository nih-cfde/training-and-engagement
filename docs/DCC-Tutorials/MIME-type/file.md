---
layout: page
title: file
---
<script src="../../javascripts/asciinema-player.js"></script>

file
======

The default option that requires **no installation** would be to use the `file` command.
For our code example we will work with 6285633006_R03C01_Red.idat file from the
[General Example Files](./Example_data_files.md).

<asciinema-player src="./mime_supplementary_files/file_screencast.cast" speed="2" theme="asciinema" font-size="medium" ></asciinema-player>

=== "Usage"

    ```
    file --mime-type <name of the file>
    ```

=== "Example output"

    ```
    file --mime-type 9969477031_R02C01_Red.idat
    9969477031_R02C01_Red.idat: application/octet-stream
    ```

Adding the `-b` flag returns only the MIME type for the selected file without the filename.

```
file --mime-type -b 9969477031_R02C01_Red.idat
application/octet-stream
```
