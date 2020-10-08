---
layout: page
title: file
---

file
======

The default option that requires **no installation** would be to use the `file` command. For our code example we will work with files from the [General Example Files](./Example_data_files.md).

<asciinema-player src="../mime_supplementary_files/file_screencast.cast" speed="2" theme="tango" font-size="medium" cols="60" rows="15" poster="data:text/plain,\x1b[1;37mTerminal Vidlet for file"></asciinema-player>

=== "Usage"

    ```
    file --mime-type <name of the file>
    ```

=== "Input"

    ```
    file --mime-type 6285633006_R03C01_Red.idat
    ```

=== "Expected Output"

    ```
    6285633006_R03C01_Red.idat: application/octet-stream
    ```

Adding the `-b` flag returns only the MIME type for the selected file without the filename.

=== "Input"

    ```
    file --mime-type -b 9969477031_R02C01_Red.idat
    ```

=== "Expected Output"

    ```
    application/octet-stream
    ```
