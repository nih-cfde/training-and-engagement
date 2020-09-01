---
layout: page
title: mimetype
---

mimetype
=========

The [`mimetype` utility](http://manpages.ubuntu.com/manpages/trusty/man1/mimetype.1p.html) closely follows the `file` command but uses MIME types instead of descriptions. It uses the shared MIME-info database and allows for options to add custom MIME types.

<asciinema-player src="../mime_supplementary_files/mimetype_screencast.cast" speed="2" theme="tango" font-size="medium" cols="60" rows="15" poster="data:text/plain,\x1b[1;37mTerminal Vidlet for mimetype"></asciinema-player>

=== "Installation"

    ```
    sudo apt install libfile-mimeinfo-perl
    ```

=== "Usage"

    ```
    mimetype <name of the file>
    ```

=== "Input"

    ```
    mimetype 6285633006_R03C01_Red.idat
    ```

=== "Expected Output"

    ```
    6285633006_R03C01_Red.idat: application/octet-stream
    ```

There are multiple options to customize the output. The `--describe` or `-d` option returns the file description instead of MIME type. The `-D` or `--debug` option prints the logic behind choosing the MIME type for the file.

=== "Usage"

    ```
    mimetype -D <name of the file>
    ```

=== "Input"

    ```
    mimetype -D 6285633006_R03C01_Red.idat
    ```

=== "Expected Output"

    ```
    > Data dirs are: /home/ubuntu/.local/share, /usr/local/share, /usr/share, /var/lib/snapd/desktop
    > Checking inode type
    > Checking globs for basename '6285633006_R03C01_Red.idat'
    > Checking for extension '.idat'
    > Checking globs for basename '6285633006_R03C01_Red.idat'
    > Checking for extension '.idat'
    > Value "^@" at offset 36 matches at /usr/share/mime/magic line 770
    > Failed nested rules
    > File exists, trying default method
    > First 10 bytes of the file contain control chars
    6285633006_R03C01_Red.idat: application/octet-stream
    ```

`mimetype` allows for addition of custom MIME types.

From our example files list, we have an `.idat` extension which is Illumina BeadArray data associated with microarray technology and solution for DNA and RNA analysis. `.idat` is an example of a file extension not listed in the [official MIME types](https://www.iana.org/assignments/media-types/media-types.xhtml). The IDAT file format varies from binary to encrypted XML depending on the array platform.

!!! warning
    For demonstration purposes we chose to apply `application/vnd.binary` as MIME type for `.idat` files. This may or may not be the appropriate choice. **The default `application/octet-stream` which is defined as arbitrary binary data is desired behavior for files without an accepted designation**. It is preferable to have unknown MIME type to prevent a file from being pushed somewhere inappropriate compared to a custom MIME type that may cause a file to be shunted into an application/pipeline where it doesn't belong.

We can add custom MIME types by creating an xml file `illumina-idat.xml` using any text editor (nano, vim, emacs etc) for the file extensions.

```
<?xml version="1.0"?>
<mime-info xmlns='http://www.freedesktop.org/standards/shared-mime-info'>
  <mime-type type="application/vnd.binary">
    <comment>Illumina proprietary IDAT format</comment>
    <glob pattern="*.idat"/>
  </mime-type>
</mime-info>
```

Three lines in the above xml code require editing for every new custom MIME type created: The `mime-type type=` sets the custom MIME type, the `comment` explains the MIME type and the `glob pattern` indicates the search pattern to be associated with the chosen MIME type.

!!! note "globs"
    The `glob` module finds all the pathnames matching a specified pattern according to the rules set by the Unix shell. For MIME type, the file extension acts as the pattern for the search across the filesystem.

Copy this file and update the mime-database.

```
sudo cp illumina-idat.xml /usr/share/mime/packages
sudo update-mime-database /usr/share/mime
```

!!! note "sudo"
    On many sudo-based distributions, it is not possible to log in as a root user. Issuing a command with `sudo` which operates on a per-command basis, helps gain administrative privileges. It is a solution to privilege-related errors. For example when copying files into directories without `w` permission for users or installation of software.

Rerunning the code now results in the updated MIME type.

=== "Input"

    ```
    mimetype 6285633006_R03C01_Red.idat
    ```

=== "Expected Output"

    ```
    6285633006_R03C01_Red.idat: application/vnd.binary
    ```

!!! note "Revert to default"
    If you decide to use the default MIME type instead of the custom association, you can delete the `.xml` file and update the mime-database.
    ```
    sudo rm /usr/share/mime/packages/illumina-idat.xml
    sudo update-mime-database /usr/share/mime
    ```
