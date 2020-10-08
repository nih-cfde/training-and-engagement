---
layout: page
title: xdg-utils
---

xdg-utils
==========

The [xdg-utils](https://www.freedesktop.org/wiki/Software/xdg-utils/) package consists of set of tools to allow for easy integration with the desktop environment and also offers options for modifying and adding new MIME types.

<asciinema-player src="../mime_supplementary_files/xdg-mime_screencast.cast" speed="2" theme="tango" font-size="medium" cols="60" rows="15" poster="data:text/plain,\x1b[1;37mTerminal Vidlet for xdg-mime"></asciinema-player>


=== "Installation"

    ```
    sudo apt update -y
    sudo apt install -y xdg-utils
    ```

=== "Usage"

    ```
    xdg-mime query filetype <name of the file>
    ```

=== "Input"
    ```
    xdg-mime query filetype 6285633006_R03C01_Red.idat
    ```

=== "Expected Output"

    ```
    application/octet-stream
    ```

Adding custom MIME types is similar to the `mimetype` utility: we can create the xml file and update the local database. Following our previous example, we will create the `illumina-idat.xml` file:

```
<?xml version="1.0"?>
<mime-info xmlns='http://www.freedesktop.org/standards/shared-mime-info'>
  <mime-type type="application/vnd.binary">
    <comment>Illumina proprietary IDAT format</comment>
    <glob pattern="*.idat"/>
  </mime-type>
</mime-info>
```

!!! note "shared-mime-info"
    The `shared-mime-info` packages contains the core database of common MIME types and is utilized by both `xdg-mime` and `mimetype`. If run consecutively, changes made in the previous section for `mimetype` will reflect with `xdg-mime` without additional steps of adding custom types.

Update the database:

```
xdg-mime install illumina-idat.xml
```

Addition of custom MIME types in `xdg-mime` are for current user by default when called by a non-root user while the changes are system wide when called by root. This behavior is controlled by the `--mode` flag with either `user` or `system` as arguments.

!!! note "Revert to default"
    To revert back to the default MIME type use the `uninstall` flag.
    ```
    xdg-mime uninstall illumina-idat.xml
    ```
    Since `xdg-mime` and `mimetypes` rely on the shared-mime-info database for assigning MIME type, system wide changes will override local user related changes. If the uninstall command above does not revert to the default MIME type, a system wide uninstall can be run:
    ```
    sudo xdg-mime uninstall --mode system illumina-idat.xml
    ```
