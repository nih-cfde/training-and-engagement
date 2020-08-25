---
layout: page
title: Siegfried
---

#### Siegfried

Another signature-based file format identification tool is Siegfried. The current installation instructions are for **64-bit systems running Ubuntu/Debian OS**. A full list of options and installation instructions for multiple platforms can be found on the [official page](https://www.itforarchivists.com/siegfried).

=== "Installation"

    ```
    wget -qO - https://bintray.com/user/downloadSubjectPublicKey?username=bintray | sudo apt-key add -
    echo "deb http://dl.bintray.com/siegfried/debian wheezy main" | sudo tee -a /etc/apt/sources.list
    sudo apt update && sudo apt install siegfried
    ```

=== "Usage"

    ```
    sf <name of the file>
    ```

=== "Example output"

    ```
    sf 9969477031_R02C01_Red.idat

    ---
    siegfried   : 1.8.0
    scandate    : 2020-08-06T20:14:22Z
    signature   : default.sig
    created     : 2020-01-21T23:30:42+01:00
    identifiers :
      - name    : 'pronom'
        details : 'DROID_SignatureFile_V96.xml; container-signature-20200121.xml'
    ---
    filename : '9969477031_R02C01_Red.idat'
    filesize : 8095228
    modified : 2020-08-05T20:51:17Z
    errors   :
    matches  :
      - ns      : 'pronom'
        id      : 'UNKNOWN'
        format  :
        version :
        mime    :
        basis   :
        warning : 'no match'
    ```

The default results are in the [National Archives UK's PRONOM](https://www.nationalarchives.gov.uk/PRONOM/Default.aspx) file format signature which is displayed in YAML format. There are built in flags such as `-csv` or `-json` to change the output format.

Modification and customization of the underlying signature database is done using [`roy` tool](https://github.com/richardlehane/siegfried/wiki/Building-a-signature-file-with-ROY). It is installed with homebrew and Ubuntu packages. [Installation instruction](https://github.com/richardlehane/siegfried/wiki/Building-a-signature-file-with-ROY) in documentation for `roy`. To build a MIME-info signature file, we can use the included signature files from [Apache Tika](https://tika.apache.org) (tika-mimetypes.xml) and [freedesktop.org](https://www.freedesktop.org/wiki/) (freedesktop.org.xml) and use the `-mi` flag. The `roy build` creates a new signature file while `roy add` adds a new identifier to an existing signature file. The changes will be reflected to the included `default.sig` file.

!!! note "MIME-info signature file"
    [Apache Tika](https://tika.apache.org) and [freedesktop.org](https://www.freedesktop.org/wiki/) are examples of software projects that contain toolkits for content detection and metadata extraction for many commonly used file types all stored in their respective MIME-info signature files. `roy` allows the ability to use any valid MIME-info file as a source when building signatures.

=== "MIME-info database"

    ```
    # Build a MIME-info database using tika identifiers
    roy build -mi tika-mimetypes.xml

    # Add freedesktop.org MIME signature list
    roy add -mi freedesktop.org.xml
    ```

=== "MIME-info output"

    ```
    # Get MIME type
    sf 9969477031_R02C01_Red.idat

    ---
    siegfried   : 1.8.0
    scandate    : 2020-08-06T18:24:38Z
    signature   : default.sig
    created     : 2020-08-06T16:48:25Z
    identifiers :
      - name    : 'tika'
        details : 'tika-mimetypes.xml (1.23, 2019-12-02)'
      - name    : 'freedesktop.org'
        details : 'freedesktop.org.xml (1.15, 2019-10-30)'
    ---
    filename : '9969477031_R02C01_Red.idat'
    filesize : 8095228
    modified : 2020-08-04T04:23:05Z
    errors   :
    matches  :
      - ns      : 'tika'
        id      : 'UNKNOWN'
        format  :
        mime    : 'UNKNOWN'
        basis   :
        warning : 'no match'
      - ns      : 'freedesktop.org'
        id      : 'UNKNOWN'
        format  :
        mime    : 'UNKNOWN'
        basis   :
        warning : 'no match'
    ```

The output has updated `identifiers` to indicate the underlying database and the `mime` fields to reflect MIME type database.

!!! note "Aliases"
    You can use aliases tika and freedektop in the above commands instead of the full `.xml` file name. For example:
    ```
    roy build -mi tika
    ```
To add custom MIME types, we can update either of the two MIME-info files. The `xml` format remains consistent with previous examples for `mimetype` and `xdg-mime` and can be added to either `tika-mimeinfo.xml` or `freedesktop.org.xml`.

```
# Open the xml file
sudoedit /usr/share/siegfried/freedektop.org.xml

# Add the entry for .idat
  <mime-type type="application/vnd.binary">
    <comment>Illumina proprietary IDAT format</comment>
    <glob pattern="*.idat"/>
  </mime-type>
```

We can now build a MIME-info database with the updated files. Instead of overwriting the `default.sig` file, it is best practice to create a different signature file with different identifier using the `-name` flag. We can also create a single signature file with multiple identifiers. In our example, the `.idat` entry was added to the `freedesktop.org.xml` file. Since it would be useful to keep original identifiers, a custom signature file is built.

=== "Custom database"

    ```
    # Builds a custom.sig file using tika identifier
    roy build -mi tika -name tika custom.sig

    # Adds the modified freedesktop signature file to custom.sig
    roy add -mi freedesktop -name freedesktop custom.sig
    ```

=== "Updated output"

    ```
    # Check the file entry using custom.sig
    sf -sig custom.sig 9969477031_R02C01_Red.idat

    ---
    siegfried   : 1.8.0
    scandate    : 2020-08-06T20:57:10Z
    signature   : custom.sig
    created     : 2020-08-06T20:49:20Z
    identifiers :
      - name    : 'tika'
        details : 'tika-mimetypes.xml (1.23, 2019-12-02)'
      - name    : 'freedesktop'
        details : 'freedesktop.org.xml (1.15, 2019-10-30)'
    ---
    filename : '9969477031_R02C01_Red.idat'
    filesize : 8095228
    modified : 2020-08-04T04:23:05Z
    errors   :
    matches  :
      - ns      : 'tika'
        id      : 'UNKNOWN'
        format  :
        mime    : 'UNKNOWN'
        basis   :
        warning : 'no match'
      - ns      : 'freedesktop'
        id      : 'application/vnd.binary'
        format  : 'Illumina proprietary IDAT format'
        mime    : 'application/vnd.binary'
        basis   : 'extension match idat'
        warning : 'match on filename only'
    ```
