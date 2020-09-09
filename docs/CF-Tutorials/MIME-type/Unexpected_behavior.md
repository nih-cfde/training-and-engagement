---
layout: page
title: Unexpected Behavior
---

Unexpected behavior for MIME type
==================================

So far in our tutorial we have used a proprietary file extension as an example to showcase the various options across different utilities/tools for file format identification. The results will be generally consistent across all utilities for the common file extensions, but there can be differences.

In our example files directory we have a `.vcf` file associated with variant-call-format in bioinformatics/genetics fields. However, `.vcf` is considered a file format standard for electronic business cards. This results in different results across the tools for the same input file.

|   Tool    |   Output   |
|-----------|------------|
| file      | text/plain |
| mimetype  | text/vcard |
| xdg-mime  | text/vcard |
| siegfried | text/vcard |

!!! note "Siegfried"
    Siegfried with the PRONOM default signature file correctly identified the format as 'Variant Call Format' but had no associated MIME type. Using the MIME-info signature databases results in the `text/vcard`.

In cases where there may be erroneous file extensions, it is useful to examine the file contents in addition to the file format. Assume there was mistake in renaming of a file and the `coatColor.pheno` was named `coatColor.png`
without any change in the file contents.

`file` which examines the contents of the file before reporting its type results in `text/plain`.

The `mimetype` and `xdg-mime` results in `image/png`. However, using the `--debug` flag in `mimetype` it is evident that the MIME type was extracted based on extension. This behavior can be overridden with the use of `-M` or `--magic-only` flag which only considers the contents of the file without accounting for the extensions or globs. The result is similar to `file` with `text/plain` MIME type.

Siegfried with the MIME-info database also results in `image/png` MIME type but includes a warning message, indicating signature error.

```
warning : 'match on filename only; byte/xml signatures for this format did not match'
```

Using `file` without the `--mime-type` flag reports additional information that could be useful for debugging differences in MIME types.
