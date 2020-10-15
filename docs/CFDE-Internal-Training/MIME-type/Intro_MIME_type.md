---
layout: page
title: Introduction to MIME types
---

MIME types
==============

Multipurpose Internet Mail Extensions (MIME) are standards for recognizing the format of a file.
MIME types follow a certain format:
```
media-type/subtype-identifier
```

`image/png` is an example of MIME type where media-type is image and png is the subtype-identifier.

The [Internet Assigned Numbers Authority (IANA)](https://www.iana.org) is the
governing body responsible for all the [official MIME types](https://www.iana.org/assignments/media-types/media-types.xhtml).
Internet programs such as Web servers and browsers work with MIME type and not file extensions to ensure consistent transfer of same types of files irrespective of the underlying operating system. Upload of Common Fund (CF) data to CFDE portal and download by users or other CF programs is a crucial example where MIME type will help determine content of the media-types and avoid erroneous file transfers.

A few general rules for the MIME types are:

- The x- prefix of a MIME subtype-identifier implies that it is non-standard i.e. not registered with IANA.
  e.g. **Adobe Flash: application/x-shockwave-flash**

- The vnd prefix of a MIME subtype-identifier means that the MIME value is vendor specific.
  e.g. **Microsoft Excel: application/vnd.ms-excel**

- MIME type for unknown file type is generally `application/octet-stream`.
