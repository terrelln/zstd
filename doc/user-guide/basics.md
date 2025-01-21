---
title: Basics
---
These are the basic helper functions that you need to interact with `zstd`.

## Error Checking

Any function that returns a `size_t` must be checked with `ZSTD_isError()`.

::: ZSTD_isError

If the function returned an error, then `ZSTD_getErrorName()` will return the error string.

::: ZSTD_getErrorName

## Compression Utilities

`ZSTD_compressBound()` returns the maximum possible compressed size of a given
input, which will be slightly larger than the input size. If you provide an
output buffer at least this large to [ZSTD_compress()][ZSTD_compress],
[ZSTD_compressCCtx()][ZSTD_compressCCtx], or [ZSTD_compress2()][ZSTD_compress2]
the compressed result is guaranteed to fit.

::: ZSTD_compressBound

## Decompression Utilities

`ZSTD_getFrameContentSize()` will return the (decompressed) content size of a
Zstandard or Skippable frame. If the content size is not known it will return
[ZSTD_CONTENTSIZE_UNKNOWN][], and if an error occurs it will return
[ZSTD_CONTENTSIZE_ERROR][].

???+ warning 

    This only returns the content size of the first compressed frame in the
    compressed source. Zstd decompression functions like
    [ZSTD_decompress()][ZSTD_decompress] will decompress zero or more
    concatenated frames. So calling this function does not guarantee that
    decompression will have enough space, unless you know a-priori that there
    is only one compressed frame in the source.

::: ZSTD_getFrameContentSize

::: ZSTD_CONTENTSIZE_UNKNOWN
    options:
        show_description: false

::: ZSTD_CONTENTSIZE_ERROR
    options:
        show_description: false

`ZSTD_findFrameCompressedSize()` will return the compressed size of the Zstd
frame. This function can be useful if:

- You don't know the compressed size of the frame.
- There may be more than one compressed frame in the compressed source.

???+ note

    This function has to traverse the entire Zstd compressed frame, and read
    each block's header. It is cheap, given that there is normally one block
    for every 128 KiB of content, but it isn't free.

::: ZSTD_findFrameCompressedSize