## Compression Context

Context lifteimte management

::: ZSTD_createCCtx

::: ZSTD_freeCCtx

## Configure Compression

Configure any parameters you like

::: ZSTD_CCtx_setParameter

::: ZSTD_cParameter

## Dictionary Compression

Load dictionaries

::: ZSTD_CCtx_loadDictionary

## Single shot compression

Then compress

::: ZSTD_compress2

## Streaming compression

You can also use streaming compression

::: ZSTD_compressStream2

::: ZSTD_EndDirective