#ifndef ZSTD_OPT2_H
#define ZSTD_OPT2_H

#include "zstd_compress_internal.h"

size_t ZSTD_compressBlock_btopt2(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        const void* src, size_t srcSize);

size_t ZSTD_compressBlock_btopt2_dictMatchState(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        const void* src, size_t srcSize);

size_t ZSTD_compressBlock_btopt2_extDict(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        const void* src, size_t srcSize);

#endif  /* ZSTD_OPT2_H */
