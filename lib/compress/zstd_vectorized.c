/*
 * Copyright (c) 2016-present, Facebook, Inc.
 * All rights reserved.
 *
 * This source code is licensed under both the BSD-style license (found in the
 * LICENSE file in the root directory of this source tree) and the GPLv2 (found
 * in the COPYING file in the root directory of this source tree).
 * You may select, at your option, one of the above-listed licenses.
 */

#include "zstd_compress_internal.h"
#include "zstd_fast.h"

void ZSTD_fillWideTable(ZSTD_matchState_t* ms,
                        void const* end, ZSTD_dictTableLoadMethod_e dtlm);

size_t ZSTD_compressBlock_vectorized_generic(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize)
{
    // Load the 8 hashes
    // Load the 8 positions and 8 copies of the current position (8 bytes)
    // for each position
    //   load the next 8 hashes
    //   load the next 8 positions
    //   compute the matches
    //   pick the longest one
    //   prepare for the next loop
}

size_t ZSTD_compressBlock_vectorized(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize);

// Fall back to extDict when this match finder is used
size_t ZSTD_compressBlock_vectorized_dictMatchState(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize);

size_t ZSTD_compressBlock_vectorized_extDict(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize);