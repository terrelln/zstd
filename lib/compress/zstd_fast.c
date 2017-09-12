/*
 * Copyright (c) 2016-present, Yann Collet, Facebook, Inc.
 * All rights reserved.
 *
 * This source code is licensed under both the BSD-style license (found in the
 * LICENSE file in the root directory of this source tree) and the GPLv2 (found
 * in the COPYING file in the root directory of this source tree).
 * You may select, at your option, one of the above-listed licenses.
 */

#include "zstd_fast.h"


void ZSTD_fillHashTable (ZSTD_CCtx* zc, const void* end, const U32 mls)
{
    U32* const hashTable = zc->hashTable;
    U32  const hBits = zc->appliedParams.cParams.hashLog;
    const BYTE* const base = zc->base;
    const BYTE* ip = base + zc->nextToUpdate;
    const BYTE* const iend = ((const BYTE*)end) - HASH_READ_SIZE;
    const size_t fastHashFillStep = 3;

    while(ip <= iend) {
        hashTable[ZSTD_hashPtr(ip, hBits, mls)] = (U32)(ip - base);
        ip += fastHashFillStep;
    }
}


FORCE_INLINE_TEMPLATE
size_t ZSTD_compressBlock_fast_generic(ZSTD_CCtx* cctx,
                                 const void* src, size_t srcSize,
                                 const U32 mls, U32 const extDict)
{
    U32* const hashTable = cctx->hashTable;
    U32  const hBits = cctx->appliedParams.cParams.hashLog;
    seqStore_t* seqStorePtr = &(cctx->seqStore);
    const BYTE* const istart = (const BYTE*)src;
    const BYTE* ip = istart;
    const BYTE* anchor = istart;
    const BYTE* const iend = istart + srcSize;
    const BYTE* const ilimit = iend - HASH_READ_SIZE;
    U32 offset_1=seqStorePtr->rep[0], offset_2=seqStorePtr->rep[1];
    U32 offsetSaved = 0;
    ZSTD_DEFINE_LIMITS(cctx, extDict);

    /* init */
    ip += (!extDict && ip == lowPrefixPtr);
    ZSTD_initOffsets(&offset_1, &offset_2, &offsetSaved, ip, extDict);

    while (ip < ilimit) {   /* < instead of <=, because repcode check at (ip+1) */
        size_t mLength;
        size_t const h = ZSTD_hashPtr(ip, hBits, mls);
        U32 const current = (U32)(ip-base);
        U32 const matchIndex = hashTable[h];
        /* Construct the match, and the lower and upper limits */
        ZSTD_DEFINE_MATCH(mstart, match, mend, matchIndex, extDict);
        /* Construct the repcode index, pointer, and the highest it can go */
        ZSTD_DEFINE_REP(repIndex, rep, rend, ip + 1, offset_1, extDict);
        hashTable[h] = current;   /* update hash table */

        /* Check if the repcode is valid and at least a 4 byte match */
        if (ZSTD_isRepValid32(repIndex, offset_1, rep, ip + 1, extDict)) {
            mLength = ZSTD_count_generic(ip + 1 + 4, rep + 4, iend, rend, extDict) + 4;
            ip++;
            ZSTD_storeSeq(seqStorePtr, ip-anchor, anchor, 0, mLength-MINMATCH);
        } else {
            U32 offset;
            if ((matchIndex <= lowestIndex) || (MEM_read32(match) != MEM_read32(ip))) {
                ip += ((ip - anchor) >> g_searchStrength) + 1;
                continue;
            }
            mLength = ZSTD_count_generic(ip + 4, match + 4, iend, mend, extDict) + 4;
            while (((ip > anchor) & (match > mstart)) && (ip[-1] == match[-1])) { ip--; match--; mLength++; }

            offset = current - matchIndex;
            offset_2 = offset_1;
            offset_1 = offset;

            ZSTD_storeSeq(seqStorePtr, ip-anchor, anchor, offset + ZSTD_REP_MOVE, mLength-MINMATCH);
        }

        /* match found */
        ip += mLength;
        anchor = ip;

        if (ip <= ilimit) {
            /* Fill Table */
            hashTable[ZSTD_hashPtr(base+current+2, hBits, mls)] = current+2;  /* here because current+2 could be > iend-8 */
            hashTable[ZSTD_hashPtr(ip-2, hBits, mls)] = (U32)(ip-2-base);
            /* check immediate repcode */
            while (ip <= ilimit) {
                U32 const current2 = (U32)(ip-base);
                /* Construct the repcode index, pointer, and the upper limit */
                ZSTD_DEFINE_REP(rep2Index, rep2, r2end, ip, offset_2, extDict);
                /* Check if the repcode is valid and at least a 4 byte match */
                if (!ZSTD_isRepValid32(rep2Index, offset_2, rep2, ip, extDict))
                    break;
                {
                    /* store sequence */
                    size_t const rLength = ZSTD_count_generic(ip + 4, rep2 + 4, iend, r2end, extDict) + 4;
                    { U32 const tmpOff = offset_2; offset_2 = offset_1; offset_1 = tmpOff; }  /* swap offset_2 <=> offset_1 */
                    hashTable[ZSTD_hashPtr(ip, hBits, mls)] = current2;
                    ZSTD_storeSeq(seqStorePtr, 0, anchor, 0, rLength-MINMATCH);
                    ip += rLength;
                    anchor = ip;
                }
    }   }   }

    /* save reps for next block */
    ZSTD_saveOffsets(seqStorePtr, offset_1, offset_2, offsetSaved);

    /* Return the last literals size */
    return iend - anchor;
}


size_t ZSTD_compressBlock_fast(ZSTD_CCtx* ctx,
                       const void* src, size_t srcSize)
{
    const U32 mls = ctx->appliedParams.cParams.searchLength;
    switch(mls)
    {
    default: /* includes case 3 */
    case 4 :
        return ZSTD_compressBlock_fast_generic(ctx, src, srcSize, 4, 0);
    case 5 :
        return ZSTD_compressBlock_fast_generic(ctx, src, srcSize, 5, 0);
    case 6 :
        return ZSTD_compressBlock_fast_generic(ctx, src, srcSize, 6, 0);
    case 7 :
        return ZSTD_compressBlock_fast_generic(ctx, src, srcSize, 7, 0);
    }
}


size_t ZSTD_compressBlock_fast_extDict(ZSTD_CCtx* ctx,
                         const void* src, size_t srcSize)
{
    U32 const mls = ctx->appliedParams.cParams.searchLength;
    switch(mls)
    {
    default: /* includes case 3 */
    case 4 :
        return ZSTD_compressBlock_fast_generic(ctx, src, srcSize, 4, 1);
    case 5 :
        return ZSTD_compressBlock_fast_generic(ctx, src, srcSize, 5, 1);
    case 6 :
        return ZSTD_compressBlock_fast_generic(ctx, src, srcSize, 6, 1);
    case 7 :
        return ZSTD_compressBlock_fast_generic(ctx, src, srcSize, 7, 1);
    }
}
