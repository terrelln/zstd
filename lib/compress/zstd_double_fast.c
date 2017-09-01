/*
 * Copyright (c) 2016-present, Yann Collet, Facebook, Inc.
 * All rights reserved.
 *
 * This source code is licensed under both the BSD-style license (found in the
 * LICENSE file in the root directory of this source tree) and the GPLv2 (found
 * in the COPYING file in the root directory of this source tree).
 * You may select, at your option, one of the above-listed licenses.
 */

#include "zstd_double_fast.h"


void ZSTD_fillDoubleHashTable(ZSTD_CCtx* cctx, const void* end, const U32 mls)
{
    U32* const hashLarge = cctx->hashTable;
    U32  const hBitsL = cctx->appliedParams.cParams.hashLog;
    U32* const hashSmall = cctx->chainTable;
    U32  const hBitsS = cctx->appliedParams.cParams.chainLog;
    const BYTE* const base = cctx->base;
    const BYTE* ip = base + cctx->nextToUpdate;
    const BYTE* const iend = ((const BYTE*)end) - HASH_READ_SIZE;
    const size_t fastHashFillStep = 3;

    while(ip <= iend) {
        hashSmall[ZSTD_hashPtr(ip, hBitsS, mls)] = (U32)(ip - base);
        hashLarge[ZSTD_hashPtr(ip, hBitsL, 8)] = (U32)(ip - base);
        ip += fastHashFillStep;
    }
}


FORCE_INLINE_TEMPLATE
size_t ZSTD_compressBlock_doubleFast_generic(ZSTD_CCtx* cctx,
                                 const void* src, size_t srcSize,
                                 U32 const mls, U32 const extDict)
{
    U32* const hashLong = cctx->hashTable;
    const U32 hBitsL = cctx->appliedParams.cParams.hashLog;
    U32* const hashSmall = cctx->chainTable;
    const U32 hBitsS = cctx->appliedParams.cParams.chainLog;
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

    /* Main Search Loop */
    while (ip < ilimit) {   /* < instead of <=, because repcode check at (ip+1) */
        size_t mLength;
        size_t const h2 = ZSTD_hashPtr(ip, hBitsL, 8);
        size_t const h = ZSTD_hashPtr(ip, hBitsS, mls);
        U32 const current = (U32)(ip-base);
        U32 const matchIndexL = hashLong[h2];
        U32 const matchIndexS = hashSmall[h];
        ZSTD_DEFINE_REP(repIndex, rep, rend, ip + 1, offset_1, extDict);
        ZSTD_DEFINE_MATCH(mstartL, matchLong, mendL, matchIndexL, extDict);
        ZSTD_DEFINE_MATCH(mstartS, match, mendS, matchIndexS, extDict);
        hashLong[h2] = hashSmall[h] = current;   /* update hash tables */

        assert(offset_1 <= current);   /* supposed guaranteed by construction */
        if (ZSTD_isRepValid32(repIndex, offset_1, rep, ip + 1, extDict)) {
            /* favor repcode */
            mLength = ZSTD_count_generic(ip + 1 + 4, rep + 4, iend, rend, extDict) + 4;
            ip++;
            ZSTD_storeSeq(seqStorePtr, ip-anchor, anchor, 0, mLength-MINMATCH);
        } else {
            U32 offset;
            if ((matchIndexL > lowestIndex) && (MEM_read64(matchLong) == MEM_read64(ip)) ) {
                mLength = ZSTD_count_generic(ip + 8, matchLong + 8, iend, mendL, extDict) + 8;
                offset = current - matchIndexL;
                while (((ip>anchor) & (matchLong>mstartL)) && (ip[-1] == matchLong[-1])) { ip--; matchLong--; mLength++; } /* catch up */
            } else if ((matchIndexS > lowestIndex) && (MEM_read32(match) == MEM_read32(ip)) ) {
                size_t const hl3 = ZSTD_hashPtr(ip+1, hBitsL, 8);
                U32 const matchIndexL3 = hashLong[hl3];
                ZSTD_DEFINE_MATCH(mstartL3, matchL3, mendL3, matchIndexL3, extDict);
                hashLong[hl3] = current + 1;
                if ((matchIndexL3 > lowestIndex) && (MEM_read64(matchL3) == MEM_read64(ip+1)) ) {
                    mLength = ZSTD_count_generic(ip + 1 + 8, matchL3 + 8, iend, mendL3, extDict) + 8;
                    ip++;
                    offset = current + 1 - matchIndexL3;
                    while (((ip>anchor) & (matchL3>mstartL3)) && (ip[-1] == matchL3[-1])) { ip--; matchL3--; mLength++; } /* catch up */
                } else {
                    mLength = ZSTD_count_generic(ip + 4, match + 4, iend, mendS, extDict) + 4;
                    offset = current - matchIndexS;
                    while (((ip>anchor) & (match>mstartS)) && (ip[-1] == match[-1])) { ip--; match--; mLength++; } /* catch up */
                }
            } else {
                ip += ((ip-anchor) >> g_searchStrength) + 1;
                continue;
            }

            offset_2 = offset_1;
            offset_1 = offset;

            ZSTD_storeSeq(seqStorePtr, ip-anchor, anchor, offset + ZSTD_REP_MOVE, mLength-MINMATCH);
        }

        /* match found */
        ip += mLength;
        anchor = ip;

        if (ip <= ilimit) {
            /* Fill Table */
            hashLong[ZSTD_hashPtr(base+current+2, hBitsL, 8)] =
                hashSmall[ZSTD_hashPtr(base+current+2, hBitsS, mls)] = current+2;  /* here because current+2 could be > iend-8 */
            hashLong[ZSTD_hashPtr(ip-2, hBitsL, 8)] =
                hashSmall[ZSTD_hashPtr(ip-2, hBitsS, mls)] = (U32)(ip-2-base);

            /* check immediate repcode */
            while (ip <= ilimit) {
                ZSTD_DEFINE_REP(rep2Index, rep2, r2end, ip, offset_2, extDict);
                U32 const current2 = (U32)(ip-base);
                if (!ZSTD_isRepValid32(rep2Index, offset_2, rep2, ip, extDict))
                    break;
                {
                    /* store sequence */
                    size_t const rLength = ZSTD_count_generic(ip + 4, rep2 + 4, iend, r2end, extDict) + 4;
                    { U32 const tmpOff = offset_2; offset_2 = offset_1; offset_1 = tmpOff; } /* swap offset_2 <=> offset_1 */
                    hashSmall[ZSTD_hashPtr(ip, hBitsS, mls)] = current2;
                    hashLong[ZSTD_hashPtr(ip, hBitsL, 8)] = current2;
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

static size_t ZSTD_compressBlock_doubleFast_4(ZSTD_CCtx* ctx, const void* src, size_t srcSize)
{
    return ZSTD_compressBlock_doubleFast_generic(ctx, src, srcSize, 4, 0);
}

static size_t ZSTD_compressBlock_doubleFast_5(ZSTD_CCtx* ctx, const void* src, size_t srcSize)
{
    return ZSTD_compressBlock_doubleFast_generic(ctx, src, srcSize, 5, 0);
}

static size_t ZSTD_compressBlock_doubleFast_6(ZSTD_CCtx* ctx, const void* src, size_t srcSize)
{
    return ZSTD_compressBlock_doubleFast_generic(ctx, src, srcSize, 6, 0);
}

static size_t ZSTD_compressBlock_doubleFast_7(ZSTD_CCtx* ctx, const void* src, size_t srcSize)
{
    return ZSTD_compressBlock_doubleFast_generic(ctx, src, srcSize, 7, 0);
}


size_t ZSTD_compressBlock_doubleFast(ZSTD_CCtx* ctx, const void* src, size_t srcSize)
{
    const U32 mls = ctx->appliedParams.cParams.searchLength;
    switch(mls)
    {
    default: /* includes case 3 */
    case 4 :
        return ZSTD_compressBlock_doubleFast_4(ctx, src, srcSize);
    case 5 :
        return ZSTD_compressBlock_doubleFast_5(ctx, src, srcSize);
    case 6 :
        return ZSTD_compressBlock_doubleFast_6(ctx, src, srcSize);
    case 7 :
        return ZSTD_compressBlock_doubleFast_7(ctx, src, srcSize);
    }
}


static size_t ZSTD_compressBlock_doubleFast_extDict_4(ZSTD_CCtx* ctx, const void* src, size_t srcSize)
{
    return ZSTD_compressBlock_doubleFast_generic(ctx, src, srcSize, 4, 1);
}

static size_t ZSTD_compressBlock_doubleFast_extDict_5(ZSTD_CCtx* ctx, const void* src, size_t srcSize)
{
    return ZSTD_compressBlock_doubleFast_generic(ctx, src, srcSize, 5, 1);
}

static size_t ZSTD_compressBlock_doubleFast_extDict_6(ZSTD_CCtx* ctx, const void* src, size_t srcSize)
{
    return ZSTD_compressBlock_doubleFast_generic(ctx, src, srcSize, 6, 1);
}

static size_t ZSTD_compressBlock_doubleFast_extDict_7(ZSTD_CCtx* ctx, const void* src, size_t srcSize)
{
    return ZSTD_compressBlock_doubleFast_generic(ctx, src, srcSize, 7, 1);
}

size_t ZSTD_compressBlock_doubleFast_extDict(ZSTD_CCtx* ctx,
                         const void* src, size_t srcSize)
{
    U32 const mls = ctx->appliedParams.cParams.searchLength;
    switch(mls)
    {
    default: /* includes case 3 */
    case 4 :
        return ZSTD_compressBlock_doubleFast_extDict_4(ctx, src, srcSize);
    case 5 :
        return ZSTD_compressBlock_doubleFast_extDict_5(ctx, src, srcSize);
    case 6 :
        return ZSTD_compressBlock_doubleFast_extDict_6(ctx, src, srcSize);
    case 7 :
        return ZSTD_compressBlock_doubleFast_extDict_7(ctx, src, srcSize);
    }
}
