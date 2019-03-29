/*
 * Copyright (c) 2016-present, Yann Collet, Facebook, Inc.
 * All rights reserved.
 *
 * This source code is licensed under both the BSD-style license (found in the
 * LICENSE file in the root directory of this source tree) and the GPLv2 (found
 * in the COPYING file in the root directory of this source tree).
 * You may select, at your option, one of the above-listed licenses.
 */

#include "zstd_compress_internal.h"
#include "zstd_fast.h"


void ZSTD_fillHashTable(ZSTD_matchState_t* ms,
                        void const* end, ZSTD_dictTableLoadMethod_e dtlm)
{
    const ZSTD_compressionParameters* const cParams = &ms->cParams;
    U32* const hashTable = ms->hashTable;
    U32  const hBits = cParams->hashLog;
    U32  const mls = cParams->minMatch;
    const BYTE* const base = ms->window.base;
    const BYTE* ip = base + ms->nextToUpdate;
    const BYTE* const iend = ((const BYTE*)end) - HASH_READ_SIZE;
    const U32 fastHashFillStep = 3;

    /* Always insert every fastHashFillStep position into the hash table.
     * Insert the other positions if their hash entry is empty.
     */
    for ( ; ip + fastHashFillStep < iend + 2; ip += fastHashFillStep) {
        U32 const current = (U32)(ip - base);
        size_t const hash0 = ZSTD_hashPtr(ip, hBits, mls);
        hashTable[hash0] = current;
        if (dtlm == ZSTD_dtlm_fast) continue;
        /* Only load extra positions for ZSTD_dtlm_full */
        {   U32 p;
            for (p = 1; p < fastHashFillStep; ++p) {
                size_t const hash = ZSTD_hashPtr(ip + p, hBits, mls);
                if (hashTable[hash] == 0) {  /* not yet filled */
                    hashTable[hash] = current + p;
    }   }   }   }
}

FORCE_INLINE_TEMPLATE
size_t ZSTD_compressBlock_fast_generic(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize,
        U32 const mls, ZSTD_dictMode_e const dictMode)
{
    const ZSTD_compressionParameters* const cParams = &ms->cParams;
    U32* const hashTable = ms->hashTable;
    U32 const hlog = cParams->hashLog;
    U32 const stepSize = MAX(cParams->targetLength, 2);
    const BYTE* const base = ms->window.base;
    const BYTE* const istart = (const BYTE*)src;
    const BYTE* ip0 = istart;
    const BYTE* ip1;
    const BYTE* anchor = istart;
    const U32   prefixStartIndex = ms->window.dictLimit;
    const BYTE* const prefixStart = base + prefixStartIndex;
    const BYTE* const iend = istart + srcSize;
    const BYTE* const ilimit = iend - HASH_READ_SIZE;
    U32 offset_1=rep[0], offset_2=rep[1];
    U32 offsetSaved = 0;

    const ZSTD_matchState_t* const dms = ms->dictMatchState;
    const ZSTD_compressionParameters* const dictCParams =
                                     dictMode == ZSTD_dictMatchState ?
                                     &dms->cParams : NULL;
    const U32* const dictHashTable = dictMode == ZSTD_dictMatchState ?
                                     dms->hashTable : NULL;
    const U32 dictStartIndex       = dictMode == ZSTD_dictMatchState ?
                                     dms->window.dictLimit : 0;
    const BYTE* const dictBase     = dictMode == ZSTD_dictMatchState ?
                                     dms->window.base : NULL;
    const BYTE* const dictStart    = dictMode == ZSTD_dictMatchState ?
                                     dictBase + dictStartIndex : NULL;
    const BYTE* const dictEnd      = dictMode == ZSTD_dictMatchState ?
                                     dms->window.nextSrc : NULL;
    const U32 dictIndexDelta       = dictMode == ZSTD_dictMatchState ?
                                     prefixStartIndex - (U32)(dictEnd - dictBase) :
                                     0;
    const U32 dictAndPrefixLength  = (U32)(ip0 - prefixStart + dictEnd - dictStart);
    const U32 dictHLog             = dictMode == ZSTD_dictMatchState ?
                                     dictCParams->hashLog : hlog;

    assert(dictMode == ZSTD_noDict || dictMode == ZSTD_dictMatchState);

    /* otherwise, we would get index underflow when translating a dict index
     * into a local index */
    assert(dictMode != ZSTD_dictMatchState
        || prefixStartIndex >= (U32)(dictEnd - dictBase));

    /* init */
    ip0 += (dictAndPrefixLength == 0);
    ip1 = ip0 + 1;
    if (dictMode == ZSTD_noDict) {
        U32 const maxRep = (U32)(ip0 - prefixStart);
        if (offset_2 > maxRep) offsetSaved = offset_2, offset_2 = 0;
        if (offset_1 > maxRep) offsetSaved = offset_1, offset_1 = 0;
    }
    if (dictMode == ZSTD_dictMatchState) {
        /* dictMatchState repCode checks don't currently handle repCode == 0
         * disabling. */
        assert(offset_1 <= dictAndPrefixLength);
        assert(offset_2 <= dictAndPrefixLength);
    }

    /* Main Search Loop */
    while (ip0 < ilimit) {   /* < instead of <=, because repcode check at (ip+1) */
        size_t mLength;
        size_t const h0 = ZSTD_hashPtr(ip0, hlog, mls);
        size_t const h1 = ZSTD_hashPtr(ip1, hlog, mls);
        U32 const val0 = MEM_read32(ip0);
        U32 const val1 = MEM_read32(ip1);
        U32 const current0 = (U32)(ip0-base);
        U32 const current1 = (U32)(ip1-base);
        U32 const matchIndex0 = hashTable[h0];
        U32 const matchIndex1 = hashTable[h1];
        const BYTE* match0 = base + matchIndex0;
        const BYTE* match1 = base + matchIndex1;
        const U32 repIndex = current1 - offset_1;
        const BYTE* repMatch = (dictMode == ZSTD_dictMatchState
                            && repIndex < prefixStartIndex) ?
                               dictBase + (repIndex - dictIndexDelta) :
                               base + repIndex;
        hashTable[h0] = current0;   /* update hash table */
        hashTable[h1] = current1;   /* update hash table */

        if ( (dictMode == ZSTD_dictMatchState)
          && ((U32)((prefixStartIndex-1) - repIndex) >= 3) /* intentional underflow : ensure repIndex isn't overlapping dict + prefix */
          && (MEM_read32(repMatch) == val1) ) {
            const BYTE* const repMatchEnd = repIndex < prefixStartIndex ? dictEnd : iend;
            mLength = ZSTD_count_2segments(ip1+4, repMatch+4, iend, repMatchEnd, prefixStart) + 4;
            ZSTD_storeSeq(seqStore, ip1-anchor, anchor, 0, mLength-MINMATCH);
            ip0++;
            goto _match;
        }
        if ( dictMode == ZSTD_noDict
                 && ((offset_1 > 0) & (MEM_read32(ip1-offset_1) == val1))) {
            mLength = ZSTD_count(ip1+4, ip1+4-offset_1, iend) + 4;
            ZSTD_storeSeq(seqStore, ip1-anchor, anchor, 0, mLength-MINMATCH);
            ip0++;
            goto _match;
        }
        if ( (matchIndex0 <= prefixStartIndex) ) {
            if (dictMode == ZSTD_dictMatchState) {
                size_t const dictHash = ZSTD_hashPtr(ip0, dictHLog, mls);
                U32 const dictMatchIndex = dictHashTable[dictHash];
                const BYTE* dictMatch = dictBase + dictMatchIndex;
                if (dictMatchIndex > dictStartIndex && MEM_read32(dictMatch) == val0) {
                    /* found a dict match */
                    U32 const offset = (U32)(current0-dictMatchIndex-dictIndexDelta);
                    mLength = ZSTD_count_2segments(ip0+4, dictMatch+4, iend, dictEnd, prefixStart) + 4;
                    while (((ip0>anchor) & (dictMatch>dictStart))
                         && (ip0[-1] == dictMatch[-1])) {
                        ip0--; dictMatch--; mLength++;
                    } /* catch up */
                    offset_2 = offset_1;
                    offset_1 = offset;
                    ZSTD_storeSeq(seqStore, ip0-anchor, anchor, offset + ZSTD_REP_MOVE, mLength-MINMATCH);
                    goto _match;
                }
            }
        } else if (MEM_read32(match0) == val0) {
            /* found a regular match */
            U32 const offset = (U32)(ip0-match0);
            mLength = ZSTD_count(ip0+4, match0+4, iend) + 4;
            while (((ip0>anchor) & (match0>prefixStart))
                 && (ip0[-1] == match0[-1])) { ip0--; match0--; mLength++; } /* catch up */
            offset_2 = offset_1;
            offset_1 = offset;
            ZSTD_storeSeq(seqStore, ip0-anchor, anchor, offset + ZSTD_REP_MOVE, mLength-MINMATCH);
            goto _match;
        }
        if ( (matchIndex1 <= prefixStartIndex) ) {
            if (dictMode == ZSTD_dictMatchState) {
                size_t const dictHash = ZSTD_hashPtr(ip1, dictHLog, mls);
                U32 const dictMatchIndex = dictHashTable[dictHash];
                const BYTE* dictMatch = dictBase + dictMatchIndex;
                if (dictMatchIndex > dictStartIndex && MEM_read32(dictMatch) == val1) {
                    /* found a dict match */
                    U32 const offset = (U32)(current1-dictMatchIndex-dictIndexDelta);
                    mLength = ZSTD_count_2segments(ip1+4, dictMatch+4, iend, dictEnd, prefixStart) + 4;
                    while (((ip1>anchor) & (dictMatch>dictStart))
                         && (ip1[-1] == dictMatch[-1])) {
                        ip1--; dictMatch--; mLength++;
                    } /* catch up */
                    offset_2 = offset_1;
                    offset_1 = offset;
                    ZSTD_storeSeq(seqStore, ip1-anchor, anchor, offset + ZSTD_REP_MOVE, mLength-MINMATCH);
                    ip0 = ip1;
                    goto _match;
                }
            }
        } else if (MEM_read32(match1) == val1) {
            /* found a regular match */
            U32 const offset = (U32)(ip1-match1);
            mLength = ZSTD_count(ip1+4, match1+4, iend) + 4;
            while (((ip1>anchor) & (match1>prefixStart))
                 && (ip1[-1] == match1[-1])) { ip1--; match1--; mLength++; } /* catch up */
            offset_2 = offset_1;
            offset_1 = offset;
            ZSTD_storeSeq(seqStore, ip1-anchor, anchor, offset + ZSTD_REP_MOVE, mLength-MINMATCH);
            ip0 = ip1;
            goto _match;
        }
        {
            size_t const step = ((ip0-anchor) >> kSearchStrength) + stepSize;
            assert(stepSize >= 2);
            ip0 += step;
            ip1 += step;
            continue;
        }
_match:
        /* match found */
        ip0 += mLength;
        anchor = ip0;
        ip1 = ip0 + 1;

        if (ip0 <= ilimit) {
            /* Fill Table */
            assert(base+current0+2 > istart);  /* check base overflow */
            hashTable[ZSTD_hashPtr(base+current0+2, hlog, mls)] = current0+2;  /* here because current+2 could be > iend-8 */
            hashTable[ZSTD_hashPtr(ip0-2, hlog, mls)] = (U32)(ip0-2-base);

            /* check immediate repcode */
            if (dictMode == ZSTD_dictMatchState) {
                while (ip0 <= ilimit) {
                    U32 const current2 = (U32)(ip0-base);
                    U32 const repIndex2 = current2 - offset_2;
                    const BYTE* repMatch2 = repIndex2 < prefixStartIndex ?
                            dictBase - dictIndexDelta + repIndex2 :
                            base + repIndex2;
                    if ( ((U32)((prefixStartIndex-1) - (U32)repIndex2) >= 3 /* intentional overflow */)
                       && (MEM_read32(repMatch2) == MEM_read32(ip0)) ) {
                        const BYTE* const repEnd2 = repIndex2 < prefixStartIndex ? dictEnd : iend;
                        size_t const repLength2 = ZSTD_count_2segments(ip0+4, repMatch2+4, iend, repEnd2, prefixStart) + 4;
                        U32 tmpOffset = offset_2; offset_2 = offset_1; offset_1 = tmpOffset;   /* swap offset_2 <=> offset_1 */
                        ZSTD_storeSeq(seqStore, 0, anchor, 0, repLength2-MINMATCH);
                        hashTable[ZSTD_hashPtr(ip0, hlog, mls)] = current2;
                        ip0 += repLength2;
                        anchor = ip0;
                        continue;
                    }
                    break;
                }
            }

            if (dictMode == ZSTD_noDict) {
                while ( (ip0 <= ilimit)
                     && ( (offset_2>0)
                        & (MEM_read32(ip0) == MEM_read32(ip0 - offset_2)) )) {
                    /* store sequence */
                    size_t const rLength = ZSTD_count(ip0+4, ip0+4-offset_2, iend) + 4;
                    U32 const tmpOff = offset_2; offset_2 = offset_1; offset_1 = tmpOff;  /* swap offset_2 <=> offset_1 */
                    hashTable[ZSTD_hashPtr(ip0, hlog, mls)] = (U32)(ip0-base);
                    ZSTD_storeSeq(seqStore, 0, anchor, 0, rLength-MINMATCH);
                    ip0 += rLength;
                    anchor = ip0;
                }
            }
            ip1 = ip0 + 1;
        }
    }

    /* save reps for next block */
    rep[0] = offset_1 ? offset_1 : offsetSaved;
    rep[1] = offset_2 ? offset_2 : offsetSaved;

    /* Return the last literals size */
    return iend - anchor;
}


size_t ZSTD_compressBlock_fast(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize)
{
    ZSTD_compressionParameters const* cParams = &ms->cParams;
    U32 const mls = cParams->minMatch;
    assert(ms->dictMatchState == NULL);
    switch(mls)
    {
    default: /* includes case 3 */
    case 4 :
        return ZSTD_compressBlock_fast_generic(ms, seqStore, rep, src, srcSize, 4, ZSTD_noDict);
    case 5 :
        return ZSTD_compressBlock_fast_generic(ms, seqStore, rep, src, srcSize, 5, ZSTD_noDict);
    case 6 :
        return ZSTD_compressBlock_fast_generic(ms, seqStore, rep, src, srcSize, 6, ZSTD_noDict);
    case 7 :
        return ZSTD_compressBlock_fast_generic(ms, seqStore, rep, src, srcSize, 7, ZSTD_noDict);
    }
}

size_t ZSTD_compressBlock_fast_dictMatchState(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize)
{
    ZSTD_compressionParameters const* cParams = &ms->cParams;
    U32 const mls = cParams->minMatch;
    assert(ms->dictMatchState != NULL);
    switch(mls)
    {
    default: /* includes case 3 */
    case 4 :
        return ZSTD_compressBlock_fast_generic(ms, seqStore, rep, src, srcSize, 4, ZSTD_dictMatchState);
    case 5 :
        return ZSTD_compressBlock_fast_generic(ms, seqStore, rep, src, srcSize, 5, ZSTD_dictMatchState);
    case 6 :
        return ZSTD_compressBlock_fast_generic(ms, seqStore, rep, src, srcSize, 6, ZSTD_dictMatchState);
    case 7 :
        return ZSTD_compressBlock_fast_generic(ms, seqStore, rep, src, srcSize, 7, ZSTD_dictMatchState);
    }
}


static size_t ZSTD_compressBlock_fast_extDict_generic(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize, U32 const mls)
{
    const ZSTD_compressionParameters* const cParams = &ms->cParams;
    U32* const hashTable = ms->hashTable;
    U32 const hlog = cParams->hashLog;
    /* support stepSize of 0 */
    U32 const stepSize = cParams->targetLength + !(cParams->targetLength);
    const BYTE* const base = ms->window.base;
    const BYTE* const dictBase = ms->window.dictBase;
    const BYTE* const istart = (const BYTE*)src;
    const BYTE* ip = istart;
    const BYTE* anchor = istart;
    const U32   dictStartIndex = ms->window.lowLimit;
    const BYTE* const dictStart = dictBase + dictStartIndex;
    const U32   prefixStartIndex = ms->window.dictLimit;
    const BYTE* const prefixStart = base + prefixStartIndex;
    const BYTE* const dictEnd = dictBase + prefixStartIndex;
    const BYTE* const iend = istart + srcSize;
    const BYTE* const ilimit = iend - 8;
    U32 offset_1=rep[0], offset_2=rep[1];

    /* Search Loop */
    while (ip < ilimit) {  /* < instead of <=, because (ip+1) */
        const size_t h = ZSTD_hashPtr(ip, hlog, mls);
        const U32    matchIndex = hashTable[h];
        const BYTE* const matchBase = matchIndex < prefixStartIndex ? dictBase : base;
        const BYTE*  match = matchBase + matchIndex;
        const U32    current = (U32)(ip-base);
        const U32    repIndex = current + 1 - offset_1;
        const BYTE* const repBase = repIndex < prefixStartIndex ? dictBase : base;
        const BYTE* const repMatch = repBase + repIndex;
        size_t mLength;
        hashTable[h] = current;   /* update hash table */
        assert(offset_1 <= current +1);   /* check repIndex */

        if ( (((U32)((prefixStartIndex-1) - repIndex) >= 3) /* intentional underflow */ & (repIndex > dictStartIndex))
           && (MEM_read32(repMatch) == MEM_read32(ip+1)) ) {
            const BYTE* repMatchEnd = repIndex < prefixStartIndex ? dictEnd : iend;
            mLength = ZSTD_count_2segments(ip+1+4, repMatch+4, iend, repMatchEnd, prefixStart) + 4;
            ip++;
            ZSTD_storeSeq(seqStore, ip-anchor, anchor, 0, mLength-MINMATCH);
        } else {
            if ( (matchIndex < dictStartIndex) ||
                 (MEM_read32(match) != MEM_read32(ip)) ) {
                assert(stepSize >= 1);
                ip += ((ip-anchor) >> kSearchStrength) + stepSize;
                continue;
            }
            {   const BYTE* matchEnd = matchIndex < prefixStartIndex ? dictEnd : iend;
                const BYTE* lowMatchPtr = matchIndex < prefixStartIndex ? dictStart : prefixStart;
                U32 offset;
                mLength = ZSTD_count_2segments(ip+4, match+4, iend, matchEnd, prefixStart) + 4;
                while (((ip>anchor) & (match>lowMatchPtr)) && (ip[-1] == match[-1])) { ip--; match--; mLength++; }   /* catch up */
                offset = current - matchIndex;
                offset_2 = offset_1;
                offset_1 = offset;
                ZSTD_storeSeq(seqStore, ip-anchor, anchor, offset + ZSTD_REP_MOVE, mLength-MINMATCH);
        }   }

        /* found a match : store it */
        ip += mLength;
        anchor = ip;

        if (ip <= ilimit) {
            /* Fill Table */
            hashTable[ZSTD_hashPtr(base+current+2, hlog, mls)] = current+2;
            hashTable[ZSTD_hashPtr(ip-2, hlog, mls)] = (U32)(ip-2-base);
            /* check immediate repcode */
            while (ip <= ilimit) {
                U32 const current2 = (U32)(ip-base);
                U32 const repIndex2 = current2 - offset_2;
                const BYTE* repMatch2 = repIndex2 < prefixStartIndex ? dictBase + repIndex2 : base + repIndex2;
                if ( (((U32)((prefixStartIndex-1) - repIndex2) >= 3) & (repIndex2 > dictStartIndex))  /* intentional overflow */
                   && (MEM_read32(repMatch2) == MEM_read32(ip)) ) {
                    const BYTE* const repEnd2 = repIndex2 < prefixStartIndex ? dictEnd : iend;
                    size_t const repLength2 = ZSTD_count_2segments(ip+4, repMatch2+4, iend, repEnd2, prefixStart) + 4;
                    U32 tmpOffset = offset_2; offset_2 = offset_1; offset_1 = tmpOffset;   /* swap offset_2 <=> offset_1 */
                    ZSTD_storeSeq(seqStore, 0, anchor, 0, repLength2-MINMATCH);
                    hashTable[ZSTD_hashPtr(ip, hlog, mls)] = current2;
                    ip += repLength2;
                    anchor = ip;
                    continue;
                }
                break;
    }   }   }

    /* save reps for next block */
    rep[0] = offset_1;
    rep[1] = offset_2;

    /* Return the last literals size */
    return iend - anchor;
}


size_t ZSTD_compressBlock_fast_extDict(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize)
{
    ZSTD_compressionParameters const* cParams = &ms->cParams;
    U32 const mls = cParams->minMatch;
    switch(mls)
    {
    default: /* includes case 3 */
    case 4 :
        return ZSTD_compressBlock_fast_extDict_generic(ms, seqStore, rep, src, srcSize, 4);
    case 5 :
        return ZSTD_compressBlock_fast_extDict_generic(ms, seqStore, rep, src, srcSize, 5);
    case 6 :
        return ZSTD_compressBlock_fast_extDict_generic(ms, seqStore, rep, src, srcSize, 6);
    case 7 :
        return ZSTD_compressBlock_fast_extDict_generic(ms, seqStore, rep, src, srcSize, 7);
    }
}
