/*
 * Copyright (c) Yann Collet, Facebook, Inc.
 * All rights reserved.
 *
 * This source code is licensed under both the BSD-style license (found in the
 * LICENSE file in the root directory of this source tree) and the GPLv2 (found
 * in the COPYING file in the root directory of this source tree).
 * You may select, at your option, one of the above-listed licenses.
 */

#include "zstd_compress_internal.h"  /* ZSTD_hashPtr, ZSTD_count, ZSTD_storeSeq */
#include "zstd_fast.h"

typedef U32 ZSTD_VecMask;

/* ZSTD_VecMask_next():
 * Starting from the LSB, returns the idx of the next non-zero bit.
 * Basically counting the nb of trailing zeroes.
 */
MEM_STATIC U32 ZSTD_VecMask_next(ZSTD_VecMask val) {
    return ZSTD_countTrailingZeros32(val);
}

FORCE_INLINE_TEMPLATE U32 ZSTD_getTagMask(void const* haystack, BYTE const needle)
{
    __m128i const validMask  = _mm_set_epi16(0, 0, 0, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF);
    __m128i const needleMask = _mm_set1_epi8((char)needle);
    __m128i const haystackMask = _mm_load_si128((const __m128i*)haystack);
    __m128i const equalMask = _mm_cmpeq_epi8(needleMask, haystackMask);
    __m128i const validEqualMask = _mm_and_si128(equalMask, validMask);
    return _mm_movemask_epi8(validEqualMask);
}

#ifdef ZSTD_USE_DICT
#include <stdio.h>

static uint32_t sortedDictPos[128 * 1024] = {};
static uint32_t rSortedDictPos[128 * 1024] = {};
static size_t loadDictPos(void)
{
    FILE* f = fopen("dict-uses.u32", "rb");
    assert(f);
    {
        size_t const nbPos = fread(sortedDictPos, sizeof(sortedDictPos[0]), 128 * 1024, f);
        assert(nbPos > 0);
        fclose(f);
        {
            size_t i;
            for (i = 0; i < nbPos; ++i) {
                rSortedDictPos[sortedDictPos[i]] = i;
            }
        }
        return nbPos;
    }
}

void ZSTD_fillHashTable(ZSTD_matchState_t* ms,
                        const void* const end,
                        ZSTD_dictTableLoadMethod_e dtlm)
{
    const ZSTD_compressionParameters* const cParams = &ms->cParams;
    U32* const hashTable = ms->hashTable;
    U32  const hBits = cParams->hashLog;
    U32  const mls = cParams->minMatch;
    const BYTE* const base = ms->window.base;
    const BYTE* const dictStart = base + ms->window.dictLimit;
    const BYTE* const dictEnd = ms->window.nextSrc;
    const BYTE* const iend = ((const BYTE*)end) - HASH_READ_SIZE;

    {
        const BYTE* ip = base + ms->nextToUpdate;
        const U32 fastHashFillStep = 3;
        for ( ; ip + fastHashFillStep < iend + 2; ip += fastHashFillStep) {
            U32 const curr = (U32)(ip - base);
            size_t const hash0 = ZSTD_hashPtr(ip, hBits, mls);
            hashTable[hash0] = curr;
            if (dtlm == ZSTD_dtlm_fast) continue;
            /* Only load extra positions for ZSTD_dtlm_full */
            {   U32 p;
                for (p = 1; p < fastHashFillStep; ++p) {
                    size_t const hash = ZSTD_hashPtr(ip + p, hBits, mls);
                    if (hashTable[hash] == 0) {  /* not yet filled */
                        hashTable[hash] = curr + p;
        }   }   }   }
    }

    size_t const nbPos = loadDictPos();
    size_t i;

    assert(ms->window.dictLimit == ms->nextToUpdate);

    for (i = 0; i < nbPos; ++i) {
        U32 pos = sortedDictPos[i];
        BYTE const* const ip = dictStart + pos;
        U32 const curr = (U32)(ip - base);
        assert(ip < iend);
        assert(ip < dictEnd);
        {
            size_t const hash = ZSTD_hashPtr(ip, hBits, mls);
            hashTable[hash] = curr;
        }
    }
    (void)dtlm;
}

static int ZSTD_insertToRow(U16* row, BYTE tag, U16 index, int replace)
{
    BYTE* tags = (BYTE*)row;
    U16* idxs = row + 5;
    size_t i;
    U16 minVal = (U16)-1;
    U16 minIdx = 0;
    assert(index != 0);
    if (replace) {
        ZSTD_VecMask const mask = ZSTD_getTagMask(tags, tag);
        if (mask != 0) {
            U32 const idx = ZSTD_VecMask_next(mask);
            assert(tags[idx] == tag);
            idxs[idx] = index;
            return 3;
        }
    }
    for (i = 0; i < 10; ++i) {
        if (idxs[i] == 0) {
            assert(tags[i] == 0);
            idxs[i] = index;
            tags[i] = tag;
            return 2;
        }
        if (rSortedDictPos[idxs[i]] < minVal) {
            minVal = rSortedDictPos[idxs[i]];
            minIdx = i;
        }
    }
    if (replace) {
        tags[minIdx] = tag;
        idxs[minIdx] = index;
        return 1;
    }
    return 0;
}
void ZSTD_fillHashTableDDS(ZSTD_matchState_t* ms,
                        const void* const end)
{
    const ZSTD_compressionParameters* const cParams = &ms->cParams;
    U16* const hashTable = (U16*)(void*)ms->hashTable;
    U32  const rowLog    = (cParams->hashLog + 2) - 5;
    U32  const tagLog    = 8;
    U32  const hashLog   = rowLog + tagLog;
    U32  const mls       = cParams->minMatch;
    const BYTE* const base = ms->window.base;
    const BYTE* const istart = base + ms->nextToUpdate;
    const BYTE* const dictStart = base + ms->window.dictLimit;
    const BYTE* const dictEnd = ms->window.nextSrc;
    const BYTE* ip = istart;
    const BYTE* const iend = ((const BYTE*)end) - HASH_READ_SIZE;
    const U32 fastHashFillStep = 3;

    assert(iend - base < (1 << 16) - 1);

    // for ( ; ip < iend - (fastHashFillStep - 1); ip += fastHashFillStep) {
    //     U16    const curr = (U16)(ip - base);
    //     size_t const hash = ZSTD_hashPtr(ip, hashLog, mls);
    //     size_t const row = (hash >> tagLog) << 4;
    //     size_t const tag = hash & 0xFF;
    //     ZSTD_insertToRow(hashTable + row, tag, curr, /* overwrite */ 1);
    // }
    {
        size_t const nbPos = loadDictPos();
        size_t i;
        size_t nKicked = 0;

        assert(ms->window.dictLimit == ms->nextToUpdate);

        for (i = 0; i < nbPos; ++i) {
            U32 pos = sortedDictPos[i];
            ip = dictStart + pos;
            U32 const curr = (U32)(ip - base);
            assert(ip < iend);
            assert(ip < dictEnd);
            {
                size_t const hash = ZSTD_hashPtr(ip, hashLog, mls);
                size_t const row = (hash >> tagLog) << 4;
                size_t const tag = hash & 0xFF;
                nKicked += 1 == ZSTD_insertToRow(hashTable + row, tag, curr, /* overwrite */ 1);
            }
        }
        DEBUGLOG(2, "nkicked = %zu", nKicked);
    }

    size_t nExtra = 0;
    size_t nEmpty = 0;
    for (ip = istart; ip < iend; ++ip) {
        U16    const curr = (U16)(ip - base);
        size_t const hash = ZSTD_hashPtr(ip, hashLog, mls);
        size_t const row = (hash >> tagLog) << 4;
        size_t const tag = hash & 0xFF;
        nExtra += 0 != ZSTD_insertToRow(hashTable + row, tag, curr, /* overwrite */ 0);
    }
    {
        U32 r;
        U32 const e = 1u << rowLog;
        U16 const dictStartIndex = ms->window.dictLimit;
        for (r = 0; r < e; ++r) {
            size_t const row = r << 4;
            while (ZSTD_insertToRow(hashTable + row, 0, dictStartIndex, /* overwrite */ 0)) {
                ++nEmpty;
            }
        }
    }
    DEBUGLOG(2, "nextra = %zu", nExtra);
    DEBUGLOG(2, "nempty = %zu", nEmpty);
}
#else
void ZSTD_fillHashTable(ZSTD_matchState_t* ms,
                        const void* const end,
                        ZSTD_dictTableLoadMethod_e dtlm)
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
        U32 const curr = (U32)(ip - base);
        size_t const hash0 = ZSTD_hashPtr(ip, hBits, mls);
        hashTable[hash0] = curr;
        if (dtlm == ZSTD_dtlm_fast) continue;
        /* Only load extra positions for ZSTD_dtlm_full */
        {   U32 p;
            for (p = 1; p < fastHashFillStep; ++p) {
                size_t const hash = ZSTD_hashPtr(ip + p, hBits, mls);
                if (hashTable[hash] == 0) {  /* not yet filled */
                    hashTable[hash] = curr + p;
    }   }   }   }
}


static int ZSTD_insertToRow(U16* row, BYTE tag, U16 index, int replace)
{
    BYTE* tags = (BYTE*)row;
    U16* idxs = row + 5;
    size_t i;
    U16 minVal = (U16)-1;
    U16 minIdx = 0;
    assert(index != 0);
    if (replace) {
        ZSTD_VecMask const mask = ZSTD_getTagMask(tags, tag);
        if (mask != 0) {
            U32 const idx = ZSTD_VecMask_next(mask);
            assert(idx < 10);
            assert(tags[idx] == tag);
            idxs[idx] = index;
            return 1;
        }
    }
    for (i = 0; i < 10; ++i) {
        if (idxs[i] == 0) {
            assert(tags[i] == 0);
            idxs[i] = index;
            tags[i] = tag;
            return 1;
        }
        if (idxs[i] < minVal) {
            minVal = idxs[i];
            minIdx = i;
        }
    }
    if (replace) {
        tags[minIdx] = tag;
        idxs[minIdx] = index;
        return 1;
    }
    return 0;
}
void ZSTD_fillHashTableDDS(ZSTD_matchState_t* ms,
                        const void* const end)
{
    const ZSTD_compressionParameters* const cParams = &ms->cParams;
    U16* const hashTable = (U16*)(void*)ms->hashTable;
    U32  const rowLog    = (cParams->hashLog + 2) - 5;
    U32  const tagLog    = 8;
    U32  const hashLog   = rowLog + tagLog;
    U32  const mls       = cParams->minMatch;
    const BYTE* const base = ms->window.base;
    const BYTE* const istart = base + ms->nextToUpdate;
    const BYTE* ip = istart;
    const BYTE* const iend = ((const BYTE*)end) - HASH_READ_SIZE;
    const U32 fastHashFillStep = 3;

    assert(iend - base < (1 << 16) - 1);

    for ( ; ip < iend - (fastHashFillStep - 1); ip += fastHashFillStep) {
        U32 u;
        {
            U16    const curr = (U16)(ip - base);
            size_t const hash = ZSTD_hashPtr(ip, hashLog, mls);
            size_t const row = (hash >> tagLog) << 4;
            size_t const tag = hash & 0xFF;
            ZSTD_insertToRow(hashTable + row, tag, curr, /* overwrite */ 1);
        }
        for (u = 1; u < fastHashFillStep; ++u) {
            U16    const curr = (U16)(ip + u - base);
            size_t const hash = ZSTD_hashPtr(ip + u, hashLog, mls);
            size_t const row = (hash >> tagLog) << 4;
            size_t const tag = hash & 0xFF;
            ZSTD_insertToRow(hashTable + row, tag, curr, /* overwrite */ 0);
        }
    }
    // for ( ip = iend - fastHashFillStep ; ip >= istart; ip -= fastHashFillStep) {
    //     U16    const curr = (U16)(ip - base);
    //     size_t const hash = ZSTD_hashPtr(ip, hashLog, mls);
    //     size_t const row = (hash >> tagLog) << 4;
    //     size_t const tag = hash & 0xFF;
    //     ZSTD_insertToRow(hashTable + row, tag, curr, /* overwrite */ 0);
    // }
    // ip += 1;
    // if (ip < istart) ip += fastHashFillStep;
    // for (; ip < iend - fastHashFillStep; ip += fastHashFillStep) {
    //     U16    const curr = (U16)(ip - base);
    //     size_t const hash = ZSTD_hashPtr(ip, hashLog, mls);
    //     size_t const row = (hash >> tagLog) << 4;
    //     size_t const tag = hash & 0xFF;
    //     ZSTD_insertToRow(hashTable + row, tag, curr, /* overwrite */ 0);
    // }
    // assert(ip == iend - fastHashFillStep + 1);
    // ip += 1;
    // for ( ; ip >= istart; ip -= fastHashFillStep) {
    //     U16    const curr = (U16)(ip - base);
    //     size_t const hash = ZSTD_hashPtr(ip, hashLog, mls);
    //     size_t const row = (hash >> tagLog) << 4;
    //     size_t const tag = hash & 0xFF;
    //     ZSTD_insertToRow(hashTable + row, tag, curr, /* overwrite */ 0);
    // }
    {
        U32 r;
        U32 const e = 1u << rowLog;
        U16 const dictStartIndex = ms->window.dictLimit;
        for (r = 0; r < e; ++r) {
            size_t const row = r << 4;
            while (ZSTD_insertToRow(hashTable + row, 0, dictStartIndex, /* overwrite */ 0)) {
            }
        }
    }
}

// Fill as close to the previous version as possible
// Basically construct the old hash table, and insert
// if it was in the old hash table.
void ZSTD_fillHashTableDDS2(ZSTD_matchState_t* ms,
                        const void* const end)
{
    const ZSTD_compressionParameters* const cParams = &ms->cParams;
    U16* const hashTable = (U16*)(void*)ms->hashTable;
    U32  const rowLog    = (cParams->hashLog + 2) - 5;
    U32  const tagLog    = 8;
    U32  const hashLog   = rowLog + tagLog;
    U32  const mls       = cParams->minMatch;
    const BYTE* const base = ms->window.base;
    const BYTE* const istart = base + ms->nextToUpdate;
    const BYTE* ip = istart;
    const BYTE* const iend = ((const BYTE*)end) - HASH_READ_SIZE;
    const U32 fastHashFillStep = 3;

    U32* tmpTable = malloc(sizeof(uint32_t) << hashLog);
    ZSTD_memset(tmpTable, 0, sizeof(uint32_t) << hashLog);

    /* Always insert every fastHashFillStep position into the hash table.
     * Insert the other positions if their hash entry is empty.
     */
    for ( ip = istart ; ip + fastHashFillStep < iend + 2; ip += fastHashFillStep) {
        U32 const curr = (U32)(ip - base);
        size_t const hash0 = ZSTD_hashPtr(ip, cParams->hashLog, mls);
        tmpTable[hash0] = curr;
        /* Only load extra positions for ZSTD_dtlm_full */
        {   U32 p;
            for (p = 1; p < fastHashFillStep; ++p) {
                size_t const hash = ZSTD_hashPtr(ip + p, cParams->hashLog, mls);
                if (tmpTable[hash] == 0) {  /* not yet filled */
                    tmpTable[hash] = curr + p;
    }   }   }   }


    assert(iend - base < (1 << 16) - 1);
    size_t missed = 0;
    size_t total = 0;
    for (ip = iend - 1; ip >= istart; --ip) {
        U32 const curr = (U32)(ip - base);
        {
            size_t const hash = ZSTD_hashPtr(ip, cParams->hashLog, mls);
            if (tmpTable[hash] != curr)
                continue;
        }
        size_t const hash = ZSTD_hashPtr(ip, hashLog, mls);
        size_t const row = (hash >> tagLog) << 4;
        size_t const tag = hash & 0xFF;
        missed += 0 == ZSTD_insertToRow(hashTable + row, tag, curr, /* overwrite */ 0);
        ++total;
    }
    DEBUGLOG(2, "missed = %zu / total = %zu", missed, total);

    // for ( ; ip < iend - (fastHashFillStep - 1); ip += fastHashFillStep) {
    // for ( ip = iend - fastHashFillStep ; ip >= istart; ip -= fastHashFillStep) {
    //     U16    const curr = (U16)(ip - base);
    //     size_t const hash = ZSTD_hashPtr(ip, hashLog, mls);
    //     size_t const row = (hash >> tagLog) << 4;
    //     size_t const tag = hash & 0xFF;
    //     ZSTD_insertToRow(hashTable + row, tag, curr, /* overwrite */ 0);
    // }
    // ip += 1;
    // if (ip < istart) ip += fastHashFillStep;
    // for (; ip < iend - fastHashFillStep; ip += fastHashFillStep) {
    //     U16    const curr = (U16)(ip - base);
    //     size_t const hash = ZSTD_hashPtr(ip, hashLog, mls);
    //     size_t const row = (hash >> tagLog) << 4;
    //     size_t const tag = hash & 0xFF;
    //     ZSTD_insertToRow(hashTable + row, tag, curr, /* overwrite */ 0);
    // }
    // assert(ip == iend - fastHashFillStep + 1);
    // ip += 1;
    // for ( ; ip >= istart; ip -= fastHashFillStep) {
    //     U16    const curr = (U16)(ip - base);
    //     size_t const hash = ZSTD_hashPtr(ip, hashLog, mls);
    //     size_t const row = (hash >> tagLog) << 4;
    //     size_t const tag = hash & 0xFF;
    //     ZSTD_insertToRow(hashTable + row, tag, curr, /* overwrite */ 0);
    // }
    {
        U32 r;
        U32 const e = 1u << rowLog;
        U16 const dictStartIndex = ms->window.dictLimit;
        size_t empty = 0;
        for (r = 0; r < e; ++r) {
            size_t const row = r << 4;
            while (ZSTD_insertToRow(hashTable + row, 0, dictStartIndex, /* overwrite */ 0)) {
                ++empty;
            }
        }
        DEBUGLOG(2, "empty = %zu", empty);
        DEBUGLOG(2, "entries = %u", e * 10);
    }
    free(tmpTable);
}
#endif


/**
 * If you squint hard enough (and ignore repcodes), the search operation at any
 * given position is broken into 4 stages:
 *
 * 1. Hash   (map position to hash value via input read)
 * 2. Lookup (map hash val to index via hashtable read)
 * 3. Load   (map index to value at that position via input read)
 * 4. Compare
 *
 * Each of these steps involves a memory read at an address which is computed
 * from the previous step. This means these steps must be sequenced and their
 * latencies are cumulative.
 *
 * Rather than do 1->2->3->4 sequentially for a single position before moving
 * onto the next, this implementation interleaves these operations across the
 * next few positions:
 *
 * R = Repcode Read & Compare
 * H = Hash
 * T = Table Lookup
 * M = Match Read & Compare
 *
 * Pos | Time -->
 * ----+-------------------
 * N   | ... M
 * N+1 | ...   TM
 * N+2 |    R H   T M
 * N+3 |         H    TM
 * N+4 |           R H   T M
 * N+5 |                H   ...
 * N+6 |                  R ...
 *
 * This is very much analogous to the pipelining of execution in a CPU. And just
 * like a CPU, we have to dump the pipeline when we find a match (i.e., take a
 * branch).
 *
 * When this happens, we throw away our current state, and do the following prep
 * to re-enter the loop:
 *
 * Pos | Time -->
 * ----+-------------------
 * N   | H T
 * N+1 |  H
 *
 * This is also the work we do at the beginning to enter the loop initially.
 */
FORCE_INLINE_TEMPLATE size_t
ZSTD_compressBlock_fast_noDict_generic(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize,
        U32 const mls, U32 const hasStep)
{
    const ZSTD_compressionParameters* const cParams = &ms->cParams;
    U32* const hashTable = ms->hashTable;
    U32 const hlog = cParams->hashLog;
    /* support stepSize of 0 */
    size_t const stepSize = hasStep ? (cParams->targetLength + !(cParams->targetLength) + 1) : 2;
    const BYTE* const base = ms->window.base;
    const BYTE* const istart = (const BYTE*)src;
    const U32   endIndex = (U32)((size_t)(istart - base) + srcSize);
    const U32   prefixStartIndex = ZSTD_getLowestPrefixIndex(ms, endIndex, cParams->windowLog);
    const BYTE* const prefixStart = base + prefixStartIndex;
    const BYTE* const iend = istart + srcSize;
    const BYTE* const ilimit = iend - HASH_READ_SIZE;

    const BYTE* anchor = istart;
    const BYTE* ip0 = istart;
    const BYTE* ip1;
    const BYTE* ip2;
    const BYTE* ip3;
    U32 current0;

    U32 rep_offset1 = rep[0];
    U32 rep_offset2 = rep[1];
    U32 offsetSaved = 0;

    size_t hash0; /* hash for ip0 */
    size_t hash1; /* hash for ip1 */
    U32 idx; /* match idx for ip0 */
    U32 mval; /* src value at match idx */

    U32 offcode;
    const BYTE* match0;
    size_t mLength;

    /* ip0 and ip1 are always adjacent. The targetLength skipping and
     * uncompressibility acceleration is applied to every other position,
     * matching the behavior of #1562. step therefore represents the gap
     * between pairs of positions, from ip0 to ip2 or ip1 to ip3. */
    size_t step;
    const BYTE* nextStep;
    const size_t kStepIncr = (1 << (kSearchStrength - 1));

    DEBUGLOG(5, "ZSTD_compressBlock_fast_generic");
    ip0 += (ip0 == prefixStart);
    {   U32 const curr = (U32)(ip0 - base);
        U32 const windowLow = ZSTD_getLowestPrefixIndex(ms, curr, cParams->windowLog);
        U32 const maxRep = curr - windowLow;
        if (rep_offset2 > maxRep) offsetSaved = rep_offset2, rep_offset2 = 0;
        if (rep_offset1 > maxRep) offsetSaved = rep_offset1, rep_offset1 = 0;
    }

    /* start each op */
_start: /* Requires: ip0 */

    step = stepSize;
    nextStep = ip0 + kStepIncr;

    /* calculate positions, ip0 - anchor == 0, so we skip step calc */
    ip1 = ip0 + 1;
    ip2 = ip0 + step;
    ip3 = ip2 + 1;

    if (ip3 >= ilimit) {
        goto _cleanup;
    }

    hash0 = ZSTD_hashPtr(ip0, hlog, mls);
    hash1 = ZSTD_hashPtr(ip1, hlog, mls);

    idx = hashTable[hash0];

    do {
        /* load repcode match for ip[2]*/
        const U32 rval = MEM_read32(ip2 - rep_offset1);

        /* write back hash table entry */
        current0 = (U32)(ip0 - base);
        hashTable[hash0] = current0;

        /* check repcode at ip[2] */
        if ((MEM_read32(ip2) == rval) & (rep_offset1 > 0)) {
            ip0 = ip2;
            match0 = ip0 - rep_offset1;
            mLength = ip0[-1] == match0[-1];
            ip0 -= mLength;
            match0 -= mLength;
            offcode = REPCODE1_TO_OFFBASE;
            mLength += 4;
            goto _match;
        }

        /* load match for ip[0] */
        if (idx >= prefixStartIndex) {
            mval = MEM_read32(base + idx);
        } else {
            mval = MEM_read32(ip0) ^ 1; /* guaranteed to not match. */
        }

        /* check match at ip[0] */
        if (MEM_read32(ip0) == mval) {
            /* found a match! */
            goto _offset;
        }

        /* lookup ip[1] */
        idx = hashTable[hash1];

        /* hash ip[2] */
        hash0 = hash1;
        hash1 = ZSTD_hashPtr(ip2, hlog, mls);

        /* advance to next positions */
        ip0 = ip1;
        ip1 = ip2;
        ip2 = ip3;

        /* write back hash table entry */
        current0 = (U32)(ip0 - base);
        hashTable[hash0] = current0;

        /* load match for ip[0] */
        if (idx >= prefixStartIndex) {
            mval = MEM_read32(base + idx);
        } else {
            mval = MEM_read32(ip0) ^ 1; /* guaranteed to not match. */
        }

        /* check match at ip[0] */
        if (MEM_read32(ip0) == mval) {
            /* found a match! */
            goto _offset;
        }

        /* lookup ip[1] */
        idx = hashTable[hash1];

        /* hash ip[2] */
        hash0 = hash1;
        hash1 = ZSTD_hashPtr(ip2, hlog, mls);

        /* advance to next positions */
        ip0 = ip1;
        ip1 = ip2;
        ip2 = ip0 + step;
        ip3 = ip1 + step;

        /* calculate step */
        if (ip2 >= nextStep) {
            step++;
            PREFETCH_L1(ip1 + 64);
            PREFETCH_L1(ip1 + 128);
            nextStep += kStepIncr;
        }
    } while (ip3 < ilimit);

_cleanup:
    /* Note that there are probably still a couple positions we could search.
     * However, it seems to be a meaningful performance hit to try to search
     * them. So let's not. */

    /* save reps for next block */
    rep[0] = rep_offset1 ? rep_offset1 : offsetSaved;
    rep[1] = rep_offset2 ? rep_offset2 : offsetSaved;

    /* Return the last literals size */
    return (size_t)(iend - anchor);

_offset: /* Requires: ip0, idx */

    /* Compute the offset code. */
    match0 = base + idx;
    rep_offset2 = rep_offset1;
    rep_offset1 = (U32)(ip0-match0);
    offcode = OFFSET_TO_OFFBASE(rep_offset1);
    mLength = 4;

    /* Count the backwards match length. */
    while (((ip0>anchor) & (match0>prefixStart)) && (ip0[-1] == match0[-1])) {
        ip0--;
        match0--;
        mLength++;
    }

_match: /* Requires: ip0, match0, offcode */

    /* Count the forward length. */
    mLength += ZSTD_count(ip0 + mLength, match0 + mLength, iend);

    ZSTD_storeSeq(seqStore, (size_t)(ip0 - anchor), anchor, iend, offcode, mLength);

    ip0 += mLength;
    anchor = ip0;

    /* write next hash table entry */
    if (ip1 < ip0) {
        hashTable[hash1] = (U32)(ip1 - base);
    }

    /* Fill table and check for immediate repcode. */
    if (ip0 <= ilimit) {
        /* Fill Table */
        assert(base+current0+2 > istart);  /* check base overflow */
        hashTable[ZSTD_hashPtr(base+current0+2, hlog, mls)] = current0+2;  /* here because current+2 could be > iend-8 */
        hashTable[ZSTD_hashPtr(ip0-2, hlog, mls)] = (U32)(ip0-2-base);

        if (rep_offset2 > 0) { /* rep_offset2==0 means rep_offset2 is invalidated */
            while ( (ip0 <= ilimit) && (MEM_read32(ip0) == MEM_read32(ip0 - rep_offset2)) ) {
                /* store sequence */
                size_t const rLength = ZSTD_count(ip0+4, ip0+4-rep_offset2, iend) + 4;
                { U32 const tmpOff = rep_offset2; rep_offset2 = rep_offset1; rep_offset1 = tmpOff; } /* swap rep_offset2 <=> rep_offset1 */
                hashTable[ZSTD_hashPtr(ip0, hlog, mls)] = (U32)(ip0-base);
                ip0 += rLength;
                ZSTD_storeSeq(seqStore, 0 /*litLen*/, anchor, iend, REPCODE1_TO_OFFBASE, rLength);
                anchor = ip0;
                continue;   /* faster when present (confirmed on gcc-8) ... (?) */
    }   }   }

    goto _start;
}

#define ZSTD_GEN_FAST_FN(dictMode, mls, step)                                                            \
    static __attribute__((noinline)) size_t ZSTD_compressBlock_fast_##dictMode##_##mls##_##step(                                      \
            ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],                    \
            void const* src, size_t srcSize)                                                       \
    {                                                                                              \
        return ZSTD_compressBlock_fast_##dictMode##_generic(ms, seqStore, rep, src, srcSize, mls, step); \
    }

ZSTD_GEN_FAST_FN(noDict, 4, 1)
ZSTD_GEN_FAST_FN(noDict, 5, 1)
ZSTD_GEN_FAST_FN(noDict, 6, 1)
ZSTD_GEN_FAST_FN(noDict, 7, 1)

ZSTD_GEN_FAST_FN(noDict, 4, 0)
ZSTD_GEN_FAST_FN(noDict, 5, 0)
ZSTD_GEN_FAST_FN(noDict, 6, 0)
ZSTD_GEN_FAST_FN(noDict, 7, 0)

size_t ZSTD_compressBlock_fast(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize)
{
    U32 const mls = ms->cParams.minMatch;
    assert(ms->dictMatchState == NULL);
    if (ms->cParams.targetLength > 1) {
        switch(mls)
        {
        default: /* includes case 3 */
        case 4 :
            return ZSTD_compressBlock_fast_noDict_4_1(ms, seqStore, rep, src, srcSize);
        case 5 :
            return ZSTD_compressBlock_fast_noDict_5_1(ms, seqStore, rep, src, srcSize);
        case 6 :
            return ZSTD_compressBlock_fast_noDict_6_1(ms, seqStore, rep, src, srcSize);
        case 7 :
            return ZSTD_compressBlock_fast_noDict_7_1(ms, seqStore, rep, src, srcSize);
        }
    } else {
        switch(mls)
        {
        default: /* includes case 3 */
        case 4 :
            return ZSTD_compressBlock_fast_noDict_4_0(ms, seqStore, rep, src, srcSize);
        case 5 :
            return ZSTD_compressBlock_fast_noDict_5_0(ms, seqStore, rep, src, srcSize);
        case 6 :
            return ZSTD_compressBlock_fast_noDict_6_0(ms, seqStore, rep, src, srcSize);
        case 7 :
            return ZSTD_compressBlock_fast_noDict_7_0(ms, seqStore, rep, src, srcSize);
        }

    }
}

#ifdef ZSTD_DICT_USES
static uint32_t dictPosUses[128 * 1024] = {};
static uint32_t sortedDictPos[128 * 1024] = {};

static void markDictIndexUsed(U32 dictStartIndex, U32 dictIndex)
{
    ++dictPosUses[dictIndex - dictStartIndex];
}

static void markRepIndexUsed(U32 dictStartIndex, U32 index, U32 dictIndexDelta, U32 prefixStartIndex)
{
    if (index < prefixStartIndex)
        ++dictPosUses[(index - dictIndexDelta) - dictStartIndex];
}

static int orderDictPos(void const* lhs, void const* rhs)
{
    uint32_t const lhsIndex = *(uint32_t const*)lhs;
    uint32_t const rhsIndex = *(uint32_t const*)rhs;
    uint32_t const lhsCount = dictPosUses[lhsIndex];
    uint32_t const rhsCount = dictPosUses[rhsIndex];
    if (lhsCount < rhsCount)
        return -1;
    if (lhsCount > rhsCount)
        return 1;
    else if (lhsIndex < rhsIndex)
        return -1;
    else {
        assert(lhsIndex > rhsIndex);
        return 1;
    }
}

#include <stdlib.h>
#include <stdio.h>

void clearDictUses(void)
{
    DEBUGLOG(2, "clear dict pos uses");
    ZSTD_memset(dictPosUses, 0, sizeof(dictPosUses));
}

void writeDictUses(ZSTD_matchState_t const* dms) {
    U32 const dictStartIndex = dms->window.dictLimit;
    BYTE const* const dictStart = dms->window.base + dictStartIndex;
    BYTE const* const dictEnd = dms->window.nextSrc;
    U32 const dictSize = (U32)(dictEnd - dictStart);
    DEBUGLOG(2, "hashTableSize = %u", 1u << dms->cParams.hashLog);
    DEBUGLOG(2, "dictStartIndex = %u", dictStartIndex);
    DEBUGLOG(2, "dictSize = %u", dictSize);
    {
        U32 used = 0;
        U32 i;
        for (i = 0; i < dictSize; ++i) {
            if (dictPosUses[i] > 0)
                ++used;
            sortedDictPos[i] = i;
        }
        DEBUGLOG(2, "usedPos = %u", used);
        qsort(sortedDictPos, dictSize, sizeof(sortedDictPos[0]), &orderDictPos);
        assert(dictPosUses[sortedDictPos[dictSize - 1]] > dictPosUses[sortedDictPos[0]]);
    }
    {
        U32 startPos;
        for (startPos = 0; dictPosUses[sortedDictPos[startPos]] == 0; ++startPos) {}
        FILE* f = fopen("dict-uses.u32", "wb");
        DEBUGLOG(2, "startPos = %u", startPos);
        size_t const w = fwrite(sortedDictPos + startPos, sizeof(sortedDictPos[0]), dictSize - startPos, f);
        assert(f);
        assert(w == dictSize - startPos);
        fclose(f);
    }
}
#else
void clearDictUses(void)
{

}
void writeDictUses(ZSTD_matchState_t const* dms) {
    (void)dms;
}
static void markDictIndexUsed(U32 dictStartIndex, U32 dictIndex)
{
    (void)dictStartIndex, (void)dictIndex;
}
static void markRepIndexUsed(U32 dictStartIndex, U32 index, U32 dictIndexDelta, U32 prefixStartIndex)
{
    (void)dictStartIndex, (void)index, (void)dictIndexDelta, (void)prefixStartIndex;
}
#endif

FORCE_INLINE_TEMPLATE
size_t ZSTD_compressBlock_fast_dictMatchState_generic(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize, U32 const mls, U32 const hasStep)
{
    const ZSTD_compressionParameters* const cParams = &ms->cParams;
    U32* const hashTable = ms->hashTable;
    U32 const hlog = cParams->hashLog;
    /* support stepSize of 0 */
    U32 const stepSize = cParams->targetLength + !(cParams->targetLength);
    const BYTE* const base = ms->window.base;
    const BYTE* const istart = (const BYTE*)src;
    const BYTE* ip = istart;
    const BYTE* anchor = istart;
    const U32   prefixStartIndex = ms->window.dictLimit;
    const BYTE* const prefixStart = base + prefixStartIndex;
    const BYTE* const iend = istart + srcSize;
    const BYTE* const ilimit = iend - HASH_READ_SIZE;
    U32 offset_1=rep[0], offset_2=rep[1];
    U32 offsetSaved = 0;

    const ZSTD_matchState_t* const dms = ms->dictMatchState;
    const ZSTD_compressionParameters* const dictCParams = &dms->cParams ;
    const U32* const dictHashTable = dms->hashTable;
    const U32 dictStartIndex       = dms->window.dictLimit;
    const BYTE* const dictBase     = dms->window.base;
    const BYTE* const dictStart    = dictBase + dictStartIndex;
    const BYTE* const dictEnd      = dms->window.nextSrc;
    const U32 dictIndexDelta       = prefixStartIndex - (U32)(dictEnd - dictBase);
    const U32 dictAndPrefixLength  = (U32)(ip - prefixStart + dictEnd - dictStart);
    const U32 dictHLog             = dictCParams->hashLog;

    /* if a dictionary is still attached, it necessarily means that
     * it is within window size. So we just check it. */
    const U32 maxDistance = 1U << cParams->windowLog;
    const U32 endIndex = (U32)((size_t)(ip - base) + srcSize);
    assert(endIndex - prefixStartIndex <= maxDistance);
    (void)maxDistance; (void)endIndex;   /* these variables are not used when assert() is disabled */

    (void)hasStep; /* not currently specialized on whether it's accelerated */

    /* ensure there will be no underflow
     * when translating a dict index into a local index */
    assert(prefixStartIndex >= (U32)(dictEnd - dictBase));

    /* init */
    DEBUGLOG(5, "ZSTD_compressBlock_fast_dictMatchState_generic");
    ip += (dictAndPrefixLength == 0);
    /* dictMatchState repCode checks don't currently handle repCode == 0
     * disabling. */
    assert(offset_1 <= dictAndPrefixLength);
    assert(offset_2 <= dictAndPrefixLength);

    /* Main Search Loop */
    while (ip < ilimit) {   /* < instead of <=, because repcode check at (ip+1) */
        size_t mLength;
        size_t const h = ZSTD_hashPtr(ip, hlog, mls);
        U32 const curr = (U32)(ip-base);
        U32 const matchIndex = hashTable[h];
        const BYTE* match = base + matchIndex;
        const U32 repIndex = curr + 1 - offset_1;
        const BYTE* repMatch = (repIndex < prefixStartIndex) ?
                               dictBase + (repIndex - dictIndexDelta) :
                               base + repIndex;
        hashTable[h] = curr;   /* update hash table */

        if ( ((U32)((prefixStartIndex-1) - repIndex) >= 3) /* intentional underflow : ensure repIndex isn't overlapping dict + prefix */
          && (MEM_read32(repMatch) == MEM_read32(ip+1)) ) {
            const BYTE* const repMatchEnd = repIndex < prefixStartIndex ? dictEnd : iend;
            mLength = ZSTD_count_2segments(ip+1+4, repMatch+4, iend, repMatchEnd, prefixStart) + 4;
            ip++;
            markRepIndexUsed(dictStartIndex, repIndex, dictIndexDelta, prefixStartIndex);
            ZSTD_storeSeq(seqStore, (size_t)(ip-anchor), anchor, iend, REPCODE1_TO_OFFBASE, mLength);
        } else if ( (matchIndex <= prefixStartIndex) ) {
            size_t const dictHash = ZSTD_hashPtr(ip, dictHLog, mls);
            U32 const dictMatchIndex = dictHashTable[dictHash];
            const BYTE* dictMatch = dictBase + dictMatchIndex;
            if (dictMatchIndex <= dictStartIndex ||
                MEM_read32(dictMatch) != MEM_read32(ip)) {
                assert(stepSize >= 1);
                ip += ((ip-anchor) >> kSearchStrength) + stepSize;
                continue;
            } else {
                /* found a dict match */
                U32 const offset = (U32)(curr-dictMatchIndex-dictIndexDelta);
                mLength = ZSTD_count_2segments(ip+4, dictMatch+4, iend, dictEnd, prefixStart) + 4;
                while (((ip>anchor) & (dictMatch>dictStart))
                     && (ip[-1] == dictMatch[-1])) {
                    ip--; dictMatch--; mLength++;
                } /* catch up */
                offset_2 = offset_1;
                offset_1 = offset;
                markDictIndexUsed(dictStartIndex, dictMatchIndex);
                ZSTD_storeSeq(seqStore, (size_t)(ip-anchor), anchor, iend, OFFSET_TO_OFFBASE(offset), mLength);
            }
        } else if (MEM_read32(match) != MEM_read32(ip)) {
            /* it's not a match, and we're not going to check the dictionary */
            assert(stepSize >= 1);
            ip += ((ip-anchor) >> kSearchStrength) + stepSize;
            continue;
        } else {
            /* found a regular match */
            U32 const offset = (U32)(ip-match);
            mLength = ZSTD_count(ip+4, match+4, iend) + 4;
            while (((ip>anchor) & (match>prefixStart))
                 && (ip[-1] == match[-1])) { ip--; match--; mLength++; } /* catch up */
            offset_2 = offset_1;
            offset_1 = offset;
            ZSTD_storeSeq(seqStore, (size_t)(ip-anchor), anchor, iend, OFFSET_TO_OFFBASE(offset), mLength);
        }

        /* match found */
        ip += mLength;
        anchor = ip;

        if (ip <= ilimit) {
            /* Fill Table */
            assert(base+curr+2 > istart);  /* check base overflow */
            hashTable[ZSTD_hashPtr(base+curr+2, hlog, mls)] = curr+2;  /* here because curr+2 could be > iend-8 */
            hashTable[ZSTD_hashPtr(ip-2, hlog, mls)] = (U32)(ip-2-base);

            /* check immediate repcode */
            while (ip <= ilimit) {
                U32 const current2 = (U32)(ip-base);
                U32 const repIndex2 = current2 - offset_2;
                const BYTE* repMatch2 = repIndex2 < prefixStartIndex ?
                        dictBase - dictIndexDelta + repIndex2 :
                        base + repIndex2;
                if ( ((U32)((prefixStartIndex-1) - (U32)repIndex2) >= 3 /* intentional overflow */)
                   && (MEM_read32(repMatch2) == MEM_read32(ip)) ) {
                    const BYTE* const repEnd2 = repIndex2 < prefixStartIndex ? dictEnd : iend;
                    size_t const repLength2 = ZSTD_count_2segments(ip+4, repMatch2+4, iend, repEnd2, prefixStart) + 4;
                    U32 tmpOffset = offset_2; offset_2 = offset_1; offset_1 = tmpOffset;   /* swap offset_2 <=> offset_1 */
                    markRepIndexUsed(dictStartIndex, repIndex2, dictIndexDelta, prefixStartIndex);
                    ZSTD_storeSeq(seqStore, 0, anchor, iend, REPCODE1_TO_OFFBASE, repLength2);
                    hashTable[ZSTD_hashPtr(ip, hlog, mls)] = current2;
                    ip += repLength2;
                    anchor = ip;
                    continue;
                }
                break;
            }
        }
    }

    /* save reps for next block */
    rep[0] = offset_1 ? offset_1 : offsetSaved;
    rep[1] = offset_2 ? offset_2 : offsetSaved;

    /* Return the last literals size */
    return (size_t)(iend - anchor);
}


ZSTD_GEN_FAST_FN(dictMatchState, 4, 0)
ZSTD_GEN_FAST_FN(dictMatchState, 5, 0)
ZSTD_GEN_FAST_FN(dictMatchState, 6, 0)
ZSTD_GEN_FAST_FN(dictMatchState, 7, 0)

FORCE_INLINE_TEMPLATE
size_t ZSTD_compressBlock_fast_dds_generic(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize, U32 const mls, U32 const hasStep)
{
    const ZSTD_compressionParameters* const cParams = &ms->cParams;
    U32* const hashTable = ms->hashTable;
    U32 const hlog = cParams->hashLog;
    /* support stepSize of 0 */
    U32 const stepSize = cParams->targetLength + !(cParams->targetLength);
    const BYTE* const base = ms->window.base;
    const BYTE* const istart = (const BYTE*)src;
    const BYTE* ip = istart;
    const BYTE* anchor = istart;
    const U32   prefixStartIndex = ms->window.dictLimit;
    const BYTE* const prefixStart = base + prefixStartIndex;
    const BYTE* const iend = istart + srcSize;
    const BYTE* const ilimit = iend - HASH_READ_SIZE;
    U32 offset_1=rep[0], offset_2=rep[1];
    U32 offsetSaved = 0;

    const ZSTD_matchState_t* const dms = ms->dictMatchState;
    const ZSTD_compressionParameters* const dictCParams = &dms->cParams ;
    const U16* const dictHashTable = (U16 const*)(void const*)dms->hashTable;
    const U16 dictStartIndex       = dms->window.dictLimit;
    const BYTE* const dictBase     = dms->window.base;
    const BYTE* const dictStart    = dictBase + dictStartIndex;
    const BYTE* const dictEnd      = dms->window.nextSrc;
    const U32 dictIndexDelta       = prefixStartIndex - (U32)(dictEnd - dictBase);
    const U32 dictAndPrefixLength  = (U32)(ip - prefixStart + dictEnd - dictStart);
    const U32 dictHLog             = (dictCParams->hashLog + 2) - 5 + 8;

    /* if a dictionary is still attached, it necessarily means that
     * it is within window size. So we just check it. */
    const U32 maxDistance = 1U << cParams->windowLog;
    const U32 endIndex = (U32)((size_t)(ip - base) + srcSize);
    assert(endIndex - prefixStartIndex <= maxDistance);
    (void)maxDistance; (void)endIndex;   /* these variables are not used when assert() is disabled */

    (void)hasStep; /* not currently specialized on whether it's accelerated */

    /* ensure there will be no underflow
     * when translating a dict index into a local index */
    assert(prefixStartIndex >= (U32)(dictEnd - dictBase));

    /* init */
    DEBUGLOG(5, "ZSTD_compressBlock_fast_dictMatchState_generic");
    ip += (dictAndPrefixLength == 0);
    /* dictMatchState repCode checks don't currently handle repCode == 0
     * disabling. */
    assert(offset_1 <= dictAndPrefixLength);
    assert(offset_2 <= dictAndPrefixLength);

    /* Main Search Loop */
    while (ip < ilimit) {   /* < instead of <=, because repcode check at (ip+1) */
        size_t mLength;
        size_t const h = ZSTD_hashPtr(ip, hlog, mls);
        U32 const curr = (U32)(ip-base);
        U32 const matchIndex = hashTable[h];
        const BYTE* match = base + matchIndex;
        const U32 repIndex = curr + 1 - offset_1;
        const BYTE* repMatch = (repIndex < prefixStartIndex) ?
                               dictBase + (repIndex - dictIndexDelta) :
                               base + repIndex;
        hashTable[h] = curr;   /* update hash table */

        if ( ((U32)((prefixStartIndex-1) - repIndex) >= 3) /* intentional underflow : ensure repIndex isn't overlapping dict + prefix */
          && (MEM_read32(repMatch) == MEM_read32(ip+1)) ) {
            const BYTE* const repMatchEnd = repIndex < prefixStartIndex ? dictEnd : iend;
            mLength = ZSTD_count_2segments(ip+1+4, repMatch+4, iend, repMatchEnd, prefixStart) + 4;
            ip++;
            markRepIndexUsed(dictStartIndex, repIndex, dictIndexDelta, prefixStartIndex);
            ZSTD_storeSeq(seqStore, (size_t)(ip-anchor), anchor, iend, REPCODE1_TO_OFFBASE, mLength);
        } else if ( (matchIndex <= prefixStartIndex) ) {
            size_t const dictHash = ZSTD_hashPtr(ip, dictHLog, mls);
            size_t const dictRow  = (dictHash >> 8) << 4;
            size_t const dictTag  = dictHash & 0xFF;
            ZSTD_VecMask const tagMask = ZSTD_getTagMask(dictHashTable + dictRow, dictTag);
            if (tagMask == 0) {
                assert(stepSize >= 1);
                ip += ((ip-anchor) >> kSearchStrength) + stepSize;
                continue;
            } else {
                U32 const dictRowIdx = ZSTD_VecMask_next(tagMask);
                U16 const dictMatchIndex = dictHashTable[dictRow + 5 + dictRowIdx];
                const BYTE* dictMatch = dictBase + dictMatchIndex;
                assert(dictMatchIndex <= curr);
                assert((U32)(dictEnd - dictBase) <= 65536);
                assert(dictMatchIndex <= (U32)(dictEnd - dictBase));
                assert(dictMatchIndex >= dictStartIndex);
                if (MEM_read32(dictMatch) != MEM_read32(ip)) {
                    assert(stepSize >= 1);
                    ip += ((ip-anchor) >> kSearchStrength) + stepSize;
                    continue;
                } else {
                    /* found a dict match */
                    U32 const offset = (U32)((curr-dictMatchIndex)-dictIndexDelta);
                    mLength = ZSTD_count_2segments(ip+4, dictMatch+4, iend, dictEnd, prefixStart) + 4;
                    while (((ip>anchor) & (dictMatch>dictStart))
                        && (ip[-1] == dictMatch[-1])) {
                        ip--; dictMatch--; mLength++;
                    } /* catch up */
                    offset_2 = offset_1;
                    offset_1 = offset;
                    markDictIndexUsed(dictStartIndex, dictMatchIndex);
                    ZSTD_storeSeq(seqStore, (size_t)(ip-anchor), anchor, iend, OFFSET_TO_OFFBASE(offset), mLength);
                }
            }
        } else if (MEM_read32(match) != MEM_read32(ip)) {
            /* it's not a match, and we're not going to check the dictionary */
            assert(stepSize >= 1);
            ip += ((ip-anchor) >> kSearchStrength) + stepSize;
            continue;
        } else {
            /* found a regular match */
            U32 const offset = (U32)(ip-match);
            mLength = ZSTD_count(ip+4, match+4, iend) + 4;
            while (((ip>anchor) & (match>prefixStart))
                 && (ip[-1] == match[-1])) { ip--; match--; mLength++; } /* catch up */
            offset_2 = offset_1;
            offset_1 = offset;
            ZSTD_storeSeq(seqStore, (size_t)(ip-anchor), anchor, iend, OFFSET_TO_OFFBASE(offset), mLength);
        }

        /* match found */
        ip += mLength;
        anchor = ip;

        if (ip <= ilimit) {
            /* Fill Table */
            assert(base+curr+2 > istart);  /* check base overflow */
            hashTable[ZSTD_hashPtr(base+curr+2, hlog, mls)] = curr+2;  /* here because curr+2 could be > iend-8 */
            hashTable[ZSTD_hashPtr(ip-2, hlog, mls)] = (U32)(ip-2-base);

            /* check immediate repcode */
            while (ip <= ilimit) {
                U32 const current2 = (U32)(ip-base);
                U32 const repIndex2 = current2 - offset_2;
                const BYTE* repMatch2 = repIndex2 < prefixStartIndex ?
                        dictBase - dictIndexDelta + repIndex2 :
                        base + repIndex2;
                if ( ((U32)((prefixStartIndex-1) - (U32)repIndex2) >= 3 /* intentional overflow */)
                   && (MEM_read32(repMatch2) == MEM_read32(ip)) ) {
                    const BYTE* const repEnd2 = repIndex2 < prefixStartIndex ? dictEnd : iend;
                    size_t const repLength2 = ZSTD_count_2segments(ip+4, repMatch2+4, iend, repEnd2, prefixStart) + 4;
                    U32 tmpOffset = offset_2; offset_2 = offset_1; offset_1 = tmpOffset;   /* swap offset_2 <=> offset_1 */
                    markRepIndexUsed(dictStartIndex, repIndex2, dictIndexDelta, prefixStartIndex);
                    ZSTD_storeSeq(seqStore, 0, anchor, iend, REPCODE1_TO_OFFBASE, repLength2);
                    hashTable[ZSTD_hashPtr(ip, hlog, mls)] = current2;
                    ip += repLength2;
                    anchor = ip;
                    continue;
                }
                break;
            }
        }
    }

    /* save reps for next block */
    rep[0] = offset_1 ? offset_1 : offsetSaved;
    rep[1] = offset_2 ? offset_2 : offsetSaved;

    /* Return the last literals size */
    return (size_t)(iend - anchor);
}

ZSTD_GEN_FAST_FN(dds, 4, 0)
ZSTD_GEN_FAST_FN(dds, 5, 0)
ZSTD_GEN_FAST_FN(dds, 6, 0)
ZSTD_GEN_FAST_FN(dds, 7, 0)

static size_t ZSTD_compressBlock_fast_dds(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize)
{
    U32 const mls = ms->cParams.minMatch;
    assert(ms->dictMatchState != NULL);
    if (ms->dictIsCold)
        PREFETCH_AREA(ms->dictMatchState->hashTable, 4 << ms->dictMatchState->cParams.hashLog);
    switch(mls)
    {
    default: /* includes case 3 */
    case 4 :
        return ZSTD_compressBlock_fast_dds_4_0(ms, seqStore, rep, src, srcSize);
    case 5 :
        return ZSTD_compressBlock_fast_dds_5_0(ms, seqStore, rep, src, srcSize);
    case 6 :
        return ZSTD_compressBlock_fast_dds_6_0(ms, seqStore, rep, src, srcSize);
    case 7 :
        return ZSTD_compressBlock_fast_dds_7_0(ms, seqStore, rep, src, srcSize);
    }
}

size_t ZSTD_compressBlock_fast_dictMatchState(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize)
{
    U32 const mls = ms->cParams.minMatch;
    if (ms->dictMatchState->ddsFast) {
        return ZSTD_compressBlock_fast_dds(ms, seqStore, rep, src, srcSize);
    }
    assert(ms->dictMatchState != NULL);
    switch(mls)
    {
    default: /* includes case 3 */
    case 4 :
        return ZSTD_compressBlock_fast_dictMatchState_4_0(ms, seqStore, rep, src, srcSize);
    case 5 :
        return ZSTD_compressBlock_fast_dictMatchState_5_0(ms, seqStore, rep, src, srcSize);
    case 6 :
        return ZSTD_compressBlock_fast_dictMatchState_6_0(ms, seqStore, rep, src, srcSize);
    case 7 :
        return ZSTD_compressBlock_fast_dictMatchState_7_0(ms, seqStore, rep, src, srcSize);
    }
}


static size_t ZSTD_compressBlock_fast_extDict_generic(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize, U32 const mls, U32 const hasStep)
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
    const U32   endIndex = (U32)((size_t)(istart - base) + srcSize);
    const U32   lowLimit = ZSTD_getLowestMatchIndex(ms, endIndex, cParams->windowLog);
    const U32   dictStartIndex = lowLimit;
    const BYTE* const dictStart = dictBase + dictStartIndex;
    const U32   dictLimit = ms->window.dictLimit;
    const U32   prefixStartIndex = dictLimit < lowLimit ? lowLimit : dictLimit;
    const BYTE* const prefixStart = base + prefixStartIndex;
    const BYTE* const dictEnd = dictBase + prefixStartIndex;
    const BYTE* const iend = istart + srcSize;
    const BYTE* const ilimit = iend - 8;
    U32 offset_1=rep[0], offset_2=rep[1];

    (void)hasStep; /* not currently specialized on whether it's accelerated */

    DEBUGLOG(5, "ZSTD_compressBlock_fast_extDict_generic (offset_1=%u)", offset_1);

    /* switch to "regular" variant if extDict is invalidated due to maxDistance */
    if (prefixStartIndex == dictStartIndex)
        return ZSTD_compressBlock_fast(ms, seqStore, rep, src, srcSize);

    /* Search Loop */
    while (ip < ilimit) {  /* < instead of <=, because (ip+1) */
        const size_t h = ZSTD_hashPtr(ip, hlog, mls);
        const U32    matchIndex = hashTable[h];
        const BYTE* const matchBase = matchIndex < prefixStartIndex ? dictBase : base;
        const BYTE*  match = matchBase + matchIndex;
        const U32    curr = (U32)(ip-base);
        const U32    repIndex = curr + 1 - offset_1;
        const BYTE* const repBase = repIndex < prefixStartIndex ? dictBase : base;
        const BYTE* const repMatch = repBase + repIndex;
        hashTable[h] = curr;   /* update hash table */
        DEBUGLOG(7, "offset_1 = %u , curr = %u", offset_1, curr);

        if ( ( ((U32)((prefixStartIndex-1) - repIndex) >= 3) /* intentional underflow */
             & (offset_1 <= curr+1 - dictStartIndex) ) /* note: we are searching at curr+1 */
           && (MEM_read32(repMatch) == MEM_read32(ip+1)) ) {
            const BYTE* const repMatchEnd = repIndex < prefixStartIndex ? dictEnd : iend;
            size_t const rLength = ZSTD_count_2segments(ip+1 +4, repMatch +4, iend, repMatchEnd, prefixStart) + 4;
            ip++;
            ZSTD_storeSeq(seqStore, (size_t)(ip-anchor), anchor, iend, REPCODE1_TO_OFFBASE, rLength);
            ip += rLength;
            anchor = ip;
        } else {
            if ( (matchIndex < dictStartIndex) ||
                 (MEM_read32(match) != MEM_read32(ip)) ) {
                assert(stepSize >= 1);
                ip += ((ip-anchor) >> kSearchStrength) + stepSize;
                continue;
            }
            {   const BYTE* const matchEnd = matchIndex < prefixStartIndex ? dictEnd : iend;
                const BYTE* const lowMatchPtr = matchIndex < prefixStartIndex ? dictStart : prefixStart;
                U32 const offset = curr - matchIndex;
                size_t mLength = ZSTD_count_2segments(ip+4, match+4, iend, matchEnd, prefixStart) + 4;
                while (((ip>anchor) & (match>lowMatchPtr)) && (ip[-1] == match[-1])) { ip--; match--; mLength++; }   /* catch up */
                offset_2 = offset_1; offset_1 = offset;  /* update offset history */
                ZSTD_storeSeq(seqStore, (size_t)(ip-anchor), anchor, iend, OFFSET_TO_OFFBASE(offset), mLength);
                ip += mLength;
                anchor = ip;
        }   }

        if (ip <= ilimit) {
            /* Fill Table */
            hashTable[ZSTD_hashPtr(base+curr+2, hlog, mls)] = curr+2;
            hashTable[ZSTD_hashPtr(ip-2, hlog, mls)] = (U32)(ip-2-base);
            /* check immediate repcode */
            while (ip <= ilimit) {
                U32 const current2 = (U32)(ip-base);
                U32 const repIndex2 = current2 - offset_2;
                const BYTE* const repMatch2 = repIndex2 < prefixStartIndex ? dictBase + repIndex2 : base + repIndex2;
                if ( (((U32)((prefixStartIndex-1) - repIndex2) >= 3) & (offset_2 <= curr - dictStartIndex))  /* intentional overflow */
                   && (MEM_read32(repMatch2) == MEM_read32(ip)) ) {
                    const BYTE* const repEnd2 = repIndex2 < prefixStartIndex ? dictEnd : iend;
                    size_t const repLength2 = ZSTD_count_2segments(ip+4, repMatch2+4, iend, repEnd2, prefixStart) + 4;
                    { U32 const tmpOffset = offset_2; offset_2 = offset_1; offset_1 = tmpOffset; }  /* swap offset_2 <=> offset_1 */
                    ZSTD_storeSeq(seqStore, 0 /*litlen*/, anchor, iend, REPCODE1_TO_OFFBASE, repLength2);
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
    return (size_t)(iend - anchor);
}

ZSTD_GEN_FAST_FN(extDict, 4, 0)
ZSTD_GEN_FAST_FN(extDict, 5, 0)
ZSTD_GEN_FAST_FN(extDict, 6, 0)
ZSTD_GEN_FAST_FN(extDict, 7, 0)

size_t ZSTD_compressBlock_fast_extDict(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize)
{
    U32 const mls = ms->cParams.minMatch;
    switch(mls)
    {
    default: /* includes case 3 */
    case 4 :
        return ZSTD_compressBlock_fast_extDict_4_0(ms, seqStore, rep, src, srcSize);
    case 5 :
        return ZSTD_compressBlock_fast_extDict_5_0(ms, seqStore, rep, src, srcSize);
    case 6 :
        return ZSTD_compressBlock_fast_extDict_6_0(ms, seqStore, rep, src, srcSize);
    case 7 :
        return ZSTD_compressBlock_fast_extDict_7_0(ms, seqStore, rep, src, srcSize);
    }
}
