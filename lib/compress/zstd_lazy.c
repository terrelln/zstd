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
#include "zstd_lazy.h"


/*-*************************************
*  Binary Tree search
***************************************/

static void
ZSTD_updateDUBT(ZSTD_matchState_t* ms,
                const BYTE* ip, const BYTE* iend,
                U32 mls)
{
    const ZSTD_compressionParameters* const cParams = &ms->cParams;
    U32* const hashTable = ms->hashTable;
    U32  const hashLog = cParams->hashLog;

    U32* const bt = ms->chainTable;
    U32  const btLog  = cParams->chainLog - 1;
    U32  const btMask = (1 << btLog) - 1;

    const BYTE* const base = ms->window.base;
    U32 const target = (U32)(ip - base);
    U32 idx = ms->nextToUpdate;

    if (idx != target)
        DEBUGLOG(7, "ZSTD_updateDUBT, from %u to %u (dictLimit:%u)",
                    idx, target, ms->window.dictLimit);
    assert(ip + 8 <= iend);   /* condition for ZSTD_hashPtr */
    (void)iend;

    assert(idx >= ms->window.dictLimit);   /* condition for valid base+idx */
    for ( ; idx < target ; idx++) {
        size_t const h  = ZSTD_hashPtr(base + idx, hashLog, mls);   /* assumption : ip + 8 <= iend */
        U32    const matchIndex = hashTable[h];

        U32*   const nextCandidatePtr = bt + 2*(idx&btMask);
        U32*   const sortMarkPtr  = nextCandidatePtr + 1;

        DEBUGLOG(8, "ZSTD_updateDUBT: insert %u", idx);
        hashTable[h] = idx;   /* Update Hash Table */
        *nextCandidatePtr = matchIndex;   /* update BT like a chain */
        *sortMarkPtr = ZSTD_DUBT_UNSORTED_MARK;
    }
    ms->nextToUpdate = target;
}


/** ZSTD_insertDUBT1() :
 *  sort one already inserted but unsorted position
 *  assumption : current >= btlow == (current - btmask)
 *  doesn't fail */
static void
ZSTD_insertDUBT1(ZSTD_matchState_t* ms,
                 U32 current, const BYTE* inputEnd,
                 U32 nbCompares, U32 btLow,
                 const ZSTD_dictMode_e dictMode)
{
    const ZSTD_compressionParameters* const cParams = &ms->cParams;
    U32* const bt = ms->chainTable;
    U32  const btLog  = cParams->chainLog - 1;
    U32  const btMask = (1 << btLog) - 1;
    size_t commonLengthSmaller=0, commonLengthLarger=0;
    const BYTE* const base = ms->window.base;
    const BYTE* const dictBase = ms->window.dictBase;
    const U32 dictLimit = ms->window.dictLimit;
    const BYTE* const ip = (current>=dictLimit) ? base + current : dictBase + current;
    const BYTE* const iend = (current>=dictLimit) ? inputEnd : dictBase + dictLimit;
    const BYTE* const dictEnd = dictBase + dictLimit;
    const BYTE* const prefixStart = base + dictLimit;
    const BYTE* match;
    U32* smallerPtr = bt + 2*(current&btMask);
    U32* largerPtr  = smallerPtr + 1;
    U32 matchIndex = *smallerPtr;   /* this candidate is unsorted : next sorted candidate is reached through *smallerPtr, while *largerPtr contains previous unsorted candidate (which is already saved and can be overwritten) */
    U32 dummy32;   /* to be nullified at the end */
    U32 const windowLow = ms->window.lowLimit;

    DEBUGLOG(8, "ZSTD_insertDUBT1(%u) (dictLimit=%u, lowLimit=%u)",
                current, dictLimit, windowLow);
    assert(current >= btLow);
    assert(ip < iend);   /* condition for ZSTD_count */

    while (nbCompares-- && (matchIndex > windowLow)) {
        U32* const nextPtr = bt + 2*(matchIndex & btMask);
        size_t matchLength = MIN(commonLengthSmaller, commonLengthLarger);   /* guaranteed minimum nb of common bytes */
        assert(matchIndex < current);
        /* note : all candidates are now supposed sorted,
         * but it's still possible to have nextPtr[1] == ZSTD_DUBT_UNSORTED_MARK
         * when a real index has the same value as ZSTD_DUBT_UNSORTED_MARK */

        if ( (dictMode != ZSTD_extDict)
          || (matchIndex+matchLength >= dictLimit)  /* both in current segment*/
          || (current < dictLimit) /* both in extDict */) {
            const BYTE* const mBase = ( (dictMode != ZSTD_extDict)
                                     || (matchIndex+matchLength >= dictLimit)) ?
                                        base : dictBase;
            assert( (matchIndex+matchLength >= dictLimit)   /* might be wrong if extDict is incorrectly set to 0 */
                 || (current < dictLimit) );
            match = mBase + matchIndex;
            matchLength += ZSTD_count(ip+matchLength, match+matchLength, iend);
        } else {
            match = dictBase + matchIndex;
            matchLength += ZSTD_count_2segments(ip+matchLength, match+matchLength, iend, dictEnd, prefixStart);
            if (matchIndex+matchLength >= dictLimit)
                match = base + matchIndex;   /* preparation for next read of match[matchLength] */
        }

        DEBUGLOG(8, "ZSTD_insertDUBT1: comparing %u with %u : found %u common bytes ",
                    current, matchIndex, (U32)matchLength);

        if (ip+matchLength == iend) {   /* equal : no way to know if inf or sup */
            break;   /* drop , to guarantee consistency ; miss a bit of compression, but other solutions can corrupt tree */
        }

        if (match[matchLength] < ip[matchLength]) {  /* necessarily within buffer */
            /* match is smaller than current */
            *smallerPtr = matchIndex;             /* update smaller idx */
            commonLengthSmaller = matchLength;    /* all smaller will now have at least this guaranteed common length */
            if (matchIndex <= btLow) { smallerPtr=&dummy32; break; }   /* beyond tree size, stop searching */
            DEBUGLOG(8, "ZSTD_insertDUBT1: %u (>btLow=%u) is smaller : next => %u",
                        matchIndex, btLow, nextPtr[1]);
            smallerPtr = nextPtr+1;               /* new "candidate" => larger than match, which was smaller than target */
            matchIndex = nextPtr[1];              /* new matchIndex, larger than previous and closer to current */
        } else {
            /* match is larger than current */
            *largerPtr = matchIndex;
            commonLengthLarger = matchLength;
            if (matchIndex <= btLow) { largerPtr=&dummy32; break; }   /* beyond tree size, stop searching */
            DEBUGLOG(8, "ZSTD_insertDUBT1: %u (>btLow=%u) is larger => %u",
                        matchIndex, btLow, nextPtr[0]);
            largerPtr = nextPtr;
            matchIndex = nextPtr[0];
    }   }

    *smallerPtr = *largerPtr = 0;
}


static size_t
ZSTD_DUBT_findBetterDictMatch (
        ZSTD_matchState_t* ms,
        const BYTE* const ip, const BYTE* const iend,
        size_t* offsetPtr,
        size_t bestLength,
        U32 nbCompares,
        U32 const mls,
        const ZSTD_dictMode_e dictMode)
{
    const ZSTD_matchState_t * const dms = ms->dictMatchState;
    const ZSTD_compressionParameters* const dmsCParams = &dms->cParams;
    const U32 * const dictHashTable = dms->hashTable;
    U32         const hashLog = dmsCParams->hashLog;
    size_t      const h  = ZSTD_hashPtr(ip, hashLog, mls);
    U32               dictMatchIndex = dictHashTable[h];

    const BYTE* const base = ms->window.base;
    const BYTE* const prefixStart = base + ms->window.dictLimit;
    U32         const current = (U32)(ip-base);
    const BYTE* const dictBase = dms->window.base;
    const BYTE* const dictEnd = dms->window.nextSrc;
    U32         const dictHighLimit = (U32)(dms->window.nextSrc - dms->window.base);
    U32         const dictLowLimit = dms->window.lowLimit;
    U32         const dictIndexDelta = ms->window.lowLimit - dictHighLimit;

    U32*        const dictBt = dms->chainTable;
    U32         const btLog  = dmsCParams->chainLog - 1;
    U32         const btMask = (1 << btLog) - 1;
    U32         const btLow = (btMask >= dictHighLimit - dictLowLimit) ? dictLowLimit : dictHighLimit - btMask;

    size_t commonLengthSmaller=0, commonLengthLarger=0;

    (void)dictMode;
    assert(dictMode == ZSTD_dictMatchState);

    while (nbCompares-- && (dictMatchIndex > dictLowLimit)) {
        U32* const nextPtr = dictBt + 2*(dictMatchIndex & btMask);
        size_t matchLength = MIN(commonLengthSmaller, commonLengthLarger);   /* guaranteed minimum nb of common bytes */
        const BYTE* match = dictBase + dictMatchIndex;
        matchLength += ZSTD_count_2segments(ip+matchLength, match+matchLength, iend, dictEnd, prefixStart);
        if (dictMatchIndex+matchLength >= dictHighLimit)
            match = base + dictMatchIndex + dictIndexDelta;   /* to prepare for next usage of match[matchLength] */

        if (matchLength > bestLength) {
            U32 matchIndex = dictMatchIndex + dictIndexDelta;
            if ( (4*(int)(matchLength-bestLength)) > (int)(ZSTD_highbit32(current-matchIndex+1) - ZSTD_highbit32((U32)offsetPtr[0]+1)) ) {
                DEBUGLOG(9, "ZSTD_DUBT_findBetterDictMatch(%u) : found better match length %u -> %u and offsetCode %u -> %u (dictMatchIndex %u, matchIndex %u)",
                    current, (U32)bestLength, (U32)matchLength, (U32)*offsetPtr, ZSTD_REP_MOVE + current - matchIndex, dictMatchIndex, matchIndex);
                bestLength = matchLength, *offsetPtr = ZSTD_REP_MOVE + current - matchIndex;
            }
            if (ip+matchLength == iend) {   /* reached end of input : ip[matchLength] is not valid, no way to know if it's larger or smaller than match */
                break;   /* drop, to guarantee consistency (miss a little bit of compression) */
            }
        }

        if (match[matchLength] < ip[matchLength]) {
            if (dictMatchIndex <= btLow) { break; }   /* beyond tree size, stop the search */
            commonLengthSmaller = matchLength;    /* all smaller will now have at least this guaranteed common length */
            dictMatchIndex = nextPtr[1];              /* new matchIndex larger than previous (closer to current) */
        } else {
            /* match is larger than current */
            if (dictMatchIndex <= btLow) { break; }   /* beyond tree size, stop the search */
            commonLengthLarger = matchLength;
            dictMatchIndex = nextPtr[0];
        }
    }

    if (bestLength >= MINMATCH) {
        U32 const mIndex = current - ((U32)*offsetPtr - ZSTD_REP_MOVE); (void)mIndex;
        DEBUGLOG(8, "ZSTD_DUBT_findBetterDictMatch(%u) : found match of length %u and offsetCode %u (pos %u)",
                    current, (U32)bestLength, (U32)*offsetPtr, mIndex);
    }
    return bestLength;

}


static size_t
ZSTD_DUBT_findBestMatch(ZSTD_matchState_t* ms,
                        const BYTE* const ip, const BYTE* const iend,
                        size_t* offsetPtr,
                        U32 const mls,
                        const ZSTD_dictMode_e dictMode)
{
    const ZSTD_compressionParameters* const cParams = &ms->cParams;
    U32*   const hashTable = ms->hashTable;
    U32    const hashLog = cParams->hashLog;
    size_t const h  = ZSTD_hashPtr(ip, hashLog, mls);
    U32          matchIndex  = hashTable[h];

    const BYTE* const base = ms->window.base;
    U32    const current = (U32)(ip-base);
    U32    const windowLow = ms->window.lowLimit;

    U32*   const bt = ms->chainTable;
    U32    const btLog  = cParams->chainLog - 1;
    U32    const btMask = (1 << btLog) - 1;
    U32    const btLow = (btMask >= current) ? 0 : current - btMask;
    U32    const unsortLimit = MAX(btLow, windowLow);

    U32*         nextCandidate = bt + 2*(matchIndex&btMask);
    U32*         unsortedMark = bt + 2*(matchIndex&btMask) + 1;
    U32          nbCompares = 1U << cParams->searchLog;
    U32          nbCandidates = nbCompares;
    U32          previousCandidate = 0;

    DEBUGLOG(7, "ZSTD_DUBT_findBestMatch (%u) ", current);
    assert(ip <= iend-8);   /* required for h calculation */

    /* reach end of unsorted candidates list */
    while ( (matchIndex > unsortLimit)
         && (*unsortedMark == ZSTD_DUBT_UNSORTED_MARK)
         && (nbCandidates > 1) ) {
        DEBUGLOG(8, "ZSTD_DUBT_findBestMatch: candidate %u is unsorted",
                    matchIndex);
        *unsortedMark = previousCandidate;  /* the unsortedMark becomes a reversed chain, to move up back to original position */
        previousCandidate = matchIndex;
        matchIndex = *nextCandidate;
        nextCandidate = bt + 2*(matchIndex&btMask);
        unsortedMark = bt + 2*(matchIndex&btMask) + 1;
        nbCandidates --;
    }

    /* nullify last candidate if it's still unsorted
     * simplification, detrimental to compression ratio, beneficial for speed */
    if ( (matchIndex > unsortLimit)
      && (*unsortedMark==ZSTD_DUBT_UNSORTED_MARK) ) {
        DEBUGLOG(7, "ZSTD_DUBT_findBestMatch: nullify last unsorted candidate %u",
                    matchIndex);
        *nextCandidate = *unsortedMark = 0;
    }

    /* batch sort stacked candidates */
    matchIndex = previousCandidate;
    while (matchIndex) {  /* will end on matchIndex == 0 */
        U32* const nextCandidateIdxPtr = bt + 2*(matchIndex&btMask) + 1;
        U32 const nextCandidateIdx = *nextCandidateIdxPtr;
        ZSTD_insertDUBT1(ms, matchIndex, iend,
                         nbCandidates, unsortLimit, dictMode);
        matchIndex = nextCandidateIdx;
        nbCandidates++;
    }

    /* find longest match */
    {   size_t commonLengthSmaller = 0, commonLengthLarger = 0;
        const BYTE* const dictBase = ms->window.dictBase;
        const U32 dictLimit = ms->window.dictLimit;
        const BYTE* const dictEnd = dictBase + dictLimit;
        const BYTE* const prefixStart = base + dictLimit;
        U32* smallerPtr = bt + 2*(current&btMask);
        U32* largerPtr  = bt + 2*(current&btMask) + 1;
        U32 matchEndIdx = current + 8 + 1;
        U32 dummy32;   /* to be nullified at the end */
        size_t bestLength = 0;

        matchIndex  = hashTable[h];
        hashTable[h] = current;   /* Update Hash Table */

        while (nbCompares-- && (matchIndex > windowLow)) {
            U32* const nextPtr = bt + 2*(matchIndex & btMask);
            size_t matchLength = MIN(commonLengthSmaller, commonLengthLarger);   /* guaranteed minimum nb of common bytes */
            const BYTE* match;

            if ((dictMode != ZSTD_extDict) || (matchIndex+matchLength >= dictLimit)) {
                match = base + matchIndex;
                matchLength += ZSTD_count(ip+matchLength, match+matchLength, iend);
            } else {
                match = dictBase + matchIndex;
                matchLength += ZSTD_count_2segments(ip+matchLength, match+matchLength, iend, dictEnd, prefixStart);
                if (matchIndex+matchLength >= dictLimit)
                    match = base + matchIndex;   /* to prepare for next usage of match[matchLength] */
            }

            if (matchLength > bestLength) {
                if (matchLength > matchEndIdx - matchIndex)
                    matchEndIdx = matchIndex + (U32)matchLength;
                if ( (4*(int)(matchLength-bestLength)) > (int)(ZSTD_highbit32(current-matchIndex+1) - ZSTD_highbit32((U32)offsetPtr[0]+1)) )
                    bestLength = matchLength, *offsetPtr = ZSTD_REP_MOVE + current - matchIndex;
                if (ip+matchLength == iend) {   /* equal : no way to know if inf or sup */
                    if (dictMode == ZSTD_dictMatchState) {
                        nbCompares = 0; /* in addition to avoiding checking any
                                         * further in this loop, make sure we
                                         * skip checking in the dictionary. */
                    }
                    break;   /* drop, to guarantee consistency (miss a little bit of compression) */
                }
            }

            if (match[matchLength] < ip[matchLength]) {
                /* match is smaller than current */
                *smallerPtr = matchIndex;             /* update smaller idx */
                commonLengthSmaller = matchLength;    /* all smaller will now have at least this guaranteed common length */
                if (matchIndex <= btLow) { smallerPtr=&dummy32; break; }   /* beyond tree size, stop the search */
                smallerPtr = nextPtr+1;               /* new "smaller" => larger of match */
                matchIndex = nextPtr[1];              /* new matchIndex larger than previous (closer to current) */
            } else {
                /* match is larger than current */
                *largerPtr = matchIndex;
                commonLengthLarger = matchLength;
                if (matchIndex <= btLow) { largerPtr=&dummy32; break; }   /* beyond tree size, stop the search */
                largerPtr = nextPtr;
                matchIndex = nextPtr[0];
        }   }

        *smallerPtr = *largerPtr = 0;

        if (dictMode == ZSTD_dictMatchState && nbCompares) {
            bestLength = ZSTD_DUBT_findBetterDictMatch(
                    ms, ip, iend,
                    offsetPtr, bestLength, nbCompares,
                    mls, dictMode);
        }

        assert(matchEndIdx > current+8); /* ensure nextToUpdate is increased */
        ms->nextToUpdate = matchEndIdx - 8;   /* skip repetitive patterns */
        if (bestLength >= MINMATCH) {
            U32 const mIndex = current - ((U32)*offsetPtr - ZSTD_REP_MOVE); (void)mIndex;
            DEBUGLOG(8, "ZSTD_DUBT_findBestMatch(%u) : found match of length %u and offsetCode %u (pos %u)",
                        current, (U32)bestLength, (U32)*offsetPtr, mIndex);
        }
        return bestLength;
    }
}


/** ZSTD_BtFindBestMatch() : Tree updater, providing best match */
FORCE_INLINE_TEMPLATE size_t
ZSTD_BtFindBestMatch( ZSTD_matchState_t* ms,
                const BYTE* const ip, const BYTE* const iLimit,
                      size_t* offsetPtr,
                const U32 mls /* template */,
                const ZSTD_dictMode_e dictMode)
{
    DEBUGLOG(7, "ZSTD_BtFindBestMatch");
    if (ip < ms->window.base + ms->nextToUpdate) return 0;   /* skipped area */
    ZSTD_updateDUBT(ms, ip, iLimit, mls);
    return ZSTD_DUBT_findBestMatch(ms, ip, iLimit, offsetPtr, mls, dictMode);
}


static size_t
ZSTD_BtFindBestMatch_selectMLS (  ZSTD_matchState_t* ms,
                            const BYTE* ip, const BYTE* const iLimit,
                                  size_t* offsetPtr)
{
    switch(ms->cParams.minMatch)
    {
    default : /* includes case 3 */
    case 4 : return ZSTD_BtFindBestMatch(ms, ip, iLimit, offsetPtr, 4, ZSTD_noDict);
    case 5 : return ZSTD_BtFindBestMatch(ms, ip, iLimit, offsetPtr, 5, ZSTD_noDict);
    case 7 :
    case 6 : return ZSTD_BtFindBestMatch(ms, ip, iLimit, offsetPtr, 6, ZSTD_noDict);
    }
}


static size_t ZSTD_BtFindBestMatch_dictMatchState_selectMLS (
                        ZSTD_matchState_t* ms,
                        const BYTE* ip, const BYTE* const iLimit,
                        size_t* offsetPtr)
{
    switch(ms->cParams.minMatch)
    {
    default : /* includes case 3 */
    case 4 : return ZSTD_BtFindBestMatch(ms, ip, iLimit, offsetPtr, 4, ZSTD_dictMatchState);
    case 5 : return ZSTD_BtFindBestMatch(ms, ip, iLimit, offsetPtr, 5, ZSTD_dictMatchState);
    case 7 :
    case 6 : return ZSTD_BtFindBestMatch(ms, ip, iLimit, offsetPtr, 6, ZSTD_dictMatchState);
    }
}


static size_t ZSTD_BtFindBestMatch_extDict_selectMLS (
                        ZSTD_matchState_t* ms,
                        const BYTE* ip, const BYTE* const iLimit,
                        size_t* offsetPtr)
{
    switch(ms->cParams.minMatch)
    {
    default : /* includes case 3 */
    case 4 : return ZSTD_BtFindBestMatch(ms, ip, iLimit, offsetPtr, 4, ZSTD_extDict);
    case 5 : return ZSTD_BtFindBestMatch(ms, ip, iLimit, offsetPtr, 5, ZSTD_extDict);
    case 7 :
    case 6 : return ZSTD_BtFindBestMatch(ms, ip, iLimit, offsetPtr, 6, ZSTD_extDict);
    }
}



/* *********************************
*  Hash Chain
***********************************/
#define NEXT_IN_CHAIN(d, mask)   chainTable[(d) & (mask)]

/* Update chains up to ip (excluded)
   Assumption : always within prefix (i.e. not within extDict) */
static U32 ZSTD_insertAndFindFirstIndex_internal(
                        ZSTD_matchState_t* ms,
                        const ZSTD_compressionParameters* const cParams,
                        const BYTE* ip, U32 const mls)
{
    U32* const hashTable  = ms->hashTable;
    const U32 hashLog = cParams->hashLog;
    U32* const chainTable = ms->chainTable;
    const U32 chainMask = (1 << cParams->chainLog) - 1;
    const BYTE* const base = ms->window.base;
    const U32 target = (U32)(ip - base);
    U32 idx = ms->nextToUpdate;

    while(idx < target) { /* catch up */
        size_t const h = ZSTD_hashPtr(base+idx, hashLog, mls);
        NEXT_IN_CHAIN(idx, chainMask) = hashTable[h];
        hashTable[h] = idx;
        idx++;
    }

    ms->nextToUpdate = target;
    return hashTable[ZSTD_hashPtr(ip, hashLog, mls)];
}

U32 ZSTD_insertAndFindFirstIndex(ZSTD_matchState_t* ms, const BYTE* ip) {
    const ZSTD_compressionParameters* const cParams = &ms->cParams;
    return ZSTD_insertAndFindFirstIndex_internal(ms, cParams, ip, ms->cParams.minMatch);
}


/* inlining is important to hardwire a hot branch (template emulation) */
FORCE_INLINE_TEMPLATE
size_t ZSTD_HcFindBestMatch_generic (
                        ZSTD_matchState_t* ms,
                        const BYTE* const ip, const BYTE* const iLimit,
                        size_t* offsetPtr,
                        const U32 mls, const ZSTD_dictMode_e dictMode)
{
    const ZSTD_compressionParameters* const cParams = &ms->cParams;
    U32* const chainTable = ms->chainTable;
    const U32 chainSize = (1 << cParams->chainLog);
    const U32 chainMask = chainSize-1;
    const BYTE* const base = ms->window.base;
    const BYTE* const dictBase = ms->window.dictBase;
    const U32 dictLimit = ms->window.dictLimit;
    const BYTE* const prefixStart = base + dictLimit;
    const BYTE* const dictEnd = dictBase + dictLimit;
    const U32 lowLimit = ms->window.lowLimit;
    const U32 current = (U32)(ip-base);
    const U32 minChain = current > chainSize ? current - chainSize : 0;
    U32 nbAttempts = 1U << cParams->searchLog;
    size_t ml=4-1;

    /* HC4 match finder */
    U32 matchIndex = ZSTD_insertAndFindFirstIndex_internal(ms, cParams, ip, mls);

    for ( ; (matchIndex>lowLimit) & (nbAttempts>0) ; nbAttempts--) {
        size_t currentMl=0;
        if ((dictMode != ZSTD_extDict) || matchIndex >= dictLimit) {
            const BYTE* const match = base + matchIndex;
            assert(matchIndex >= dictLimit);   /* ensures this is true if dictMode != ZSTD_extDict */
            if (match[ml] == ip[ml])   /* potentially better */
                currentMl = ZSTD_count(ip, match, iLimit);
        } else {
            const BYTE* const match = dictBase + matchIndex;
            assert(match+4 <= dictEnd);
            if (MEM_read32(match) == MEM_read32(ip))   /* assumption : matchIndex <= dictLimit-4 (by table construction) */
                currentMl = ZSTD_count_2segments(ip+4, match+4, iLimit, dictEnd, prefixStart) + 4;
        }

        /* save best solution */
        if (currentMl > ml) {
            ml = currentMl;
            *offsetPtr = current - matchIndex + ZSTD_REP_MOVE;
            if (ip+currentMl == iLimit) break; /* best possible, avoids read overflow on next attempt */
        }

        if (matchIndex <= minChain) break;
        matchIndex = NEXT_IN_CHAIN(matchIndex, chainMask);
    }

    if (dictMode == ZSTD_dictMatchState) {
        const ZSTD_matchState_t* const dms = ms->dictMatchState;
        const U32* const dmsChainTable = dms->chainTable;
        const U32 dmsChainSize         = (1 << dms->cParams.chainLog);
        const U32 dmsChainMask         = dmsChainSize - 1;
        const U32 dmsLowestIndex       = dms->window.dictLimit;
        const BYTE* const dmsBase      = dms->window.base;
        const BYTE* const dmsEnd       = dms->window.nextSrc;
        const U32 dmsSize              = (U32)(dmsEnd - dmsBase);
        const U32 dmsIndexDelta        = dictLimit - dmsSize;
        const U32 dmsMinChain = dmsSize > dmsChainSize ? dmsSize - dmsChainSize : 0;

        matchIndex = dms->hashTable[ZSTD_hashPtr(ip, dms->cParams.hashLog, mls)];

        for ( ; (matchIndex>dmsLowestIndex) & (nbAttempts>0) ; nbAttempts--) {
            size_t currentMl=0;
            const BYTE* const match = dmsBase + matchIndex;
            assert(match+4 <= dmsEnd);
            if (MEM_read32(match) == MEM_read32(ip))   /* assumption : matchIndex <= dictLimit-4 (by table construction) */
                currentMl = ZSTD_count_2segments(ip+4, match+4, iLimit, dmsEnd, prefixStart) + 4;

            /* save best solution */
            if (currentMl > ml) {
                ml = currentMl;
                *offsetPtr = current - (matchIndex + dmsIndexDelta) + ZSTD_REP_MOVE;
                if (ip+currentMl == iLimit) break; /* best possible, avoids read overflow on next attempt */
            }

            if (matchIndex <= dmsMinChain) break;
            matchIndex = dmsChainTable[matchIndex & dmsChainMask];
        }
    }

    return ml;
}

#include <immintrin.h>

static void updateTable(U32* bucket, __m128i const vOld, U32 newValue)
{
    __m128i const vShifted = _mm_slli_si128(vOld, 4);
    __m128i const vNew = _mm_set_epi32(0, 0, 0, newValue);
    __m128i const vCombined = _mm_or_si128(vShifted, vNew);
    _mm_store_si128((__m128i*)bucket, vCombined);
}

#if defined(DEBUGLEVEL) && DEBUGLEVEL >= 2
#define P32_128(x) do { \
    U32 y[4]; \
    _mm_store_si128((__m128i*)y, x); \
    DEBUGLOG(2, "L%u %s = {%u, %u, %u, %u}", __LINE__, #x, y[0], y[1], y[2], y[3]); \
    } while (0)
#define P32_256(x) do { \
    U32 y[8]; \
    _mm256_store_si256((__m256i*)y, x); \
    DEBUGLOG(2, "L%u %s = {%u, %u, %u, %u} {%u, %u, %u, %u}", __LINE__, #x, y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7]); \
    } while (0)
#else
#define P32_128(x)
#define P32_256(x)
#endif

/* inlining is important to hardwire a hot branch (template emulation) */
FORCE_INLINE_TEMPLATE
size_t ZSTD_VectFindBestMatch_generic (
                        ZSTD_matchState_t* ms,
                        const BYTE* const ip, const BYTE* const iLimit,
                        size_t* offsetPtr,
                        const U32 mls)
{
    const ZSTD_compressionParameters* const cParams = &ms->cParams;
    // TODO: Could use the chain-table for a 8-byte hash like double fast
    U32* const hashTableL = ms->hashTable;
    U32* const hashTableS = ms->chainTable;
    U32 const searchLog = 2;
    const U32 hashLogL = cParams->hashLog - 2; /* We have buckets of 4 entries */
    const U32 hashLogS = cParams->chainLog - 2;
    const BYTE* const base = ms->window.base;
    const BYTE* const dictBase = ms->window.dictBase;
    const U32 dictLimit = ms->window.dictLimit;
    const BYTE* const prefixStart = base + dictLimit;
    const BYTE* const dictEnd = dictBase + dictLimit;
    const U32 lowLimit = ms->window.lowLimit;
    
    size_t const hashL = ZSTD_hashPtr(ip, hashLogL, 8) << searchLog; /* Find the bucket */
    size_t const hashS = ZSTD_hashPtr(ip, hashLogS, mls) << searchLog;
    const U32 current = (U32)(ip-base);
    const U32 searches = 1 << searchLog;



    assert(cParams->searchLog == 3); /* We do 8 searches, 4 in each table */
#if 1
    const U32 target = current;
    U32 idx = ms->nextToUpdate;

    while(idx < target) { /* catch up */
        size_t const hL = ZSTD_hashPtr(base+idx, hashLogL, 8) << searchLog;
        size_t const hS = ZSTD_hashPtr(base+idx, hashLogS, mls) << searchLog;
        __m128i const vMatchIndicesL = _mm_loadu_si128((__m128i const*)(hashTableL + hL));
        __m128i const vMatchIndicesS = _mm_loadu_si128((__m128i const*)(hashTableS + hS));
        updateTable(hashTableL + hL, vMatchIndicesL, idx);
        updateTable(hashTableS + hS, vMatchIndicesS, idx);
        idx++;
    }
    ms->nextToUpdate = target;
    /* Load the match candidates. Out of bounds candidates get mapped to
     * values that are 0-byte matches.
     */
    U64      const value          = MEM_read64(ip);
    __m256i const vNotCurrent    = _mm256_set1_epi64x(~value);
    P32_256(vNotCurrent);

    __m128i const vMatchIndicesL = _mm_loadu_si128((__m128i const*)(hashTableL + hashL));
    P32_128(vMatchIndicesL);
    __m128i const vMatchIndicesS = _mm_loadu_si128((__m128i const*)(hashTableS + hashS));
    P32_128(vMatchIndicesS);

    updateTable(hashTableL + hashL, vMatchIndicesL, current);
    updateTable(hashTableS + hashS, vMatchIndicesS, current);

    __m128i const vLowLimit      = _mm_set1_epi32(lowLimit);
    P32_128(vLowLimit);
    __m256i const vValidL        = _mm256_cvtepi32_epi64(_mm_cmpgt_epi32(vMatchIndicesL, vLowLimit));
    __m256i const vValidS        = _mm256_cvtepi32_epi64(_mm_cmpgt_epi32(vMatchIndicesS, vLowLimit));
    __m256i const vMatchesL      = _mm256_mask_i32gather_epi64(vNotCurrent, (long long int const*)base, vMatchIndicesL, vValidL, 1);
    P32_256(vMatchesL);
    __m256i const vMatchesS      = _mm256_mask_i32gather_epi64(vNotCurrent, (long long int const*)base, vMatchIndicesS, vValidS, 1);
    P32_256(vMatchesS);
    /* Compute the match lengths. */
    __m256i const vCurrent       = _mm256_set1_epi64x(value);
    P32_256(vCurrent);
    __m256i const vMatchXorL      = _mm256_xor_si256(vMatchesL, vCurrent);
    __m256i const vMatchXorS      = _mm256_xor_si256(vMatchesS, vCurrent);
    __m256i const vZero           = _mm256_setzero_si256();
    __m256i const vIsByteMatchL   = _mm256_cmpeq_epi8(vMatchXorL, vZero);
    P32_256(vIsByteMatchL);
    __m256i const vIsByteMatchS   = _mm256_cmpeq_epi8(vMatchXorS, vZero);
    P32_256(vIsByteMatchS);
    __m256i const vShiftFinish    = _mm256_set1_epi64x(0x8040201008040201ULL);
    P32_256(vShiftFinish);
    __m256i const vMatchSpreadL   = _mm256_and_si256(vIsByteMatchL, vShiftFinish);
    P32_256(vMatchSpreadL);
    __m256i const vMatchSpreadS   = _mm256_and_si256(vIsByteMatchS, vShiftFinish);
    P32_256(vMatchSpreadS);
    __m256i const vMatchMask256L  = _mm256_sad_epu8(vMatchSpreadL, vZero);
    P32_256(vMatchMask256L);
    __m256i const vMatchMask256S  = _mm256_sad_epu8(vMatchSpreadS, vZero);
    P32_256(vMatchMask256S);
    // __m256i const vDownConvert    = _mm256_set_epi32(0xffffff00, 0xffffff08, -1, -1, -1, -1, 0xffffff00, 0xffffff08);
    __m256i const vDownConvert    = _mm256_set_epi32(0xffffff08, 0xffffff00, -1, -1, -1, -1, 0xffffff08, 0xffffff00);
    P32_256(vDownConvert);
    __m256i const vMatchMask128L  = _mm256_shuffle_epi8(vMatchMask256L, vDownConvert);
    P32_256(vMatchMask128L);
    __m256i const vMatchMask128S  = _mm256_shuffle_epi8(vMatchMask256S, vDownConvert);
    P32_256(vMatchMask128S);

    __m128i const vMatchMaskL0    = _mm256_extracti128_si256(vMatchMask128L, 0);
    __m128i const vMatchMaskL1    = _mm256_extracti128_si256(vMatchMask128L, 1);
    __m128i const vMatchMaskL     = _mm_or_si128(vMatchMaskL0, vMatchMaskL1);
    __m128i const vMatchMaskS0    = _mm256_extracti128_si256(vMatchMask128S, 0);
    __m128i const vMatchMaskS1    = _mm256_extracti128_si256(vMatchMask128S, 1);
    __m128i const vMatchMaskS     = _mm_or_si128(vMatchMaskS0, vMatchMaskS1);
    __m256i const vMatchMask      = _mm256_set_m128i(vMatchMaskL, vMatchMaskS);
    P32_256(vMatchMask);

    __m256i const vMatchMaskHigh  = _mm256_srli_epi32(vMatchMask, 4);
    P32_256(vMatchMaskHigh);
    __m256i const vLowNibbleMask  = _mm256_set1_epi32(0xf);
    __m256i const vMatchMaskLow   = _mm256_and_si256(vMatchMask, vLowNibbleMask);
    P32_256(vMatchMaskLow);

    __m256i const vMatchLowMax    = _mm256_set1_epi32(0x4);
    __m128i const vMatchCountP128 = _mm_set_epi8(
    0x04, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x03, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00);
    // 0x0, 0x1, 0x0, 0x2, 0x0, 0x1, 0x0, 0x3, 0x0, 0x1, 0x0, 0x2, 0x0, 0x1, 0x0, 0x4);
    P32_128(vMatchCountP128);
    __m256i const vMatchCountPerm = _mm256_broadcastsi128_si256(vMatchCountP128);
    P32_256(vMatchCountPerm);

    __m256i const vMatchCountHigh = _mm256_shuffle_epi8(vMatchCountPerm, vMatchMaskHigh);
    P32_256(vMatchCountHigh);
    __m256i const vMatchCountLow  = _mm256_shuffle_epi8(vMatchCountPerm, vMatchMaskLow);
    P32_256(vMatchCountLow);

    __m256i const vMatchedLowMax  = _mm256_cmpeq_epi8(vMatchCountLow, vMatchLowMax);
    P32_256(vMatchedLowMax);
    __m256i const vMatchCountHigh2= _mm256_and_si256(vMatchCountHigh, vMatchedLowMax);
    P32_256(vMatchCountHigh2);

    // _mm256i const vMatchIndices;
    __m256i const vMatchLengths   = _mm256_add_epi32(vMatchCountLow, vMatchCountHigh2);
    P32_256(vMatchLengths);
    U32 matchLengths[8];
    U32 matchIndices[8];
    _mm256_store_si256((__m256i*)matchLengths, vMatchLengths);
    _mm_store_si128((__m128i*)matchIndices, vMatchIndicesS);
    _mm_store_si128((__m128i*)(matchIndices + 4), vMatchIndicesL);

    size_t longest = 4-1;
    for (int i = 0; i < 8; ++i) {
        U32 const matchIndex = matchIndices[i];
        const BYTE* const match = base + matchIndex;
        size_t length = matchLengths[i];
        if (length != 0) {
            assert(matchIndex > lowLimit);
            assert(matchIndex > dictLimit);
            DEBUGLOG(2, "Length %zu vs %zu", length, ZSTD_count(ip, match, iLimit));
            assert(ZSTD_count(ip, match, iLimit) >= length);
        }
        assert(length <= 8);
        if (length == 8) {
            length = ZSTD_count(ip+8, match+8, iLimit)+8;
        }
        if (length > longest) {
            longest = length;
            *offsetPtr = current - matchIndex + ZSTD_REP_MOVE;
        }
    }
    // Update the hash table
    return longest;
#else
    const U32 target = (U32)(ip - base);
    U32 idx = ms->nextToUpdate;

    while(idx < target) { /* catch up */
        size_t const hL = ZSTD_hashPtr(base+idx, hashLogL, 8) << searchLog;
        size_t const hS = ZSTD_hashPtr(base+idx, hashLogS, mls) << searchLog;
        U32* bucketL = hashTableL + hL;
        U32* bucketS = hashTableS + hS;
        for (int i = searches - 1; i > 0; --i) {
            bucketL[i] = bucketL[i-1];
        }
        bucketL[0] = idx;
        for (int i = searches - 1; i > 0; --i) {
            bucketS[i] = bucketS[i-1];
        }
        bucketS[0] = idx;
        idx++;
    }
    ms->nextToUpdate = target;

    size_t ml = 4-1;
    for (int i = 0; i < searches; ++i) {
        U32 const matchIndex = hashTableL[hashL + i];
        const BYTE* const match = base + matchIndex;
        size_t currentMl = 0;
        /* We store matches in order so stop once we're out of bounds */
        if (matchIndex <= lowLimit)
            break;
        assert(matchIndex >= dictLimit);   /* ensures this is true if dictMode != ZSTD_extDict */
        if (match[ml] == ip[ml])   /* potentially better */
            currentMl = ZSTD_count(ip, match, iLimit);
        if (currentMl > ml) {
            ml = currentMl;
            *offsetPtr = current - matchIndex + ZSTD_REP_MOVE;
            if (ip+currentMl == iLimit)
                break; /* best possible, avoids read overflow on next attempt */
        }
    }
    if (ip+ml < iLimit)
        for (int i = 0; i < searches; ++i) {
            U32 const matchIndex = hashTableS[hashS + i];
            const BYTE* const match = base + matchIndex;
            size_t currentMl = 0;
            /* We store matches in order so stop once we're out of bounds */
            if (matchIndex <= lowLimit)
                break;
            assert(matchIndex >= dictLimit);   /* ensures this is true if dictMode != ZSTD_extDict */
            if (match[ml] == ip[ml])   /* potentially better */
                currentMl = ZSTD_count(ip, match, iLimit);
            if (currentMl > ml) {
                ml = currentMl;
                *offsetPtr = current - matchIndex + ZSTD_REP_MOVE;
                if (ip+currentMl == iLimit)
                    break; /* best possible, avoids read overflow on next attempt */
            }
        }
    for (int i = searches - 1; i > 0; --i) {
        hashTableL[hashL + i] = hashTableL[hashL + i-1];
    }
    hashTableL[hashL] = current;
    for (int i = searches - 1; i > 0; --i) {
        hashTableS[hashS + i] = hashTableS[hashS + i-1];
    }
    hashTableS[hashS] = current;
    return ml;
#endif
}

FORCE_INLINE_TEMPLATE size_t ZSTD_VectFindBestMatch_selectMLS (
                        ZSTD_matchState_t* ms,
                        const BYTE* ip, const BYTE* const iLimit,
                        size_t* offsetPtr)
{
    // assert(ip + 1 <= iLimit); /* We want to precompute hashes for ip + 1 for the common case */
    switch(ms->cParams.minMatch)
    {
    default : /* includes case 3 */
    case 4 : return ZSTD_VectFindBestMatch_generic(ms, ip, iLimit, offsetPtr, 4);
    case 5 : return ZSTD_VectFindBestMatch_generic(ms, ip, iLimit, offsetPtr, 5);
    case 7 :
    case 6 : return ZSTD_VectFindBestMatch_generic(ms, ip, iLimit, offsetPtr, 6);
    }
}


FORCE_INLINE_TEMPLATE size_t ZSTD_HcFindBestMatch_selectMLS (
                        ZSTD_matchState_t* ms,
                        const BYTE* ip, const BYTE* const iLimit,
                        size_t* offsetPtr)
{
    switch(ms->cParams.minMatch)
    {
    default : /* includes case 3 */
    case 4 : return ZSTD_HcFindBestMatch_generic(ms, ip, iLimit, offsetPtr, 4, ZSTD_noDict);
    case 5 : return ZSTD_HcFindBestMatch_generic(ms, ip, iLimit, offsetPtr, 5, ZSTD_noDict);
    case 7 :
    case 6 : return ZSTD_HcFindBestMatch_generic(ms, ip, iLimit, offsetPtr, 6, ZSTD_noDict);
    }
}


static size_t ZSTD_HcFindBestMatch_dictMatchState_selectMLS (
                        ZSTD_matchState_t* ms,
                        const BYTE* ip, const BYTE* const iLimit,
                        size_t* offsetPtr)
{
    switch(ms->cParams.minMatch)
    {
    default : /* includes case 3 */
    case 4 : return ZSTD_HcFindBestMatch_generic(ms, ip, iLimit, offsetPtr, 4, ZSTD_dictMatchState);
    case 5 : return ZSTD_HcFindBestMatch_generic(ms, ip, iLimit, offsetPtr, 5, ZSTD_dictMatchState);
    case 7 :
    case 6 : return ZSTD_HcFindBestMatch_generic(ms, ip, iLimit, offsetPtr, 6, ZSTD_dictMatchState);
    }
}


FORCE_INLINE_TEMPLATE size_t ZSTD_HcFindBestMatch_extDict_selectMLS (
                        ZSTD_matchState_t* ms,
                        const BYTE* ip, const BYTE* const iLimit,
                        size_t* offsetPtr)
{
    switch(ms->cParams.minMatch)
    {
    default : /* includes case 3 */
    case 4 : return ZSTD_HcFindBestMatch_generic(ms, ip, iLimit, offsetPtr, 4, ZSTD_extDict);
    case 5 : return ZSTD_HcFindBestMatch_generic(ms, ip, iLimit, offsetPtr, 5, ZSTD_extDict);
    case 7 :
    case 6 : return ZSTD_HcFindBestMatch_generic(ms, ip, iLimit, offsetPtr, 6, ZSTD_extDict);
    }
}


/* *******************************
*  Common parser - lazy strategy
*********************************/
FORCE_INLINE_TEMPLATE
size_t ZSTD_compressBlock_lazy_generic(
                        ZSTD_matchState_t* ms, seqStore_t* seqStore,
                        U32 rep[ZSTD_REP_NUM],
                        const void* src, size_t srcSize,
                        const U32 searchMethod, const U32 depth,
                        ZSTD_dictMode_e const dictMode)
{
    const BYTE* const istart = (const BYTE*)src;
    const BYTE* ip = istart;
    const BYTE* anchor = istart;
    const BYTE* const iend = istart + srcSize;
    const BYTE* const ilimit = iend - 8;
    const BYTE* const base = ms->window.base;
    const U32 prefixLowestIndex = ms->window.dictLimit;
    const BYTE* const prefixLowest = base + prefixLowestIndex;

    typedef size_t (*searchMax_f)(
                        ZSTD_matchState_t* ms,
                        const BYTE* ip, const BYTE* iLimit, size_t* offsetPtr);
    // searchMax_f const searchMax = dictMode == ZSTD_dictMatchState ?
    //     (searchMethod ? ZSTD_BtFindBestMatch_dictMatchState_selectMLS : ZSTD_HcFindBestMatch_dictMatchState_selectMLS) :
    //     (searchMethod ? ZSTD_BtFindBestMatch_selectMLS : ZSTD_HcFindBestMatch_selectMLS);
    searchMax_f const searchMax = ZSTD_VectFindBestMatch_selectMLS;
    U32 offset_1 = rep[0], offset_2 = rep[1], savedOffset=0;

    const ZSTD_matchState_t* const dms = ms->dictMatchState;
    const U32 dictLowestIndex      = dictMode == ZSTD_dictMatchState ?
                                     dms->window.dictLimit : 0;
    const BYTE* const dictBase     = dictMode == ZSTD_dictMatchState ?
                                     dms->window.base : NULL;
    const BYTE* const dictLowest   = dictMode == ZSTD_dictMatchState ?
                                     dictBase + dictLowestIndex : NULL;
    const BYTE* const dictEnd      = dictMode == ZSTD_dictMatchState ?
                                     dms->window.nextSrc : NULL;
    const U32 dictIndexDelta       = dictMode == ZSTD_dictMatchState ?
                                     prefixLowestIndex - (U32)(dictEnd - dictBase) :
                                     0;
    const U32 dictAndPrefixLength = (U32)(ip - prefixLowest + dictEnd - dictLowest);

    /* init */
    ip += (dictAndPrefixLength == 0);
    ms->nextToUpdate3 = ms->nextToUpdate;
    if (dictMode == ZSTD_noDict) {
        U32 const maxRep = (U32)(ip - prefixLowest);
        if (offset_2 > maxRep) savedOffset = offset_2, offset_2 = 0;
        if (offset_1 > maxRep) savedOffset = offset_1, offset_1 = 0;
    }
    if (dictMode == ZSTD_dictMatchState) {
        /* dictMatchState repCode checks don't currently handle repCode == 0
         * disabling. */
        assert(offset_1 <= dictAndPrefixLength);
        assert(offset_2 <= dictAndPrefixLength);
    }

    /* Match Loop */
    while (ip < ilimit) {
        size_t matchLength=0;
        size_t offset=0;
        const BYTE* start=ip+1;

        /* check repCode */
        if (dictMode == ZSTD_dictMatchState) {
            const U32 repIndex = (U32)(ip - base) + 1 - offset_1;
            const BYTE* repMatch = (dictMode == ZSTD_dictMatchState
                                && repIndex < prefixLowestIndex) ?
                                   dictBase + (repIndex - dictIndexDelta) :
                                   base + repIndex;
            if (((U32)((prefixLowestIndex-1) - repIndex) >= 3 /* intentional underflow */)
                && (MEM_read32(repMatch) == MEM_read32(ip+1)) ) {
                const BYTE* repMatchEnd = repIndex < prefixLowestIndex ? dictEnd : iend;
                matchLength = ZSTD_count_2segments(ip+1+4, repMatch+4, iend, repMatchEnd, prefixLowest) + 4;
                if (depth==0) goto _storeSequence;
            }
        }
        if ( dictMode == ZSTD_noDict
          && ((offset_1 > 0) & (MEM_read32(ip+1-offset_1) == MEM_read32(ip+1)))) {
            matchLength = ZSTD_count(ip+1+4, ip+1+4-offset_1, iend) + 4;
            if (depth==0) goto _storeSequence;
        }

        /* first search (depth 0) */
        {   size_t offsetFound = 999999999;
            size_t const ml2 = searchMax(ms, ip, iend, &offsetFound);
            if (ml2 > matchLength)
                matchLength = ml2, start = ip, offset=offsetFound;
        }

        if (matchLength < 4) {
            ip += ((ip-anchor) >> kSearchStrength) + 1;   /* jump faster over incompressible sections */
            continue;
        }

        /* let's try to find a better solution */
        if (depth>=1)
        while (ip<ilimit) {
            ip ++;
            if ( (dictMode == ZSTD_noDict)
              && (offset) && ((offset_1>0) & (MEM_read32(ip) == MEM_read32(ip - offset_1)))) {
                size_t const mlRep = ZSTD_count(ip+4, ip+4-offset_1, iend) + 4;
                int const gain2 = (int)(mlRep * 3);
                int const gain1 = (int)(matchLength*3 - ZSTD_highbit32((U32)offset+1) + 1);
                if ((mlRep >= 4) && (gain2 > gain1))
                    matchLength = mlRep, offset = 0, start = ip;
            }
            if (dictMode == ZSTD_dictMatchState) {
                const U32 repIndex = (U32)(ip - base) - offset_1;
                const BYTE* repMatch = repIndex < prefixLowestIndex ?
                               dictBase + (repIndex - dictIndexDelta) :
                               base + repIndex;
                if (((U32)((prefixLowestIndex-1) - repIndex) >= 3 /* intentional underflow */)
                    && (MEM_read32(repMatch) == MEM_read32(ip)) ) {
                    const BYTE* repMatchEnd = repIndex < prefixLowestIndex ? dictEnd : iend;
                    size_t const mlRep = ZSTD_count_2segments(ip+4, repMatch+4, iend, repMatchEnd, prefixLowest) + 4;
                    int const gain2 = (int)(mlRep * 3);
                    int const gain1 = (int)(matchLength*3 - ZSTD_highbit32((U32)offset+1) + 1);
                    if ((mlRep >= 4) && (gain2 > gain1))
                        matchLength = mlRep, offset = 0, start = ip;
                }
            }
            {   size_t offset2=999999999;
                size_t const ml2 = searchMax(ms, ip, iend, &offset2);
                int const gain2 = (int)(ml2*4 - ZSTD_highbit32((U32)offset2+1));   /* raw approx */
                int const gain1 = (int)(matchLength*4 - ZSTD_highbit32((U32)offset+1) + 4);
                if ((ml2 >= 4) && (gain2 > gain1)) {
                    matchLength = ml2, offset = offset2, start = ip;
                    continue;   /* search a better one */
            }   }

            /* let's find an even better one */
            if ((depth==2) && (ip<ilimit)) {
                ip ++;
                if ( (dictMode == ZSTD_noDict)
                  && (offset) && ((offset_1>0) & (MEM_read32(ip) == MEM_read32(ip - offset_1)))) {
                    size_t const mlRep = ZSTD_count(ip+4, ip+4-offset_1, iend) + 4;
                    int const gain2 = (int)(mlRep * 4);
                    int const gain1 = (int)(matchLength*4 - ZSTD_highbit32((U32)offset+1) + 1);
                    if ((mlRep >= 4) && (gain2 > gain1))
                        matchLength = mlRep, offset = 0, start = ip;
                }
                if (dictMode == ZSTD_dictMatchState) {
                    const U32 repIndex = (U32)(ip - base) - offset_1;
                    const BYTE* repMatch = repIndex < prefixLowestIndex ?
                                   dictBase + (repIndex - dictIndexDelta) :
                                   base + repIndex;
                    if (((U32)((prefixLowestIndex-1) - repIndex) >= 3 /* intentional underflow */)
                        && (MEM_read32(repMatch) == MEM_read32(ip)) ) {
                        const BYTE* repMatchEnd = repIndex < prefixLowestIndex ? dictEnd : iend;
                        size_t const mlRep = ZSTD_count_2segments(ip+4, repMatch+4, iend, repMatchEnd, prefixLowest) + 4;
                        int const gain2 = (int)(mlRep * 4);
                        int const gain1 = (int)(matchLength*4 - ZSTD_highbit32((U32)offset+1) + 1);
                        if ((mlRep >= 4) && (gain2 > gain1))
                            matchLength = mlRep, offset = 0, start = ip;
                    }
                }
                {   size_t offset2=999999999;
                    size_t const ml2 = searchMax(ms, ip, iend, &offset2);
                    int const gain2 = (int)(ml2*4 - ZSTD_highbit32((U32)offset2+1));   /* raw approx */
                    int const gain1 = (int)(matchLength*4 - ZSTD_highbit32((U32)offset+1) + 7);
                    if ((ml2 >= 4) && (gain2 > gain1)) {
                        matchLength = ml2, offset = offset2, start = ip;
                        continue;
            }   }   }
            break;  /* nothing found : store previous solution */
        }

        /* NOTE:
         * start[-offset+ZSTD_REP_MOVE-1] is undefined behavior.
         * (-offset+ZSTD_REP_MOVE-1) is unsigned, and is added to start, which
         * overflows the pointer, which is undefined behavior.
         */
        /* catch up */
        if (offset) {
            if (dictMode == ZSTD_noDict) {
                while ( ((start > anchor) & (start - (offset-ZSTD_REP_MOVE) > prefixLowest))
                     && (start[-1] == (start-(offset-ZSTD_REP_MOVE))[-1]) )  /* only search for offset within prefix */
                    { start--; matchLength++; }
            }
            if (dictMode == ZSTD_dictMatchState) {
                U32 const matchIndex = (U32)((start-base) - (offset - ZSTD_REP_MOVE));
                const BYTE* match = (matchIndex < prefixLowestIndex) ? dictBase + matchIndex - dictIndexDelta : base + matchIndex;
                const BYTE* const mStart = (matchIndex < prefixLowestIndex) ? dictLowest : prefixLowest;
                while ((start>anchor) && (match>mStart) && (start[-1] == match[-1])) { start--; match--; matchLength++; }  /* catch up */
            }
            offset_2 = offset_1; offset_1 = (U32)(offset - ZSTD_REP_MOVE);
        }
        /* store sequence */
_storeSequence:
        {   size_t const litLength = start - anchor;
            ZSTD_storeSeq(seqStore, litLength, anchor, (U32)offset, matchLength-MINMATCH);
            anchor = ip = start + matchLength;
        }

        /* check immediate repcode */
        if (dictMode == ZSTD_dictMatchState) {
            while (ip <= ilimit) {
                U32 const current2 = (U32)(ip-base);
                U32 const repIndex = current2 - offset_2;
                const BYTE* repMatch = dictMode == ZSTD_dictMatchState
                    && repIndex < prefixLowestIndex ?
                        dictBase - dictIndexDelta + repIndex :
                        base + repIndex;
                if ( ((U32)((prefixLowestIndex-1) - (U32)repIndex) >= 3 /* intentional overflow */)
                   && (MEM_read32(repMatch) == MEM_read32(ip)) ) {
                    const BYTE* const repEnd2 = repIndex < prefixLowestIndex ? dictEnd : iend;
                    matchLength = ZSTD_count_2segments(ip+4, repMatch+4, iend, repEnd2, prefixLowest) + 4;
                    offset = offset_2; offset_2 = offset_1; offset_1 = (U32)offset;   /* swap offset_2 <=> offset_1 */
                    ZSTD_storeSeq(seqStore, 0, anchor, 0, matchLength-MINMATCH);
                    ip += matchLength;
                    anchor = ip;
                    continue;
                }
                break;
            }
        }

        if (dictMode == ZSTD_noDict) {
            while ( ((ip <= ilimit) & (offset_2>0))
                 && (MEM_read32(ip) == MEM_read32(ip - offset_2)) ) {
                /* store sequence */
                matchLength = ZSTD_count(ip+4, ip+4-offset_2, iend) + 4;
                offset = offset_2; offset_2 = offset_1; offset_1 = (U32)offset; /* swap repcodes */
                ZSTD_storeSeq(seqStore, 0, anchor, 0, matchLength-MINMATCH);
                ip += matchLength;
                anchor = ip;
                continue;   /* faster when present ... (?) */
    }   }   }

    /* Save reps for next block */
    rep[0] = offset_1 ? offset_1 : savedOffset;
    rep[1] = offset_2 ? offset_2 : savedOffset;

    /* Return the last literals size */
    return iend - anchor;
}


size_t ZSTD_compressBlock_btlazy2(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize)
{
    return ZSTD_compressBlock_lazy_generic(ms, seqStore, rep, src, srcSize, 1, 2, ZSTD_noDict);
}

size_t ZSTD_compressBlock_lazy2(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize)
{
    return ZSTD_compressBlock_lazy_generic(ms, seqStore, rep, src, srcSize, 0, 2, ZSTD_noDict);
}

size_t ZSTD_compressBlock_lazy(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize)
{
    return ZSTD_compressBlock_lazy_generic(ms, seqStore, rep, src, srcSize, 0, 1, ZSTD_noDict);
}

size_t ZSTD_compressBlock_greedy(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize)
{
    return ZSTD_compressBlock_lazy_generic(ms, seqStore, rep, src, srcSize, 0, 0, ZSTD_noDict);
}

size_t ZSTD_compressBlock_btlazy2_dictMatchState(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize)
{
    assert(0);
    return ZSTD_compressBlock_lazy_generic(ms, seqStore, rep, src, srcSize, 1, 2, ZSTD_dictMatchState);
}

size_t ZSTD_compressBlock_lazy2_dictMatchState(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize)
{
    assert(0);
    return ZSTD_compressBlock_lazy_generic(ms, seqStore, rep, src, srcSize, 0, 2, ZSTD_dictMatchState);
}

size_t ZSTD_compressBlock_lazy_dictMatchState(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize)
{
    assert(0);
    return ZSTD_compressBlock_lazy_generic(ms, seqStore, rep, src, srcSize, 0, 1, ZSTD_dictMatchState);
}

size_t ZSTD_compressBlock_greedy_dictMatchState(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize)
{
    assert(0);
    return ZSTD_compressBlock_lazy_generic(ms, seqStore, rep, src, srcSize, 0, 0, ZSTD_dictMatchState);
}


FORCE_INLINE_TEMPLATE
size_t ZSTD_compressBlock_lazy_extDict_generic(
                        ZSTD_matchState_t* ms, seqStore_t* seqStore,
                        U32 rep[ZSTD_REP_NUM],
                        const void* src, size_t srcSize,
                        const U32 searchMethod, const U32 depth)
{
    const BYTE* const istart = (const BYTE*)src;
    const BYTE* ip = istart;
    const BYTE* anchor = istart;
    const BYTE* const iend = istart + srcSize;
    const BYTE* const ilimit = iend - 8;
    const BYTE* const base = ms->window.base;
    const U32 dictLimit = ms->window.dictLimit;
    const U32 lowestIndex = ms->window.lowLimit;
    const BYTE* const prefixStart = base + dictLimit;
    const BYTE* const dictBase = ms->window.dictBase;
    const BYTE* const dictEnd  = dictBase + dictLimit;
    const BYTE* const dictStart  = dictBase + lowestIndex;

    typedef size_t (*searchMax_f)(
                        ZSTD_matchState_t* ms,
                        const BYTE* ip, const BYTE* iLimit, size_t* offsetPtr);
    searchMax_f searchMax = searchMethod ? ZSTD_BtFindBestMatch_extDict_selectMLS : ZSTD_HcFindBestMatch_extDict_selectMLS;

    U32 offset_1 = rep[0], offset_2 = rep[1];

    /* init */
    ms->nextToUpdate3 = ms->nextToUpdate;
    ip += (ip == prefixStart);

    /* Match Loop */
    while (ip < ilimit) {
        size_t matchLength=0;
        size_t offset=0;
        const BYTE* start=ip+1;
        U32 current = (U32)(ip-base);

        /* check repCode */
        {   const U32 repIndex = (U32)(current+1 - offset_1);
            const BYTE* const repBase = repIndex < dictLimit ? dictBase : base;
            const BYTE* const repMatch = repBase + repIndex;
            if (((U32)((dictLimit-1) - repIndex) >= 3) & (repIndex > lowestIndex))   /* intentional overflow */
            if (MEM_read32(ip+1) == MEM_read32(repMatch)) {
                /* repcode detected we should take it */
                const BYTE* const repEnd = repIndex < dictLimit ? dictEnd : iend;
                matchLength = ZSTD_count_2segments(ip+1+4, repMatch+4, iend, repEnd, prefixStart) + 4;
                if (depth==0) goto _storeSequence;
        }   }

        /* first search (depth 0) */
        {   size_t offsetFound = 999999999;
            size_t const ml2 = searchMax(ms, ip, iend, &offsetFound);
            if (ml2 > matchLength)
                matchLength = ml2, start = ip, offset=offsetFound;
        }

         if (matchLength < 4) {
            ip += ((ip-anchor) >> kSearchStrength) + 1;   /* jump faster over incompressible sections */
            continue;
        }

        /* let's try to find a better solution */
        if (depth>=1)
        while (ip<ilimit) {
            ip ++;
            current++;
            /* check repCode */
            if (offset) {
                const U32 repIndex = (U32)(current - offset_1);
                const BYTE* const repBase = repIndex < dictLimit ? dictBase : base;
                const BYTE* const repMatch = repBase + repIndex;
                if (((U32)((dictLimit-1) - repIndex) >= 3) & (repIndex > lowestIndex))  /* intentional overflow */
                if (MEM_read32(ip) == MEM_read32(repMatch)) {
                    /* repcode detected */
                    const BYTE* const repEnd = repIndex < dictLimit ? dictEnd : iend;
                    size_t const repLength = ZSTD_count_2segments(ip+4, repMatch+4, iend, repEnd, prefixStart) + 4;
                    int const gain2 = (int)(repLength * 3);
                    int const gain1 = (int)(matchLength*3 - ZSTD_highbit32((U32)offset+1) + 1);
                    if ((repLength >= 4) && (gain2 > gain1))
                        matchLength = repLength, offset = 0, start = ip;
            }   }

            /* search match, depth 1 */
            {   size_t offset2=999999999;
                size_t const ml2 = searchMax(ms, ip, iend, &offset2);
                int const gain2 = (int)(ml2*4 - ZSTD_highbit32((U32)offset2+1));   /* raw approx */
                int const gain1 = (int)(matchLength*4 - ZSTD_highbit32((U32)offset+1) + 4);
                if ((ml2 >= 4) && (gain2 > gain1)) {
                    matchLength = ml2, offset = offset2, start = ip;
                    continue;   /* search a better one */
            }   }

            /* let's find an even better one */
            if ((depth==2) && (ip<ilimit)) {
                ip ++;
                current++;
                /* check repCode */
                if (offset) {
                    const U32 repIndex = (U32)(current - offset_1);
                    const BYTE* const repBase = repIndex < dictLimit ? dictBase : base;
                    const BYTE* const repMatch = repBase + repIndex;
                    if (((U32)((dictLimit-1) - repIndex) >= 3) & (repIndex > lowestIndex))  /* intentional overflow */
                    if (MEM_read32(ip) == MEM_read32(repMatch)) {
                        /* repcode detected */
                        const BYTE* const repEnd = repIndex < dictLimit ? dictEnd : iend;
                        size_t const repLength = ZSTD_count_2segments(ip+4, repMatch+4, iend, repEnd, prefixStart) + 4;
                        int const gain2 = (int)(repLength * 4);
                        int const gain1 = (int)(matchLength*4 - ZSTD_highbit32((U32)offset+1) + 1);
                        if ((repLength >= 4) && (gain2 > gain1))
                            matchLength = repLength, offset = 0, start = ip;
                }   }

                /* search match, depth 2 */
                {   size_t offset2=999999999;
                    size_t const ml2 = searchMax(ms, ip, iend, &offset2);
                    int const gain2 = (int)(ml2*4 - ZSTD_highbit32((U32)offset2+1));   /* raw approx */
                    int const gain1 = (int)(matchLength*4 - ZSTD_highbit32((U32)offset+1) + 7);
                    if ((ml2 >= 4) && (gain2 > gain1)) {
                        matchLength = ml2, offset = offset2, start = ip;
                        continue;
            }   }   }
            break;  /* nothing found : store previous solution */
        }

        /* catch up */
        if (offset) {
            U32 const matchIndex = (U32)((start-base) - (offset - ZSTD_REP_MOVE));
            const BYTE* match = (matchIndex < dictLimit) ? dictBase + matchIndex : base + matchIndex;
            const BYTE* const mStart = (matchIndex < dictLimit) ? dictStart : prefixStart;
            while ((start>anchor) && (match>mStart) && (start[-1] == match[-1])) { start--; match--; matchLength++; }  /* catch up */
            offset_2 = offset_1; offset_1 = (U32)(offset - ZSTD_REP_MOVE);
        }

        /* store sequence */
_storeSequence:
        {   size_t const litLength = start - anchor;
            ZSTD_storeSeq(seqStore, litLength, anchor, (U32)offset, matchLength-MINMATCH);
            anchor = ip = start + matchLength;
        }

        /* check immediate repcode */
        while (ip <= ilimit) {
            const U32 repIndex = (U32)((ip-base) - offset_2);
            const BYTE* const repBase = repIndex < dictLimit ? dictBase : base;
            const BYTE* const repMatch = repBase + repIndex;
            if (((U32)((dictLimit-1) - repIndex) >= 3) & (repIndex > lowestIndex))  /* intentional overflow */
            if (MEM_read32(ip) == MEM_read32(repMatch)) {
                /* repcode detected we should take it */
                const BYTE* const repEnd = repIndex < dictLimit ? dictEnd : iend;
                matchLength = ZSTD_count_2segments(ip+4, repMatch+4, iend, repEnd, prefixStart) + 4;
                offset = offset_2; offset_2 = offset_1; offset_1 = (U32)offset;   /* swap offset history */
                ZSTD_storeSeq(seqStore, 0, anchor, 0, matchLength-MINMATCH);
                ip += matchLength;
                anchor = ip;
                continue;   /* faster when present ... (?) */
            }
            break;
    }   }

    /* Save reps for next block */
    rep[0] = offset_1;
    rep[1] = offset_2;

    /* Return the last literals size */
    return iend - anchor;
}


size_t ZSTD_compressBlock_greedy_extDict(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize)
{
    assert(0);
    return ZSTD_compressBlock_lazy_extDict_generic(ms, seqStore, rep, src, srcSize, 0, 0);
}

size_t ZSTD_compressBlock_lazy_extDict(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize)

{
    assert(0);
    return ZSTD_compressBlock_lazy_extDict_generic(ms, seqStore, rep, src, srcSize, 0, 1);
}

size_t ZSTD_compressBlock_lazy2_extDict(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize)

{
    assert(0);
    return ZSTD_compressBlock_lazy_extDict_generic(ms, seqStore, rep, src, srcSize, 0, 2);
}

size_t ZSTD_compressBlock_btlazy2_extDict(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize)

{
    assert(0);
    return ZSTD_compressBlock_lazy_extDict_generic(ms, seqStore, rep, src, srcSize, 1, 2);
}
