/* ******************************************************************
   Huffman decoder, part of New Generation Entropy library
   Copyright (C) 2013-2016, Yann Collet.

   BSD 2-Clause License (http://www.opensource.org/licenses/bsd-license.php)

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are
   met:

       * Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
       * Redistributions in binary form must reproduce the above
   copyright notice, this list of conditions and the following disclaimer
   in the documentation and/or other materials provided with the
   distribution.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
   OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

    You can contact the author at :
    - FSE+HUF source repository : https://github.com/Cyan4973/FiniteStateEntropy
    - Public forum : https://groups.google.com/forum/#!forum/lz4c
****************************************************************** */

/* **************************************************************
*  Compiler specifics
****************************************************************/
#if defined (__cplusplus) || (defined (__STDC_VERSION__) && (__STDC_VERSION__ >= 199901L) /* C99 */)
/* inline is defined */
#elif defined(_MSC_VER) || defined(__GNUC__)
#  define inline __inline
#else
#  define inline /* disable inline */
#endif

#ifdef _MSC_VER    /* Visual Studio */
#  pragma warning(disable : 4127)        /* disable: C4127: conditional expression is constant */
#endif


/* **************************************************************
*  Dependencies
****************************************************************/
#include <string.h>     /* memcpy, memset */
#include "bitstream.h"  /* BIT_* */
#include "fse.h"        /* header compression */
#define HUF_STATIC_LINKING_ONLY
#include "huf.h"


/* **************************************************************
*  Error Management
****************************************************************/
#define HUF_STATIC_ASSERT(c) { enum { HUF_static_assert = 1/(int)(!!(c)) }; }   /* use only *after* variable declarations */


/*-***************************/
/*  generic DTableDesc       */
/*-***************************/

typedef struct { BYTE maxTableLog; BYTE tableType; BYTE tableLog; BYTE reserved; } DTableDesc;

static DTableDesc HUF_getDTableDesc(const HUF_DTable* table)
{
    DTableDesc dtd;
    memcpy(&dtd, table, sizeof(dtd));
    return dtd;
}


/*-***************************/
/*  single-symbol decoding   */
/*-***************************/


typedef struct { BYTE byte; BYTE nbBits; } HUF_DEltX2;   /* single-symbol decoding */


size_t HUF_readDTableX2 (HUF_DTable* DTable, const void* src, size_t srcSize)
{
    BYTE huffWeight[HUF_SYMBOLVALUE_MAX + 1];
    U32 rankVal[HUF_TABLELOG_ABSOLUTEMAX + 1];   /* large enough for values from 0 to 16 */
    U32 tableLog = 0;
    U32 nbSymbols = 0;
    size_t iSize;
    void* const dtPtr = DTable + 1;
    HUF_DEltX2* const dt = (HUF_DEltX2*)dtPtr;

    HUF_STATIC_ASSERT(sizeof(DTableDesc) == sizeof(HUF_DTable));
    /* memset(huffWeight, 0, sizeof(huffWeight)); */   /* is not necessary, even though some analyzer complain ... */

    iSize = HUF_readStats(huffWeight, HUF_SYMBOLVALUE_MAX + 1, rankVal, &nbSymbols, &tableLog, src, srcSize);
    if (HUF_isError(iSize)) return iSize;

    /* Table header */
    {   DTableDesc dtd = HUF_getDTableDesc(DTable);
        if (tableLog > (U32)(dtd.maxTableLog+1)) return ERROR(tableLog_tooLarge);   /* DTable too small, huffman tree cannot fit in */
        dtd.tableType = 0;
        dtd.tableLog = (BYTE)tableLog;
        memcpy(DTable, &dtd, sizeof(dtd));
    }

    /* Prepare ranks */
    {   U32 n, nextRankStart = 0;
        for (n=1; n<tableLog+1; n++) {
            U32 current = nextRankStart;
            nextRankStart += (rankVal[n] << (n-1));
            rankVal[n] = current;
    }   }

    /* fill DTable */
    {   U32 n;
        for (n=0; n<nbSymbols; n++) {
            U32 const w = huffWeight[n];
            U32 const length = (1 << w) >> 1;
            U32 i;
            HUF_DEltX2 D;
            D.byte = (BYTE)n; D.nbBits = (BYTE)(tableLog + 1 - w);
            for (i = rankVal[w]; i < rankVal[w] + length; i++)
                dt[i] = D;
            rankVal[w] += length;
    }   }

    return iSize;
}

/* *************************/
/* double-symbols decoding */
/* *************************/


typedef struct { U16 sequence; BYTE nbBits; BYTE length; } HUF_DEltX4;  /* double-symbols decoding */

typedef struct { BYTE symbol; BYTE weight; } sortedSymbol_t;


static void HUF_fillDTableX4Level2(HUF_DEltX4* DTable, U32 sizeLog, const U32 consumed,
                           const U32* rankValOrigin, const int minWeight,
                           const sortedSymbol_t* sortedSymbols, const U32 sortedListSize,
                           U32 nbBitsBaseline, U16 baseSeq)
{
    HUF_DEltX4 DElt;
    U32 rankVal[HUF_TABLELOG_ABSOLUTEMAX + 1];

    /* get pre-calculated rankVal */
    memcpy(rankVal, rankValOrigin, sizeof(rankVal));

    /* fill skipped values */
    if (minWeight>1) {
        U32 i, skipSize = rankVal[minWeight];
        MEM_writeLE16(&(DElt.sequence), baseSeq);
        DElt.nbBits   = (BYTE)(consumed);
        DElt.length   = 1;
        for (i = 0; i < skipSize; i++)
            DTable[i] = DElt;
    }

    /* fill DTable */
    {   U32 s; for (s=0; s<sortedListSize; s++) {   /* note : sortedSymbols already skipped */
            const U32 symbol = sortedSymbols[s].symbol;
            const U32 weight = sortedSymbols[s].weight;
            const U32 nbBits = nbBitsBaseline - weight;
            const U32 length = 1 << (sizeLog-nbBits);
            const U32 start = rankVal[weight];
            U32 i = start;
            const U32 end = start + length;

            MEM_writeLE16(&(DElt.sequence), (U16)(baseSeq + (symbol << 8)));
            DElt.nbBits = (BYTE)(nbBits + consumed);
            DElt.length = 2;
            do { DTable[i++] = DElt; } while (i<end);   /* since length >= 1 */

            rankVal[weight] += length;
    }   }
}

typedef U32 rankVal_t[HUF_TABLELOG_ABSOLUTEMAX][HUF_TABLELOG_ABSOLUTEMAX + 1];

static void HUF_fillDTableX4(HUF_DEltX4* DTable, const U32 targetLog,
                           const sortedSymbol_t* sortedList, const U32 sortedListSize,
                           const U32* rankStart, rankVal_t rankValOrigin, const U32 maxWeight,
                           const U32 nbBitsBaseline)
{
    U32 rankVal[HUF_TABLELOG_ABSOLUTEMAX + 1];
    const int scaleLog = nbBitsBaseline - targetLog;   /* note : targetLog >= srcLog, hence scaleLog <= 1 */
    const U32 minBits  = nbBitsBaseline - maxWeight;
    U32 s;

    memcpy(rankVal, rankValOrigin, sizeof(rankVal));

    /* fill DTable */
    for (s=0; s<sortedListSize; s++) {
        const U16 symbol = sortedList[s].symbol;
        const U32 weight = sortedList[s].weight;
        const U32 nbBits = nbBitsBaseline - weight;
        const U32 start = rankVal[weight];
        const U32 length = 1 << (targetLog-nbBits);

        if (targetLog-nbBits >= minBits) {   /* enough room for a second symbol */
            U32 sortedRank;
            int minWeight = nbBits + scaleLog;
            if (minWeight < 1) minWeight = 1;
            sortedRank = rankStart[minWeight];
            HUF_fillDTableX4Level2(DTable+start, targetLog-nbBits, nbBits,
                           rankValOrigin[nbBits], minWeight,
                           sortedList+sortedRank, sortedListSize-sortedRank,
                           nbBitsBaseline, symbol);
        } else {
            HUF_DEltX4 DElt;
            MEM_writeLE16(&(DElt.sequence), symbol);
            DElt.nbBits = (BYTE)(nbBits);
            DElt.length = 1;
            {   U32 const end = start + length;
                U32 u;
                for (u = start; u < end; u++) DTable[u] = DElt;
        }   }
        rankVal[weight] += length;
    }
}

size_t HUF_readDTableX4 (HUF_DTable* DTable, const void* src, size_t srcSize)
{
    BYTE weightList[HUF_SYMBOLVALUE_MAX + 1];
    sortedSymbol_t sortedSymbol[HUF_SYMBOLVALUE_MAX + 1];
    U32 rankStats[HUF_TABLELOG_ABSOLUTEMAX + 1] = { 0 };
    U32 rankStart0[HUF_TABLELOG_ABSOLUTEMAX + 2] = { 0 };
    U32* const rankStart = rankStart0+1;
    rankVal_t rankVal;
    U32 tableLog, maxW, sizeOfSort, nbSymbols;
    DTableDesc dtd = HUF_getDTableDesc(DTable);
    U32 const maxTableLog = dtd.maxTableLog;
    size_t iSize;
    void* dtPtr = DTable+1;   /* force compiler to avoid strict-aliasing */
    HUF_DEltX4* const dt = (HUF_DEltX4*)dtPtr;

    HUF_STATIC_ASSERT(sizeof(HUF_DEltX4) == sizeof(HUF_DTable));   /* if compilation fails here, assertion is false */
    if (maxTableLog > HUF_TABLELOG_ABSOLUTEMAX) return ERROR(tableLog_tooLarge);
    /* memset(weightList, 0, sizeof(weightList)); */  /* is not necessary, even though some analyzer complain ... */

    iSize = HUF_readStats(weightList, HUF_SYMBOLVALUE_MAX + 1, rankStats, &nbSymbols, &tableLog, src, srcSize);
    if (HUF_isError(iSize)) return iSize;

    /* check result */
    if (tableLog > maxTableLog) return ERROR(tableLog_tooLarge);   /* DTable can't fit code depth */

    /* find maxWeight */
    for (maxW = tableLog; rankStats[maxW]==0; maxW--) {}  /* necessarily finds a solution before 0 */

    /* Get start index of each weight */
    {   U32 w, nextRankStart = 0;
        for (w=1; w<maxW+1; w++) {
            U32 current = nextRankStart;
            nextRankStart += rankStats[w];
            rankStart[w] = current;
        }
        rankStart[0] = nextRankStart;   /* put all 0w symbols at the end of sorted list*/
        sizeOfSort = nextRankStart;
    }

    /* sort symbols by weight */
    {   U32 s;
        for (s=0; s<nbSymbols; s++) {
            U32 const w = weightList[s];
            U32 const r = rankStart[w]++;
            sortedSymbol[r].symbol = (BYTE)s;
            sortedSymbol[r].weight = (BYTE)w;
        }
        rankStart[0] = 0;   /* forget 0w symbols; this is beginning of weight(1) */
    }

    /* Build rankVal */
    {   U32* const rankVal0 = rankVal[0];
        {   int const rescale = (maxTableLog-tableLog) - 1;   /* tableLog <= maxTableLog */
            U32 nextRankVal = 0;
            U32 w;
            for (w=1; w<maxW+1; w++) {
                U32 current = nextRankVal;
                nextRankVal += rankStats[w] << (w+rescale);
                rankVal0[w] = current;
        }   }
        {   U32 const minBits = tableLog+1 - maxW;
            U32 consumed;
            for (consumed = minBits; consumed < maxTableLog - minBits + 1; consumed++) {
                U32* const rankValPtr = rankVal[consumed];
                U32 w;
                for (w = 1; w < maxW+1; w++) {
                    rankValPtr[w] = rankVal0[w] >> consumed;
    }   }   }   }

    HUF_fillDTableX4(dt, maxTableLog,
                   sortedSymbol, sizeOfSort,
                   rankStart0, rankVal, maxW,
                   tableLog+1);

    dtd.tableLog = (BYTE)maxTableLog;
    dtd.tableType = 1;
    memcpy(DTable, &dtd, sizeof(dtd));
    return iSize;
}


/* ********************************/
/* Generic decompression selector */
/* ********************************/

size_t HUF_decompress1X_usingDTable(void* dst, size_t maxDstSize,
                                    const void* cSrc, size_t cSrcSize,
                                    const HUF_DTable* DTable)
{
    DTableDesc const dtd = HUF_getDTableDesc(DTable);
    return dtd.tableType ? HUF_decompress1X4_usingDTable_internal(dst, maxDstSize, cSrc, cSrcSize, DTable) :
                           HUF_decompress1X2_usingDTable_internal(dst, maxDstSize, cSrc, cSrcSize, DTable);
}

size_t HUF_decompress4X_usingDTable(void* dst, size_t maxDstSize,
                                    const void* cSrc, size_t cSrcSize,
                                    const HUF_DTable* DTable)
{
    DTableDesc const dtd = HUF_getDTableDesc(DTable);
    return dtd.tableType ? HUF_decompress4X4_usingDTable_internal(dst, maxDstSize, cSrc, cSrcSize, DTable) :
                           HUF_decompress4X2_usingDTable_internal(dst, maxDstSize, cSrc, cSrcSize, DTable);
}


typedef struct { U32 tableTime; U32 decode256Time; } algo_time_t;
static const algo_time_t algoTime[16 /* Quantization */][3 /* single, double, quad */] =
{
    /* single, double, quad */
    {{0,0}, {1,1}, {2,2}},  /* Q==0 : impossible */
    {{0,0}, {1,1}, {2,2}},  /* Q==1 : impossible */
    {{  38,130}, {1313, 74}, {2151, 38}},   /* Q == 2 : 12-18% */
    {{ 448,128}, {1353, 74}, {2238, 41}},   /* Q == 3 : 18-25% */
    {{ 556,128}, {1353, 74}, {2238, 47}},   /* Q == 4 : 25-32% */
    {{ 714,128}, {1418, 74}, {2436, 53}},   /* Q == 5 : 32-38% */
    {{ 883,128}, {1437, 74}, {2464, 61}},   /* Q == 6 : 38-44% */
    {{ 897,128}, {1515, 75}, {2622, 68}},   /* Q == 7 : 44-50% */
    {{ 926,128}, {1613, 75}, {2730, 75}},   /* Q == 8 : 50-56% */
    {{ 947,128}, {1729, 77}, {3359, 77}},   /* Q == 9 : 56-62% */
    {{1107,128}, {2083, 81}, {4006, 84}},   /* Q ==10 : 62-69% */
    {{1177,128}, {2379, 87}, {4785, 88}},   /* Q ==11 : 69-75% */
    {{1242,128}, {2415, 93}, {5155, 84}},   /* Q ==12 : 75-81% */
    {{1349,128}, {2644,106}, {5260,106}},   /* Q ==13 : 81-87% */
    {{1455,128}, {2422,124}, {4174,124}},   /* Q ==14 : 87-93% */
    {{ 722,128}, {1891,145}, {1936,146}},   /* Q ==15 : 93-99% */
};

/** HUF_selectDecoder() :
*   Tells which decoder is likely to decode faster,
*   based on a set of pre-determined metrics.
*   @return : 0==HUF_decompress4X2, 1==HUF_decompress4X4 .
*   Assumption : 0 < cSrcSize < dstSize <= 128 KB */
U32 HUF_selectDecoder (size_t dstSize, size_t cSrcSize)
{
    /* decoder timing evaluation */
    U32 const Q = (U32)(cSrcSize * 16 / dstSize);   /* Q < 16 since dstSize > cSrcSize */
    U32 const D256 = (U32)(dstSize >> 8);
    U32 const DTime0 = algoTime[Q][0].tableTime + (algoTime[Q][0].decode256Time * D256);
    U32 DTime1 = algoTime[Q][1].tableTime + (algoTime[Q][1].decode256Time * D256);
    DTime1 += DTime1 >> 3;  /* advantage to algorithm using less memory, for cache eviction */

    return DTime1 < DTime0;
}


typedef size_t (*decompressionAlgo)(void* dst, size_t dstSize, const void* cSrc, size_t cSrcSize);

size_t HUF_decompress (void* dst, size_t dstSize, const void* cSrc, size_t cSrcSize)
{
    static const decompressionAlgo decompress[2] = { HUF_decompress4X2, HUF_decompress4X4 };

    /* validation checks */
    if (dstSize == 0) return ERROR(dstSize_tooSmall);
    if (cSrcSize > dstSize) return ERROR(corruption_detected);   /* invalid */
    if (cSrcSize == dstSize) { memcpy(dst, cSrc, dstSize); return dstSize; }   /* not compressed */
    if (cSrcSize == 1) { memset(dst, *(const BYTE*)cSrc, dstSize); return dstSize; }   /* RLE */

    {   U32 const algoNb = HUF_selectDecoder(dstSize, cSrcSize);
        return decompress[algoNb](dst, dstSize, cSrc, cSrcSize);
    }
}

size_t HUF_decompress4X_DCtx (HUF_DTable* dctx, void* dst, size_t dstSize, const void* cSrc, size_t cSrcSize)
{
    /* validation checks */
    if (dstSize == 0) return ERROR(dstSize_tooSmall);
    if (cSrcSize > dstSize) return ERROR(corruption_detected);   /* invalid */
    if (cSrcSize == dstSize) { memcpy(dst, cSrc, dstSize); return dstSize; }   /* not compressed */
    if (cSrcSize == 1) { memset(dst, *(const BYTE*)cSrc, dstSize); return dstSize; }   /* RLE */

    {   U32 const algoNb = HUF_selectDecoder(dstSize, cSrcSize);
        return algoNb ? HUF_decompress4X4_DCtx(dctx, dst, dstSize, cSrc, cSrcSize) :
                        HUF_decompress4X2_DCtx(dctx, dst, dstSize, cSrc, cSrcSize) ;
    }
}

size_t HUF_decompress4X_hufOnly (HUF_DTable* dctx, void* dst, size_t dstSize, const void* cSrc, size_t cSrcSize)
{
    /* validation checks */
    if (dstSize == 0) return ERROR(dstSize_tooSmall);
    if ((cSrcSize >= dstSize) || (cSrcSize <= 1)) return ERROR(corruption_detected);   /* invalid */

    {   U32 const algoNb = HUF_selectDecoder(dstSize, cSrcSize);
        return algoNb ? HUF_decompress4X4_DCtx(dctx, dst, dstSize, cSrc, cSrcSize) :
                        HUF_decompress4X2_DCtx(dctx, dst, dstSize, cSrc, cSrcSize) ;
    }
}

size_t HUF_decompress1X_DCtx (HUF_DTable* dctx, void* dst, size_t dstSize, const void* cSrc, size_t cSrcSize)
{
    /* validation checks */
    if (dstSize == 0) return ERROR(dstSize_tooSmall);
    if (cSrcSize > dstSize) return ERROR(corruption_detected);   /* invalid */
    if (cSrcSize == dstSize) { memcpy(dst, cSrc, dstSize); return dstSize; }   /* not compressed */
    if (cSrcSize == 1) { memset(dst, *(const BYTE*)cSrc, dstSize); return dstSize; }   /* RLE */

    {   U32 const algoNb = HUF_selectDecoder(dstSize, cSrcSize);
        return algoNb ? HUF_decompress1X4_DCtx(dctx, dst, dstSize, cSrc, cSrcSize) :
                        HUF_decompress1X2_DCtx(dctx, dst, dstSize, cSrc, cSrcSize) ;
    }
}
