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


typedef struct { U16 sequence; BYTE nbBits; BYTE length; } HUF_DEltX4;  /* double-symbols decoding */

typedef struct { BYTE symbol; BYTE weight; } sortedSymbol_t;


static U32 HUF_decodeSymbolX4(void* op, BIT_DStream_t* DStream, const HUF_DEltX4* dt, const U32 dtLog)
{
    size_t const val = BIT_lookBitsFast(DStream, dtLog);   /* note : dtLog >= 1 */
    memcpy(op, dt+val, 2);
    BIT_skipBits(DStream, dt[val].nbBits);
    return dt[val].length;
}

static U32 HUF_decodeLastSymbolX4(void* op, BIT_DStream_t* DStream, const HUF_DEltX4* dt, const U32 dtLog)
{
    size_t const val = BIT_lookBitsFast(DStream, dtLog);   /* note : dtLog >= 1 */
    memcpy(op, dt+val, 1);
    if (dt[val].length==1) BIT_skipBits(DStream, dt[val].nbBits);
    else {
        if (DStream->bitsConsumed < (sizeof(DStream->bitContainer)*8)) {
            BIT_skipBits(DStream, dt[val].nbBits);
            if (DStream->bitsConsumed > (sizeof(DStream->bitContainer)*8))
                DStream->bitsConsumed = (sizeof(DStream->bitContainer)*8);   /* ugly hack; works only because it's the last symbol. Note : can't easily extract nbBits from just this symbol */
    }   }
    return 1;
}


#define HUF_DECODE_SYMBOLX4_0(ptr, DStreamPtr) \
    ptr += HUF_decodeSymbolX4(ptr, DStreamPtr, dt, dtLog)

#define HUF_DECODE_SYMBOLX4_1(ptr, DStreamPtr) \
    if (MEM_64bits() || (HUF_TABLELOG_MAX<=12)) \
        ptr += HUF_decodeSymbolX4(ptr, DStreamPtr, dt, dtLog)

#define HUF_DECODE_SYMBOLX4_2(ptr, DStreamPtr) \
    if (MEM_64bits()) \
        ptr += HUF_decodeSymbolX4(ptr, DStreamPtr, dt, dtLog)

static inline size_t HUF_decodeStreamX4(BYTE* p, BIT_DStream_t* bitDPtr, BYTE* const pEnd, const HUF_DEltX4* const dt, const U32 dtLog)
{
    BYTE* const pStart = p;

    /* up to 8 symbols at a time */
    while ((BIT_reloadDStream(bitDPtr) == BIT_DStream_unfinished) & (p < pEnd-(sizeof(bitDPtr->bitContainer)-1))) {
        HUF_DECODE_SYMBOLX4_2(p, bitDPtr);
        HUF_DECODE_SYMBOLX4_1(p, bitDPtr);
        HUF_DECODE_SYMBOLX4_2(p, bitDPtr);
        HUF_DECODE_SYMBOLX4_0(p, bitDPtr);
    }

    /* closer to end : up to 2 symbols at a time */
    while ((BIT_reloadDStream(bitDPtr) == BIT_DStream_unfinished) & (p <= pEnd-2))
        HUF_DECODE_SYMBOLX4_0(p, bitDPtr);

    while (p <= pEnd-2)
        HUF_DECODE_SYMBOLX4_0(p, bitDPtr);   /* no need to reload : reached the end of DStream */

    if (p < pEnd)
        p += HUF_decodeLastSymbolX4(p, bitDPtr, dt, dtLog);

    return p-pStart;
}


size_t HUF_decompress1X4_usingDTable_internal(
          void* dst,  size_t dstSize,
    const void* cSrc, size_t cSrcSize,
    const HUF_DTable* DTable)
{
    BIT_DStream_t bitD;

    /* Init */
    {   size_t const errorCode = BIT_initDStream(&bitD, cSrc, cSrcSize);
        if (HUF_isError(errorCode)) return errorCode;
    }

    /* decode */
    {   BYTE* const ostart = (BYTE*) dst;
        BYTE* const oend = ostart + dstSize;
        const void* const dtPtr = DTable+1;   /* force compiler to not use strict-aliasing */
        const HUF_DEltX4* const dt = (const HUF_DEltX4*)dtPtr;
        DTableDesc const dtd = HUF_getDTableDesc(DTable);
        HUF_decodeStreamX4(ostart, &bitD, oend, dt, dtd.tableLog);
    }

    /* check */
    if (!BIT_endOfDStream(&bitD)) return ERROR(corruption_detected);

    /* decoded size */
    return dstSize;
}

size_t HUF_decompress1X4_usingDTable(
          void* dst,  size_t dstSize,
    const void* cSrc, size_t cSrcSize,
    const HUF_DTable* DTable)
{
    DTableDesc dtd = HUF_getDTableDesc(DTable);
    if (dtd.tableType != 1) return ERROR(GENERIC);
    return HUF_decompress1X4_usingDTable_internal(dst, dstSize, cSrc, cSrcSize, DTable);
}

size_t HUF_decompress1X4_DCtx (HUF_DTable* DCtx, void* dst, size_t dstSize, const void* cSrc, size_t cSrcSize)
{
    const BYTE* ip = (const BYTE*) cSrc;

    size_t const hSize = HUF_readDTableX4 (DCtx, cSrc, cSrcSize);
    if (HUF_isError(hSize)) return hSize;
    if (hSize >= cSrcSize) return ERROR(srcSize_wrong);
    ip += hSize; cSrcSize -= hSize;

    return HUF_decompress1X4_usingDTable_internal (dst, dstSize, ip, cSrcSize, DCtx);
}

size_t HUF_decompress1X4 (void* dst, size_t dstSize, const void* cSrc, size_t cSrcSize)
{
    HUF_CREATE_STATIC_DTABLEX4(DTable, HUF_TABLELOG_MAX);
    return HUF_decompress1X4_DCtx(DTable, dst, dstSize, cSrc, cSrcSize);
}

size_t HUF_decompress4X4_usingDTable_internal(
          void* dst,  size_t dstSize,
    const void* cSrc, size_t cSrcSize,
    const HUF_DTable* DTable)
{
    if (cSrcSize < 10) return ERROR(corruption_detected);   /* strict minimum : jump table + 1 byte per stream */

    {   const BYTE* const istart = (const BYTE*) cSrc;
        BYTE* const ostart = (BYTE*) dst;
        BYTE* const oend = ostart + dstSize;
        const void* const dtPtr = DTable+1;
        const HUF_DEltX4* const dt = (const HUF_DEltX4*)dtPtr;

        /* Init */
        BIT_DStream_t bitD1;
        BIT_DStream_t bitD2;
        BIT_DStream_t bitD3;
        BIT_DStream_t bitD4;
        size_t const length1 = MEM_readLE16(istart);
        size_t const length2 = MEM_readLE16(istart+2);
        size_t const length3 = MEM_readLE16(istart+4);
        size_t const length4 = cSrcSize - (length1 + length2 + length3 + 6);
        const BYTE* const istart1 = istart + 6;  /* jumpTable */
        const BYTE* const istart2 = istart1 + length1;
        const BYTE* const istart3 = istart2 + length2;
        const BYTE* const istart4 = istart3 + length3;
        size_t const segmentSize = (dstSize+3) / 4;
        BYTE* const opStart2 = ostart + segmentSize;
        BYTE* const opStart3 = opStart2 + segmentSize;
        BYTE* const opStart4 = opStart3 + segmentSize;
        BYTE* op1 = ostart;
        BYTE* op2 = opStart2;
        BYTE* op3 = opStart3;
        BYTE* op4 = opStart4;
        U32 endSignal;
        DTableDesc const dtd = HUF_getDTableDesc(DTable);
        U32 const dtLog = dtd.tableLog;

        if (length4 > cSrcSize) return ERROR(corruption_detected);   /* overflow */
        { size_t const errorCode = BIT_initDStream(&bitD1, istart1, length1);
          if (HUF_isError(errorCode)) return errorCode; }
        { size_t const errorCode = BIT_initDStream(&bitD2, istart2, length2);
          if (HUF_isError(errorCode)) return errorCode; }
        { size_t const errorCode = BIT_initDStream(&bitD3, istart3, length3);
          if (HUF_isError(errorCode)) return errorCode; }
        { size_t const errorCode = BIT_initDStream(&bitD4, istart4, length4);
          if (HUF_isError(errorCode)) return errorCode; }

        /* 16-32 symbols per loop (4-8 symbols per stream) */
        endSignal = BIT_reloadDStream(&bitD1) | BIT_reloadDStream(&bitD2) | BIT_reloadDStream(&bitD3) | BIT_reloadDStream(&bitD4);
        for ( ; (endSignal==BIT_DStream_unfinished) & (op4<(oend-(sizeof(bitD4.bitContainer)-1))) ; ) {
            HUF_DECODE_SYMBOLX4_2(op1, &bitD1);
            HUF_DECODE_SYMBOLX4_2(op2, &bitD2);
            HUF_DECODE_SYMBOLX4_2(op3, &bitD3);
            HUF_DECODE_SYMBOLX4_2(op4, &bitD4);
            HUF_DECODE_SYMBOLX4_1(op1, &bitD1);
            HUF_DECODE_SYMBOLX4_1(op2, &bitD2);
            HUF_DECODE_SYMBOLX4_1(op3, &bitD3);
            HUF_DECODE_SYMBOLX4_1(op4, &bitD4);
            HUF_DECODE_SYMBOLX4_2(op1, &bitD1);
            HUF_DECODE_SYMBOLX4_2(op2, &bitD2);
            HUF_DECODE_SYMBOLX4_2(op3, &bitD3);
            HUF_DECODE_SYMBOLX4_2(op4, &bitD4);
            HUF_DECODE_SYMBOLX4_0(op1, &bitD1);
            HUF_DECODE_SYMBOLX4_0(op2, &bitD2);
            HUF_DECODE_SYMBOLX4_0(op3, &bitD3);
            HUF_DECODE_SYMBOLX4_0(op4, &bitD4);

            endSignal = BIT_reloadDStream(&bitD1) | BIT_reloadDStream(&bitD2) | BIT_reloadDStream(&bitD3) | BIT_reloadDStream(&bitD4);
        }

        /* check corruption */
        if (op1 > opStart2) return ERROR(corruption_detected);
        if (op2 > opStart3) return ERROR(corruption_detected);
        if (op3 > opStart4) return ERROR(corruption_detected);
        /* note : op4 already verified within main loop */

        /* finish bitStreams one by one */
        HUF_decodeStreamX4(op1, &bitD1, opStart2, dt, dtLog);
        HUF_decodeStreamX4(op2, &bitD2, opStart3, dt, dtLog);
        HUF_decodeStreamX4(op3, &bitD3, opStart4, dt, dtLog);
        HUF_decodeStreamX4(op4, &bitD4, oend,     dt, dtLog);

        /* check */
        { U32 const endCheck = BIT_endOfDStream(&bitD1) & BIT_endOfDStream(&bitD2) & BIT_endOfDStream(&bitD3) & BIT_endOfDStream(&bitD4);
          if (!endCheck) return ERROR(corruption_detected); }

        /* decoded size */
        return dstSize;
    }
}


size_t HUF_decompress4X4_usingDTable(
          void* dst,  size_t dstSize,
    const void* cSrc, size_t cSrcSize,
    const HUF_DTable* DTable)
{
    DTableDesc dtd = HUF_getDTableDesc(DTable);
    if (dtd.tableType != 1) return ERROR(GENERIC);
    return HUF_decompress4X4_usingDTable_internal(dst, dstSize, cSrc, cSrcSize, DTable);
}


size_t HUF_decompress4X4_DCtx (HUF_DTable* dctx, void* dst, size_t dstSize, const void* cSrc, size_t cSrcSize)
{
    const BYTE* ip = (const BYTE*) cSrc;

    size_t hSize = HUF_readDTableX4 (dctx, cSrc, cSrcSize);
    if (HUF_isError(hSize)) return hSize;
    if (hSize >= cSrcSize) return ERROR(srcSize_wrong);
    ip += hSize; cSrcSize -= hSize;

    return HUF_decompress4X4_usingDTable_internal(dst, dstSize, ip, cSrcSize, dctx);
}

size_t HUF_decompress4X4 (void* dst, size_t dstSize, const void* cSrc, size_t cSrcSize)
{
    HUF_CREATE_STATIC_DTABLEX4(DTable, HUF_TABLELOG_MAX);
    return HUF_decompress4X4_DCtx(DTable, dst, dstSize, cSrc, cSrcSize);
}
