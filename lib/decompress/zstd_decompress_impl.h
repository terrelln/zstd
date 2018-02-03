/*
 * Copyright (c) 2018-present, Facebook, Inc.
 * All rights reserved.
 *
 * This source code is licensed under both the BSD-style license (found in the
 * LICENSE file in the root directory of this source tree) and the GPLv2 (found
 * in the COPYING file in the root directory of this source tree).
 * You may select, at your option, one of the above-listed licenses.
 */

#ifndef FUNCTION_NAME
#  error "FUNCTION_NAME(name) must be defined"
#endif

#ifndef TARGET
#  error "TARGET must be defined"
#endif

/* We need to add at most (ZSTD_WINDOWLOG_MAX_32 - 1) bits to read the maximum
 * offset bits. But we can only read at most (STREAM_ACCUMULATOR_MIN_32 - 1)
 * bits before reloading. This value is the maximum number of bytes we read
 * after reloading when we are decoding long offets.
 */
#define LONG_OFFSETS_MAX_EXTRA_BITS_32                                         \
    (ZSTD_WINDOWLOG_MAX_32 > STREAM_ACCUMULATOR_MIN_32                         \
        ? ZSTD_WINDOWLOG_MAX_32 - STREAM_ACCUMULATOR_MIN_32                    \
        : 0)

static TARGET seq_t
FUNCTION_NAME(ZSTD_decodeSequence)(seqState_t* seqState, const ZSTD_longOffset_e longOffsets)
{
    seq_t seq;

    U32 const llCode = FSE_peekSymbol(&seqState->stateLL);
    U32 const mlCode = FSE_peekSymbol(&seqState->stateML);
    U32 const ofCode = FSE_peekSymbol(&seqState->stateOffb);   /* <= MaxOff, by table construction */

    U32 const llBits = LL_bits[llCode];
    U32 const mlBits = ML_bits[mlCode];
    U32 const ofBits = ofCode;
    U32 const totalBits = llBits+mlBits+ofBits;

    static const U32 LL_base[MaxLL+1] = {
                             0,    1,    2,     3,     4,     5,     6,      7,
                             8,    9,   10,    11,    12,    13,    14,     15,
                            16,   18,   20,    22,    24,    28,    32,     40,
                            48,   64, 0x80, 0x100, 0x200, 0x400, 0x800, 0x1000,
                            0x2000, 0x4000, 0x8000, 0x10000 };

    static const U32 ML_base[MaxML+1] = {
                             3,  4,  5,    6,     7,     8,     9,    10,
                            11, 12, 13,   14,    15,    16,    17,    18,
                            19, 20, 21,   22,    23,    24,    25,    26,
                            27, 28, 29,   30,    31,    32,    33,    34,
                            35, 37, 39,   41,    43,    47,    51,    59,
                            67, 83, 99, 0x83, 0x103, 0x203, 0x403, 0x803,
                            0x1003, 0x2003, 0x4003, 0x8003, 0x10003 };

    static const U32 OF_base[MaxOff+1] = {
                     0,        1,       1,       5,     0xD,     0x1D,     0x3D,     0x7D,
                     0xFD,   0x1FD,   0x3FD,   0x7FD,   0xFFD,   0x1FFD,   0x3FFD,   0x7FFD,
                     0xFFFD, 0x1FFFD, 0x3FFFD, 0x7FFFD, 0xFFFFD, 0x1FFFFD, 0x3FFFFD, 0x7FFFFD,
                     0xFFFFFD, 0x1FFFFFD, 0x3FFFFFD, 0x7FFFFFD, 0xFFFFFFD, 0x1FFFFFFD, 0x3FFFFFFD, 0x7FFFFFFD };

    /* sequence */
    {   size_t offset;
        if (!ofCode)
            offset = 0;
        else {
            ZSTD_STATIC_ASSERT(ZSTD_lo_isLongOffset == 1);
            ZSTD_STATIC_ASSERT(LONG_OFFSETS_MAX_EXTRA_BITS_32 == 5);
            assert(ofBits <= MaxOff);
            if (MEM_32bits() && longOffsets) {
                U32 const extraBits = ofBits - MIN(ofBits, STREAM_ACCUMULATOR_MIN_32-1);
                offset = OF_base[ofCode] + (BIT_readBitsFast(&seqState->DStream, ofBits - extraBits) << extraBits);
                if (MEM_32bits() || extraBits) BIT_reloadDStream(&seqState->DStream);
                if (extraBits) offset += BIT_readBitsFast(&seqState->DStream, extraBits);
            } else {
                offset = OF_base[ofCode] + BIT_readBitsFast(&seqState->DStream, ofBits);   /* <=  (ZSTD_WINDOWLOG_MAX-1) bits */
                if (MEM_32bits()) BIT_reloadDStream(&seqState->DStream);
            }
        }

        if (ofCode <= 1) {
            offset += (llCode==0);
            if (offset) {
                size_t temp = (offset==3) ? seqState->prevOffset[0] - 1 : seqState->prevOffset[offset];
                temp += !temp;   /* 0 is not valid; input is corrupted; force offset to 1 */
                if (offset != 1) seqState->prevOffset[2] = seqState->prevOffset[1];
                seqState->prevOffset[1] = seqState->prevOffset[0];
                seqState->prevOffset[0] = offset = temp;
            } else {
                offset = seqState->prevOffset[0];
            }
        } else {
            seqState->prevOffset[2] = seqState->prevOffset[1];
            seqState->prevOffset[1] = seqState->prevOffset[0];
            seqState->prevOffset[0] = offset;
        }
        seq.offset = offset;
    }

    seq.matchLength = ML_base[mlCode]
                    + ((mlCode>31) ? BIT_readBitsFast(&seqState->DStream, mlBits) : 0);  /* <=  16 bits */
    if (MEM_32bits() && (mlBits+llBits >= STREAM_ACCUMULATOR_MIN_32-LONG_OFFSETS_MAX_EXTRA_BITS_32))
        BIT_reloadDStream(&seqState->DStream);
    if (MEM_64bits() && (totalBits >= STREAM_ACCUMULATOR_MIN_64-(LLFSELog+MLFSELog+OffFSELog)))
        BIT_reloadDStream(&seqState->DStream);
    /* Verify that there is enough bits to read the rest of the data in 64-bit mode. */
    ZSTD_STATIC_ASSERT(16+LLFSELog+MLFSELog+OffFSELog < STREAM_ACCUMULATOR_MIN_64);

    seq.litLength = LL_base[llCode]
                  + ((llCode>15) ? BIT_readBitsFast(&seqState->DStream, llBits) : 0);    /* <=  16 bits */
    if (MEM_32bits())
        BIT_reloadDStream(&seqState->DStream);

    DEBUGLOG(6, "seq: litL=%u, matchL=%u, offset=%u",
                (U32)seq.litLength, (U32)seq.matchLength, (U32)seq.offset);

    /* ANS state update */
    FSE_updateState(&seqState->stateLL, &seqState->DStream);    /* <=  9 bits */
    FSE_updateState(&seqState->stateML, &seqState->DStream);    /* <=  9 bits */
    if (MEM_32bits()) BIT_reloadDStream(&seqState->DStream);    /* <= 18 bits */
    FSE_updateState(&seqState->stateOffb, &seqState->DStream);  /* <=  8 bits */

    return seq;
}


HINT_INLINE TARGET seq_t
FUNCTION_NAME(ZSTD_decodeSequenceLong)(seqState_t* seqState, ZSTD_longOffset_e const longOffsets)
{
    seq_t seq;

    U32 const llCode = FSE_peekSymbol(&seqState->stateLL);
    U32 const mlCode = FSE_peekSymbol(&seqState->stateML);
    U32 const ofCode = FSE_peekSymbol(&seqState->stateOffb);   /* <= MaxOff, by table construction */

    U32 const llBits = LL_bits[llCode];
    U32 const mlBits = ML_bits[mlCode];
    U32 const ofBits = ofCode;
    U32 const totalBits = llBits+mlBits+ofBits;

    static const U32 LL_base[MaxLL+1] = {
                             0,  1,    2,     3,     4,     5,     6,      7,
                             8,  9,   10,    11,    12,    13,    14,     15,
                            16, 18,   20,    22,    24,    28,    32,     40,
                            48, 64, 0x80, 0x100, 0x200, 0x400, 0x800, 0x1000,
                            0x2000, 0x4000, 0x8000, 0x10000 };

    static const U32 ML_base[MaxML+1] = {
                             3,  4,  5,    6,     7,     8,     9,    10,
                            11, 12, 13,   14,    15,    16,    17,    18,
                            19, 20, 21,   22,    23,    24,    25,    26,
                            27, 28, 29,   30,    31,    32,    33,    34,
                            35, 37, 39,   41,    43,    47,    51,    59,
                            67, 83, 99, 0x83, 0x103, 0x203, 0x403, 0x803,
                            0x1003, 0x2003, 0x4003, 0x8003, 0x10003 };

    static const U32 OF_base[MaxOff+1] = {
                     0,        1,       1,       5,     0xD,     0x1D,     0x3D,     0x7D,
                     0xFD,   0x1FD,   0x3FD,   0x7FD,   0xFFD,   0x1FFD,   0x3FFD,   0x7FFD,
                     0xFFFD, 0x1FFFD, 0x3FFFD, 0x7FFFD, 0xFFFFD, 0x1FFFFD, 0x3FFFFD, 0x7FFFFD,
                     0xFFFFFD, 0x1FFFFFD, 0x3FFFFFD, 0x7FFFFFD, 0xFFFFFFD, 0x1FFFFFFD, 0x3FFFFFFD, 0x7FFFFFFD };

    /* sequence */
    {   size_t offset;
        if (!ofCode)
            offset = 0;
        else {
            ZSTD_STATIC_ASSERT(ZSTD_lo_isLongOffset == 1);
            ZSTD_STATIC_ASSERT(LONG_OFFSETS_MAX_EXTRA_BITS_32 == 5);
            assert(ofBits <= MaxOff);
            if (MEM_32bits() && longOffsets) {
                U32 const extraBits = ofBits - MIN(ofBits, STREAM_ACCUMULATOR_MIN_32-1);
                offset = OF_base[ofCode] + (BIT_readBitsFast(&seqState->DStream, ofBits - extraBits) << extraBits);
                if (MEM_32bits() || extraBits) BIT_reloadDStream(&seqState->DStream);
                if (extraBits) offset += BIT_readBitsFast(&seqState->DStream, extraBits);
            } else {
                offset = OF_base[ofCode] + BIT_readBitsFast(&seqState->DStream, ofBits);   /* <=  (ZSTD_WINDOWLOG_MAX-1) bits */
                if (MEM_32bits()) BIT_reloadDStream(&seqState->DStream);
            }
        }

        if (ofCode <= 1) {
            offset += (llCode==0);
            if (offset) {
                size_t temp = (offset==3) ? seqState->prevOffset[0] - 1 : seqState->prevOffset[offset];
                temp += !temp;   /* 0 is not valid; input is corrupted; force offset to 1 */
                if (offset != 1) seqState->prevOffset[2] = seqState->prevOffset[1];
                seqState->prevOffset[1] = seqState->prevOffset[0];
                seqState->prevOffset[0] = offset = temp;
            } else {
                offset = seqState->prevOffset[0];
            }
        } else {
            seqState->prevOffset[2] = seqState->prevOffset[1];
            seqState->prevOffset[1] = seqState->prevOffset[0];
            seqState->prevOffset[0] = offset;
        }
        seq.offset = offset;
    }

    seq.matchLength = ML_base[mlCode] + ((mlCode>31) ? BIT_readBitsFast(&seqState->DStream, mlBits) : 0);  /* <=  16 bits */
    if (MEM_32bits() && (mlBits+llBits >= STREAM_ACCUMULATOR_MIN_32-LONG_OFFSETS_MAX_EXTRA_BITS_32))
        BIT_reloadDStream(&seqState->DStream);
    if (MEM_64bits() && (totalBits >= STREAM_ACCUMULATOR_MIN_64-(LLFSELog+MLFSELog+OffFSELog)))
        BIT_reloadDStream(&seqState->DStream);
    /* Verify that there is enough bits to read the rest of the data in 64-bit mode. */
    ZSTD_STATIC_ASSERT(16+LLFSELog+MLFSELog+OffFSELog < STREAM_ACCUMULATOR_MIN_64);

    seq.litLength = LL_base[llCode] + ((llCode>15) ? BIT_readBitsFast(&seqState->DStream, llBits) : 0);    /* <=  16 bits */
    if (MEM_32bits())
        BIT_reloadDStream(&seqState->DStream);

    {   size_t const pos = seqState->pos + seq.litLength;
        const BYTE* const matchBase = (seq.offset > pos) ? seqState->dictEnd : seqState->prefixStart;
        seq.match = matchBase + pos - seq.offset;  /* note : this operation can overflow when seq.offset is really too large, which can only happen when input is corrupted.
                                                    * No consequence though : no memory access will occur, overly large offset will be detected in ZSTD_execSequenceLong() */
        seqState->pos = pos + seq.matchLength;
    }

    /* ANS state update */
    FSE_updateState(&seqState->stateLL, &seqState->DStream);    /* <=  9 bits */
    FSE_updateState(&seqState->stateML, &seqState->DStream);    /* <=  9 bits */
    if (MEM_32bits()) BIT_reloadDStream(&seqState->DStream);    /* <= 18 bits */
    FSE_updateState(&seqState->stateOffb, &seqState->DStream);  /* <=  8 bits */

    return seq;
}

static TARGET
size_t FUNCTION_NAME(ZSTD_decompressSequences)(
                               ZSTD_DCtx* dctx,
                               void* dst, size_t maxDstSize,
                         const void* seqStart, size_t seqSize,
                         const ZSTD_longOffset_e isLongOffset)
{
    const BYTE* ip = (const BYTE*)seqStart;
    const BYTE* const iend = ip + seqSize;
    BYTE* const ostart = (BYTE* const)dst;
    BYTE* const oend = ostart + maxDstSize;
    BYTE* op = ostart;
    const BYTE* litPtr = dctx->litPtr;
    const BYTE* const litEnd = litPtr + dctx->litSize;
    const BYTE* const base = (const BYTE*) (dctx->base);
    const BYTE* const vBase = (const BYTE*) (dctx->vBase);
    const BYTE* const dictEnd = (const BYTE*) (dctx->dictEnd);
    int nbSeq;
    DEBUGLOG(5, "ZSTD_decompressSequences");

    /* Build Decoding Tables */
    {   size_t const seqHSize = ZSTD_decodeSeqHeaders(dctx, &nbSeq, ip, seqSize);
        DEBUGLOG(5, "ZSTD_decodeSeqHeaders: size=%u, nbSeq=%i",
                    (U32)seqHSize, nbSeq);
        if (ZSTD_isError(seqHSize)) return seqHSize;
        ip += seqHSize;
    }

    /* Regen sequences */
    if (nbSeq) {
        seqState_t seqState;
        dctx->fseEntropy = 1;
        { U32 i; for (i=0; i<ZSTD_REP_NUM; i++) seqState.prevOffset[i] = dctx->entropy.rep[i]; }
        CHECK_E(BIT_initDStream(&seqState.DStream, ip, iend-ip), corruption_detected);
        FSE_initDState(&seqState.stateLL, &seqState.DStream, dctx->LLTptr);
        FSE_initDState(&seqState.stateOffb, &seqState.DStream, dctx->OFTptr);
        FSE_initDState(&seqState.stateML, &seqState.DStream, dctx->MLTptr);

        for ( ; (BIT_reloadDStream(&(seqState.DStream)) <= BIT_DStream_completed) && nbSeq ; ) {
            nbSeq--;
            {   seq_t const sequence = FUNCTION_NAME(ZSTD_decodeSequence)(&seqState, isLongOffset);
                size_t const oneSeqSize = ZSTD_execSequence(op, oend, sequence, &litPtr, litEnd, base, vBase, dictEnd);
                DEBUGLOG(6, "regenerated sequence size : %u", (U32)oneSeqSize);
                if (ZSTD_isError(oneSeqSize)) return oneSeqSize;
                op += oneSeqSize;
        }   }

        /* check if reached exact end */
        DEBUGLOG(5, "after decode loop, remaining nbSeq : %i", nbSeq);
        if (nbSeq) return ERROR(corruption_detected);
        /* save reps for next block */
        { U32 i; for (i=0; i<ZSTD_REP_NUM; i++) dctx->entropy.rep[i] = (U32)(seqState.prevOffset[i]); }
    }

    /* last literal segment */
    {   size_t const lastLLSize = litEnd - litPtr;
        if (lastLLSize > (size_t)(oend-op)) return ERROR(dstSize_tooSmall);
        memcpy(op, litPtr, lastLLSize);
        op += lastLLSize;
    }

    return op-ostart;
}

static TARGET
size_t FUNCTION_NAME(ZSTD_decompressSequencesLong)(
                               ZSTD_DCtx* dctx,
                               void* dst, size_t maxDstSize,
                         const void* seqStart, size_t seqSize,
                         const ZSTD_longOffset_e isLongOffset)
{
    const BYTE* ip = (const BYTE*)seqStart;
    const BYTE* const iend = ip + seqSize;
    BYTE* const ostart = (BYTE* const)dst;
    BYTE* const oend = ostart + maxDstSize;
    BYTE* op = ostart;
    const BYTE* litPtr = dctx->litPtr;
    const BYTE* const litEnd = litPtr + dctx->litSize;
    const BYTE* const prefixStart = (const BYTE*) (dctx->base);
    const BYTE* const dictStart = (const BYTE*) (dctx->vBase);
    const BYTE* const dictEnd = (const BYTE*) (dctx->dictEnd);
    int nbSeq;

    /* Build Decoding Tables */
    {   size_t const seqHSize = ZSTD_decodeSeqHeaders(dctx, &nbSeq, ip, seqSize);
        if (ZSTD_isError(seqHSize)) return seqHSize;
        ip += seqHSize;
    }

    /* Regen sequences */
    if (nbSeq) {
#define STORED_SEQS 4
#define STOSEQ_MASK (STORED_SEQS-1)
#define ADVANCED_SEQS 4
        seq_t sequences[STORED_SEQS];
        int const seqAdvance = MIN(nbSeq, ADVANCED_SEQS);
        seqState_t seqState;
        int seqNb;
        dctx->fseEntropy = 1;
        { U32 i; for (i=0; i<ZSTD_REP_NUM; i++) seqState.prevOffset[i] = dctx->entropy.rep[i]; }
        seqState.prefixStart = prefixStart;
        seqState.pos = (size_t)(op-prefixStart);
        seqState.dictEnd = dictEnd;
        CHECK_E(BIT_initDStream(&seqState.DStream, ip, iend-ip), corruption_detected);
        FSE_initDState(&seqState.stateLL, &seqState.DStream, dctx->LLTptr);
        FSE_initDState(&seqState.stateOffb, &seqState.DStream, dctx->OFTptr);
        FSE_initDState(&seqState.stateML, &seqState.DStream, dctx->MLTptr);

        /* prepare in advance */
        for (seqNb=0; (BIT_reloadDStream(&seqState.DStream) <= BIT_DStream_completed) && seqNb<seqAdvance; seqNb++) {
            sequences[seqNb] = FUNCTION_NAME(ZSTD_decodeSequenceLong)(&seqState, isLongOffset);
        }
        if (seqNb<seqAdvance) return ERROR(corruption_detected);

        /* decode and decompress */
        for ( ; (BIT_reloadDStream(&(seqState.DStream)) <= BIT_DStream_completed) && seqNb<nbSeq ; seqNb++) {
            seq_t const sequence = FUNCTION_NAME(ZSTD_decodeSequenceLong)(&seqState, isLongOffset);
            size_t const oneSeqSize = ZSTD_execSequenceLong(op, oend, sequences[(seqNb-ADVANCED_SEQS) & STOSEQ_MASK], &litPtr, litEnd, prefixStart, dictStart, dictEnd);
            if (ZSTD_isError(oneSeqSize)) return oneSeqSize;
            PREFETCH(sequence.match);  /* note : it's safe to invoke PREFETCH() on any memory address, including invalid ones */
            sequences[seqNb&STOSEQ_MASK] = sequence;
            op += oneSeqSize;
        }
        if (seqNb<nbSeq) return ERROR(corruption_detected);

        /* finish queue */
        seqNb -= seqAdvance;
        for ( ; seqNb<nbSeq ; seqNb++) {
            size_t const oneSeqSize = ZSTD_execSequenceLong(op, oend, sequences[seqNb&STOSEQ_MASK], &litPtr, litEnd, prefixStart, dictStart, dictEnd);
            if (ZSTD_isError(oneSeqSize)) return oneSeqSize;
            op += oneSeqSize;
        }

        /* save reps for next block */
        { U32 i; for (i=0; i<ZSTD_REP_NUM; i++) dctx->entropy.rep[i] = (U32)(seqState.prevOffset[i]); }
#undef STORED_SEQS
#undef STOSEQ_MASK
#undef ADVANCED_SEQS
    }

    /* last literal segment */
    {   size_t const lastLLSize = litEnd - litPtr;
        if (lastLLSize > (size_t)(oend-op)) return ERROR(dstSize_tooSmall);
        memcpy(op, litPtr, lastLLSize);
        op += lastLLSize;
    }

    return op-ostart;
}
