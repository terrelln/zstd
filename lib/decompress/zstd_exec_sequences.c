#define ZSTD_STATIC_LINKING_ONLY
#include "zstd.h"
#include "zstd_internal.h"

static void ZSTD_copy4(void* dst, const void* src) { memcpy(dst, src, 4); }


FORCE_NOINLINE
size_t ZSTD_execSequenceLast7(BYTE* op,
                              BYTE* const oend, seq_t sequence,
                              const BYTE** litPtr, const BYTE* const litLimit,
                              const BYTE* const base, const BYTE* const vBase, const BYTE* const dictEnd)
{
    BYTE* const oLitEnd = op + sequence.litLength;
    size_t const sequenceLength = sequence.litLength + sequence.matchLength;
    BYTE* const oMatchEnd = op + sequenceLength;   /* risk : address space overflow (32-bits) */
    BYTE* const oend_w = oend - WILDCOPY_OVERLENGTH;
    const BYTE* const iLitEnd = *litPtr + sequence.litLength;
    const BYTE* match = oLitEnd - sequence.offset;

    /* check */
    if (oMatchEnd>oend) return ERROR(dstSize_tooSmall); /* last match must start at a minimum distance of WILDCOPY_OVERLENGTH from oend */
    if (iLitEnd > litLimit) return ERROR(corruption_detected);   /* over-read beyond lit buffer */
    if (oLitEnd <= oend_w) return ERROR(GENERIC);   /* Precondition */

    /* copy literals */
    if (op < oend_w) {
        ZSTD_wildcopy(op, *litPtr, oend_w - op);
        *litPtr += oend_w - op;
        op = oend_w;
    }
    while (op < oLitEnd) *op++ = *(*litPtr)++;

    /* copy Match */
    if (sequence.offset > (size_t)(oLitEnd - base)) {
        /* offset beyond prefix */
        if (sequence.offset > (size_t)(oLitEnd - vBase)) return ERROR(corruption_detected);
        match = dictEnd - (base-match);
        if (match + sequence.matchLength <= dictEnd) {
            memmove(oLitEnd, match, sequence.matchLength);
            return sequenceLength;
        }
        /* span extDict & currentPrefixSegment */
        {   size_t const length1 = dictEnd - match;
            memmove(oLitEnd, match, length1);
            op = oLitEnd + length1;
            sequence.matchLength -= length1;
            match = base;
    }   }
    while (op < oMatchEnd) *op++ = *match++;
    return sequenceLength;
}


size_t ZSTD_execSequence(BYTE* op,
                                BYTE* const oend, seq_t sequence,
                                const BYTE** litPtr, const BYTE* const litLimit,
                                const BYTE* const base, const BYTE* const vBase, const BYTE* const dictEnd)
{
    BYTE* const oLitEnd = op + sequence.litLength;
    size_t const sequenceLength = sequence.litLength + sequence.matchLength;
    BYTE* const oMatchEnd = op + sequenceLength;   /* risk : address space overflow (32-bits) */
    BYTE* const oend_w = oend - WILDCOPY_OVERLENGTH;
    const BYTE* const iLitEnd = *litPtr + sequence.litLength;
    const BYTE* match = oLitEnd - sequence.offset;

    /* check */
    if (oMatchEnd>oend) return ERROR(dstSize_tooSmall); /* last match must start at a minimum distance of WILDCOPY_OVERLENGTH from oend */
    if (iLitEnd > litLimit) return ERROR(corruption_detected);   /* over-read beyond lit buffer */
    if (oLitEnd>oend_w) return ZSTD_execSequenceLast7(op, oend, sequence, litPtr, litLimit, base, vBase, dictEnd);

    /* copy Literals */
    ZSTD_copy8(op, *litPtr);
    if (sequence.litLength > 8)
        ZSTD_wildcopy(op+8, (*litPtr)+8, sequence.litLength - 8);   /* note : since oLitEnd <= oend-WILDCOPY_OVERLENGTH, no risk of overwrite beyond oend */
    op = oLitEnd;
    *litPtr = iLitEnd;   /* update for next sequence */

    /* copy Match */
    if (sequence.offset > (size_t)(oLitEnd - base)) {
        /* offset beyond prefix */
        if (sequence.offset > (size_t)(oLitEnd - vBase)) return ERROR(corruption_detected);
        match = dictEnd - (base-match);
        if (match + sequence.matchLength <= dictEnd) {
            memmove(oLitEnd, match, sequence.matchLength);
            return sequenceLength;
        }
        /* span extDict & currentPrefixSegment */
        {   size_t const length1 = dictEnd - match;
            memmove(oLitEnd, match, length1);
            op = oLitEnd + length1;
            sequence.matchLength -= length1;
            match = base;
            if (op > oend_w) {
              U32 i;
              for (i = 0; i < sequence.matchLength; ++i) op[i] = match[i];
              return sequenceLength;
            }
    }   }
    /* Requirement: op <= oend_w */

    /* match within prefix */
    if (sequence.offset < 8) {
        /* close range match, overlap */
        static const U32 dec32table[] = { 0, 1, 2, 1, 4, 4, 4, 4 };   /* added */
        static const int dec64table[] = { 8, 8, 8, 7, 8, 9,10,11 };   /* substracted */
        int const sub2 = dec64table[sequence.offset];
        op[0] = match[0];
        op[1] = match[1];
        op[2] = match[2];
        op[3] = match[3];
        match += dec32table[sequence.offset];
        ZSTD_copy4(op+4, match);
        match -= sub2;
    } else {
        ZSTD_copy8(op, match);
    }
    op += 8; match += 8;

    if (oMatchEnd > oend-(16-MINMATCH)) {
        if (op < oend_w) {
            ZSTD_wildcopy(op, match, oend_w - op);
            match += oend_w - op;
            op = oend_w;
        }
        while (op < oMatchEnd) *op++ = *match++;
    } else {
        ZSTD_wildcopy(op, match, sequence.matchLength-8);   /* works even if matchLength < 8 */
    }
    return sequenceLength;
}
