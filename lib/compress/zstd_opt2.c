#if 1

#include "hist.h"
#include "zstd_internal.h"
#include "zstd_compress_internal.h"
#include "zstd_opt2.h"

#ifndef ZSTD_OPT_MAX_CHUNK
#define ZSTD_OPT_MAX_CHUNK ZSTD_OPT_NUM
#endif

#ifndef ZSTD_OPT_MIN_MATCH
#define ZSTD_OPT_MIN_MATCH 3
#endif

#define ZSTD_OPT_MAX_PRICE ((uint32_t)(1u << 30))

#define NT 3

// typedef struct repcodes_s {
//     U32 rep[ZSTD_REP_NUM];
// } ZSTD_OptRepcode;

// typedef enum {
//   ZSTD_OCT_literals,
//   ZSTD_OCT_match,
// } ZSTD_OptCommandType;

// typedef struct {
//   uint32_t length;
//   uint32_t offset;
//   char type;
// } ZSTD_OptCommand;

// typedef struct {
//   ZSTD_OptCommand const* begin;
//   ZSTD_OptCommand const* end;
// } ZSTD_OptCommands;

// typedef struct {
//   ZSTD_OptCommand cmd;
//   uint32_t price;
//   ZSTD_OptRepcode reps;
// } ZSTD_OptState;

// typedef struct {
//   // PARAMS:MATCHFINDING
//   ZSTD_OptParams params;
//   U32 nextToUpdate3;
//   ZSTD_OptCommand cmds[ZSTD_OPT_MAX_CHUNK + 1];
//   // PARAMS:PRICING
//   unsigned litFreq[256];             /* table of literals statistics, of size 256 */
//   unsigned litLengthFreq[MaxLL+1];   /* table of litLength statistics, of size (MaxLL+1) */
//   unsigned matchLengthFreq[MaxML+1]; /* table of matchLength statistics, of size (MaxML+1) */
//   unsigned offCodeFreq[MaxOff+1];    /* table of offCode statistics, of size (MaxOff+1) */
//   ZSTD_OptState opt[ZSTD_OPT_MAX_CHUNK];

//   U32  litSum;                 /* nb of literals */
//   U32  litLengthSum;           /* nb of litLength codes */
//   U32  matchLengthSum;         /* nb of matchLength codes */
//   U32  offCodeSum;             /* nb of offset codes */
//   U32  litSumBasePrice;        /* to compare to log2(litfreq) */
//   U32  litLengthSumBasePrice;  /* to compare to log2(llfreq)  */
//   U32  matchLengthSumBasePrice;/* to compare to log2(mlfreq)  */
//   U32  offCodeSumBasePrice;    /* to compare to log2(offreq)  */
//   ZSTD_OptPrice_e priceType;   /* prices can be determined dynamically, or follow a pre-defined cost structure */
//   const ZSTD_entropyCTables_t* symbolCosts;  /* pre-calculated dictionary statistics */
//   ZSTD_literalCompressionMode_e literalCompressionMode;

//   ZSTD_matchState_t* ms;
//   ZSTD_OptRepcode initialReps;
// } ZSTD_OptContext;

// typedef struct {
//   U32* hash;
//   U32* bt;
//   U32 hashLog;
//   U32 btLog;
//   U32 btMask;
// } ZSTD_OptBt;

// typedef struct {
//   BYTE const* base;
//   BYTE const* end;
//   U32 highLimit;
//   U32 lowLimit;
//   ZSTD_OptBt bt;
//   U32 btLow;
// } ZSTD_OptDms;

// typedef struct {
//   ZSTD_dictMode_e dictMode; // Templated
//   U32 minMatch; // Templated
//   U32 mls;
//   U32 sufficientLen;
//   U32 windowLog;
//   BYTE const* base;
//   BYTE const* dictBase;
//   U32 dictLimit;
//   BYTE const* dictEnd;
//   BYTE const* prefixStart;
//   U32 nbCompares;
//   ZSTD_OptBt bt;
//   ZSTD_OptDms dms;
//   int optLevel;
// } ZSTD_OptParams;

void ZSTD_OptContext_init(ZSTD_OptContext* ctx, ZSTD_matchState_t* ms) {
  memset(ctx, 0, sizeof(*ctx));
  ctx->ms = ms;
  // ctx->symbolCosts = symbolCosts;
  // ctx->literalCompressionMode = lcm;
}

////// SECTION:MATCHFINDING

/* ZSTD_readMINMATCH() :
 * function safe only for comparisons
 * assumption : memPtr must be at least 4 bytes before end of buffer */
MEM_STATIC U32 ZSTD_readMINMATCH(const void* memPtr, U32 length)
{
    switch (length)
    {
    default :
    case 4 : return MEM_read32(memPtr);
    case 3 : if (MEM_isLittleEndian())
                return MEM_read32(memPtr)<<8;
             else
                return MEM_read32(memPtr)>>8;
    }
}

/* Update hashTable3 up to ip (excluded)
   Assumption : always within prefix (i.e. not within extDict) */
static U32 ZSTD_insertAndFindFirstIndexHash3 (ZSTD_matchState_t* ms,
                                              U32* nextToUpdate3,
                                              const BYTE* const ip)
{
    U32* const hashTable3 = ms->hashTable3;
    U32 const hashLog3 = ms->hashLog3;
    const BYTE* const base = ms->window.base;
    U32 idx = *nextToUpdate3;
    U32 const target = (U32)(ip - base);
    size_t const hash3 = ZSTD_hash3Ptr(ip, hashLog3);
    assert(hashLog3 > 0);

    while(idx < target) {
        hashTable3[ZSTD_hash3Ptr(base+idx, hashLog3)] = idx;
        idx++;
    }

    *nextToUpdate3 = target;
    return hashTable3[hash3];
}

static void fillOptParams(ZSTD_OptParams* params, ZSTD_matchState_t const* ms, ZSTD_dictMode_e dictMode, int optLevel) {
  memset(params, 0, sizeof(*params));
  params->optLevel = optLevel;
  params->dictMode = dictMode;
  params->minMatch = ms->cParams.minMatch == 3 ? 3 : 4;
  params->mls = ms->cParams.minMatch;
  params->windowLog = ms->cParams.windowLog;
  params->sufficientLen = MIN(ms->cParams.targetLength, ZSTD_OPT_MAX_CHUNK);
  params->base = ms->window.base;
  params->dictBase = ms->window.dictBase;
  params->dictLimit = ms->window.dictLimit;
  params->dictEnd = params->dictBase + params->dictLimit;
  params->prefixStart = params->base + params->dictLimit;
  params->nbCompares = 1u << ms->cParams.searchLog;

  params->bt.hash = ms->hashTable;
  params->bt.bt = ms->chainTable;
  params->bt.hashLog = ms->cParams.hashLog;
  params->bt.btLog = ms->cParams.chainLog - 1;
  params->bt.btMask = (1u << params->bt.btLog) - 1;

  if (dictMode == ZSTD_dictMatchState) {
    ZSTD_matchState_t const* dms = ms->dictMatchState;
    params->dms.base             = dms->window.base;
    params->dms.end              = dms->window.nextSrc;
    params->dms.highLimit        = (U32)(params->dms.end - params->dms.base);
    params->dms.lowLimit         = dms->window.lowLimit;
    params->dms.bt.hash          = dms->hashTable;
    params->dms.bt.bt            = dms->chainTable;
    params->dms.bt.hashLog       = dms->cParams.hashLog;
    params->dms.bt.btLog         = dms->cParams.chainLog - 1;
    params->dms.bt.btMask        = (1U << params->dms.bt.btLog) - 1;
    params->dms.btLow            =
        params->dms.bt.btMask < params->dms.highLimit - params->dms.lowLimit
            ? params->dms.highLimit - params->dms.bt.btMask
            : params->dms.lowLimit;
  }
}

static size_t ZSTD_Opt_getLiteralCommand(ZSTD_OptCommand* cmds) {
  cmds[0].length = 1;
  cmds[0].type = ZSTD_OCT_literals;
  return 1;
}

FORCE_INLINE_TEMPLATE size_t ZSTD_Opt_getRepCommands(
    ZSTD_OptContext const* ctx,
    ZSTD_OptCommand* cmds,
    size_t cur,
    uint8_t const* ip,
    uint8_t const* iLimit,
    size_t* bestLength,
    U32 const minMatch /* template */,
    ZSTD_dictMode_e const dictMode /* template */) {
  ZSTD_OptParams const* p = &ctx->params;
  U32 const current = (U32)(ip - p->base);
  U32 const windowLow = ZSTD_getLowestMatchIndex(ctx->ms, current, p->windowLog);
  U32 const dmsIndexDelta = dictMode == ZSTD_dictMatchState ? windowLow - p->dms.highLimit : 0;

  int const ll0 = ctx->opt[cur].cmd.type != ZSTD_OCT_literals || (cur == 0 && ctx->opt[cur].cmd.length == 0);
  U32 const* rep = ctx->opt[cur].reps.rep;
  // TODO: Handle position 0 better
  int const repEnd = ZSTD_REP_NUM + ll0;
  int repCode;
  size_t cnum = 0;
  for (repCode = ll0; repCode < repEnd; ++repCode) {
    U32 const repOffset = (repCode==ZSTD_REP_NUM) ? (rep[0] - 1) : rep[repCode];
    U32 const repIndex = current - repOffset;
    U32 repLen = 0;
    assert(current >= p->dictLimit);
    if (repOffset-1 /* intentional overflow, discards 0 and -1 */ < current - p->dictLimit) {  /* equivalent to `current > repIndex >= dictLimit` */
      if (ZSTD_readMINMATCH(ip, minMatch) == ZSTD_readMINMATCH(ip - repOffset, minMatch)) {
        repLen = (U32)ZSTD_count(ip+minMatch, ip+minMatch-repOffset, iLimit) + minMatch;
      }
    } else {  /* repIndex < dictLimit || repIndex >= current */
      const BYTE* const repMatch = dictMode == ZSTD_dictMatchState ?
                                    p->dms.base + repIndex - dmsIndexDelta :
                                    p->dictBase + repIndex;
      assert(current >= windowLow);
      if ( dictMode == ZSTD_extDict
        && ( ((repOffset-1) /*intentional overflow*/ < current - windowLow)  /* equivalent to `current > repIndex >= windowLow` */
            & (((U32)((p->dictLimit-1) - repIndex) >= 3) ) /* intentional overflow : do not test positions overlapping 2 memory segments */)
        && (ZSTD_readMINMATCH(ip, minMatch) == ZSTD_readMINMATCH(repMatch, minMatch)) ) {
        repLen = (U32)ZSTD_count_2segments(ip+minMatch, repMatch+minMatch, iLimit, p->dictEnd, p->prefixStart) + minMatch;
      }
      if (dictMode == ZSTD_dictMatchState
        && ( ((repOffset-1) /*intentional overflow*/ < current - (p->dms.lowLimit + dmsIndexDelta))  /* equivalent to `current > repIndex >= dmsLowLimit` */
            & ((U32)((p->dictLimit-1) - repIndex) >= 3) ) /* intentional overflow : do not test positions overlapping 2 memory segments */
        && (ZSTD_readMINMATCH(ip, minMatch) == ZSTD_readMINMATCH(repMatch, minMatch)) ) {
        repLen = (U32)ZSTD_count_2segments(ip+minMatch, repMatch+minMatch, iLimit, p->dms.end, p->prefixStart) + minMatch;
      }
    }
    /* save longer solution */
    if (repLen > *bestLength) {
      DEBUGLOG(8, "found repCode %u (ll0:%u, offset:%u) of length %u",
                  repCode, ll0, repOffset, repLen);
      *bestLength = repLen;
      cmds[cnum].offset = repCode - ll0;
      cmds[cnum].length = (U32)repLen;
      cmds[cnum].type = ZSTD_OCT_match;
      cnum++;
      if ( (repLen > p->sufficientLen)
          | (ip+repLen == iLimit) ) {  /* best possible */
        break;
      }
    }
  }
  return cnum;
}

static ZSTD_OptRepcode ZSTD_updateRep(U32 const rep[3], U32 const offset, U32 const ll0)
{
    ZSTD_OptRepcode newReps;
    if (offset >= ZSTD_REP_NUM) {  /* full offset */
        newReps.rep[2] = rep[1];
        newReps.rep[1] = rep[0];
        newReps.rep[0] = offset - ZSTD_REP_MOVE;
    } else {   /* repcode */
        U32 const repCode = offset + ll0;
        if (repCode > 0) {  /* note : if repCode==0, no change */
            U32 const currentOffset = (repCode==ZSTD_REP_NUM) ? (rep[0] - 1) : rep[repCode];
            newReps.rep[2] = (repCode >= 2) ? rep[1] : rep[2];
            newReps.rep[1] = rep[0];
            newReps.rep[0] = currentOffset;
        } else {   /* repCode == 0 */
            memcpy(&newReps, rep, sizeof(newReps));
        }
    }
    return newReps;
}

FORCE_INLINE_TEMPLATE size_t ZSTD_Opt_getMLRCommands(
    ZSTD_OptContext const* ctx,
    ZSTD_OptCommand* cmds,
    size_t cur,
    uint8_t const* const istart,
    uint8_t const* iLimit,
    size_t* bestLength,
    U32 const minMatch /* template */,
    ZSTD_dictMode_e const dictMode /* template */) {
  const size_t kMaxLits = 8;
  ZSTD_OptParams const* p = &ctx->params;
  ZSTD_OptCommand const* cmd;
  ZSTD_OptCommand const* const cmdEnd = cmds;
  int const curll0 = ctx->opt[cur].cmd.type != ZSTD_OCT_literals || (cur == 0 && ctx->opt[cur].cmd.length == 0);
  size_t cnum = 0;
  for (cmd = ctx->cmds; cmd != cmdEnd; ++cmd) {
    uint8_t const* ip = istart + cmd->length;
    if (ip > iLimit - kMaxLits - minMatch) {
      return 0;
    }
    if (cmd->type == ZSTD_OCT_literals) {
      continue;
    }
    ZSTD_OptRepcode reps = ZSTD_updateRep(ctx->opt[cur].reps.rep, cmd->offset, curll0);
    for (size_t lits = 0; lits < kMaxLits + 1; ++lits, ++ip) {
      U32 const current = (U32)(ip - p->base);
      U32 const windowLow = ZSTD_getLowestMatchIndex(ctx->ms, current, p->windowLog);
      U32 const dmsIndexDelta = dictMode == ZSTD_dictMatchState ? windowLow - p->dms.highLimit : 0;
      int const ll0 = lits == 0;
      U32 const* rep = reps.rep;
      // TODO: Handle position 0 better
      int const repEnd = ZSTD_REP_NUM + ll0;
      int repCode;
      for (repCode = ll0; repCode < repEnd; ++repCode) {
        U32 const repOffset = (repCode==ZSTD_REP_NUM) ? (rep[0] - 1) : rep[repCode];
        U32 const repIndex = current - repOffset;
        U32 repLen = 0;
        assert(current >= p->dictLimit);
        if (repOffset-1 /* intentional overflow, discards 0 and -1 */ < current - p->dictLimit) {  /* equivalent to `current > repIndex >= dictLimit` */
          if (ZSTD_readMINMATCH(ip, minMatch) == ZSTD_readMINMATCH(ip - repOffset, minMatch)) {
            repLen = (U32)ZSTD_count(ip+minMatch, ip+minMatch-repOffset, iLimit) + minMatch;
          }
        } else {  /* repIndex < dictLimit || repIndex >= current */
          const BYTE* const repMatch = dictMode == ZSTD_dictMatchState ?
                                        p->dms.base + repIndex - dmsIndexDelta :
                                        p->dictBase + repIndex;
          assert(current >= windowLow);
          if ( dictMode == ZSTD_extDict
            && ( ((repOffset-1) /*intentional overflow*/ < current - windowLow)  /* equivalent to `current > repIndex >= windowLow` */
                & (((U32)((p->dictLimit-1) - repIndex) >= 3) ) /* intentional overflow : do not test positions overlapping 2 memory segments */)
            && (ZSTD_readMINMATCH(ip, minMatch) == ZSTD_readMINMATCH(repMatch, minMatch)) ) {
            repLen = (U32)ZSTD_count_2segments(ip+minMatch, repMatch+minMatch, iLimit, p->dictEnd, p->prefixStart) + minMatch;
          }
          if (dictMode == ZSTD_dictMatchState
            && ( ((repOffset-1) /*intentional overflow*/ < current - (p->dms.lowLimit + dmsIndexDelta))  /* equivalent to `current > repIndex >= dmsLowLimit` */
                & ((U32)((p->dictLimit-1) - repIndex) >= 3) ) /* intentional overflow : do not test positions overlapping 2 memory segments */
            && (ZSTD_readMINMATCH(ip, minMatch) == ZSTD_readMINMATCH(repMatch, minMatch)) ) {
            repLen = (U32)ZSTD_count_2segments(ip+minMatch, repMatch+minMatch, iLimit, p->dms.end, p->prefixStart) + minMatch;
          }
        }
        /* save longer solution */
        repLen = MIN(255, repLen);
        size_t const totalLen = cmd->length + lits + repLen;
        if (repLen >= minMatch && totalLen > *bestLength) {
          assert(repLen >= MINMATCH);
          DEBUGLOG(8, "found repCode %u (ll0:%u, offset:%u) of length %u",
                      repCode, ll0, repOffset, repLen);
          *bestLength = totalLen;
          cmds[cnum] = *cmd;
          cmds[cnum].type = ZSTD_OCT_mlr;
          cmds[cnum].lits = lits;
          cmds[cnum].rep = repCode - ll0;
          cmds[cnum].repLen = repLen;
          cmds[cnum].length = totalLen;
          cnum++;
          if ( (totalLen > p->sufficientLen) | (repLen >= 255)
              | (ip+repLen == iLimit) ) {  /* best possible */
            break;
          }
        }
      }
    }
  }
  return cnum;
}

FORCE_INLINE_TEMPLATE size_t ZSTD_Opt_getHC3Commands(
    ZSTD_OptContext* ctx,
    ZSTD_OptCommand* cmds,
    uint8_t const* ip,
    uint8_t const* iLimit,
    size_t* bestLength,
    ZSTD_dictMode_e const dictMode /* template */) {
  ZSTD_OptParams const* p = &ctx->params;
  U32 const matchIndex3 = ZSTD_insertAndFindFirstIndexHash3(ctx->ms, &ctx->nextToUpdate3, ip);
  U32 const current = (U32)(ip - p->base);
  U32 const windowLow = ZSTD_getLowestMatchIndex(ctx->ms, current, p->windowLog);
  U32 const matchLow = windowLow ? windowLow : 1;

  assert(p->minMatch == 3);
  assert(*bestLength < 3);

  // Heuristic: Ignore longer distances, likely too expensive
  if ((matchIndex3 >= matchLow) & (current - matchIndex3 < (1<<18))) {
    size_t mlen;
    if ((dictMode == ZSTD_noDict) || (dictMode == ZSTD_dictMatchState) || (matchIndex3 >= p->dictLimit)) {
      BYTE const* const match = p->base + matchIndex3;
      mlen = ZSTD_count(ip, match, iLimit);
    } else {
      BYTE const* const match = p->dictBase + matchIndex3;
      mlen = ZSTD_count_2segments(ip, match, iLimit, p->dictEnd, p->prefixStart);
    }

    /* save best solution */
    if (mlen >= 3) {
      DEBUGLOG(8, "found small match with hlog3, of length %u",
                  (U32)mlen);
      *bestLength = mlen;
      assert(current > matchIndex3);
      cmds[0].offset = (current - matchIndex3) + ZSTD_REP_MOVE;
      cmds[0].length = (U32)mlen;
      cmds[0].type = ZSTD_OCT_match;
      if ( (mlen > p->sufficientLen) |
            (ip+mlen == iLimit) ) {  /* best possible length */
        // TODO: How do we want to handle this?
        // ms->nextToUpdate = current+1;  /* skip insertion */
      }
      return 1;
    }
  }
  /* no dictMatchState lookup: dicts don't have a populated HC3 table */
  return 0;
}

static size_t ZSTD_Opt_getBtCommands(
    ZSTD_OptContext const* ctx,
    ZSTD_OptCommand* cmds,
    uint8_t const* ip,
    uint8_t const* iLimit,
    size_t* bestLength,
    ZSTD_dictMode_e const dictMode /* template */) {
  ZSTD_OptParams const* p = &ctx->params;
  U32 const current = (U32)(ip - p->base);
  U32 const windowLow = ZSTD_getLowestMatchIndex(ctx->ms, current, p->windowLog);
  U32 const matchLow = windowLow ? windowLow : 1;
  size_t const h  = ZSTD_hashPtr(ip, p->bt.hashLog, p->mls);
  U32 matchIndex  = p->bt.hash[h];
  U32 nbCompares  = p->nbCompares;
  size_t commonLengthSmaller = 0;
  size_t commonLengthLarger = 0;
  U32 matchEndIdx = current+8+1;   /* farthest referenced position of any match => detects repetitive patterns */
  size_t cnum = 0;
  U32* const bt = p->bt.bt;
  U32 const btMask = p->bt.btMask;
  U32 const btLow = (btMask >= current) ? 0 : current - btMask;
  U32* smallerPtr = bt + 2*(current&btMask);
  U32* largerPtr  = bt + 2*(current&btMask) + 1;
  U32 dummy32;
  DEBUGLOG(NT, "current = %u", current);

  p->bt.hash[h] = current;   /* Update Hash Table */

  while (nbCompares-- && (matchIndex >= matchLow)) {
    U32* const nextPtr = bt + 2*(matchIndex & btMask);
    const BYTE* match;
    size_t matchLength = MIN(commonLengthSmaller, commonLengthLarger);   /* guaranteed minimum nb of common bytes */
    assert(current > matchIndex);

    if ((dictMode == ZSTD_noDict) || (dictMode == ZSTD_dictMatchState) || (matchIndex+matchLength >= p->dictLimit)) {
        assert(matchIndex+matchLength >= p->dictLimit);  /* ensure the condition is correct when !extDict */
        match = p->base + matchIndex;
        if (matchIndex >= p->dictLimit) assert(memcmp(match, ip, matchLength) == 0);  /* ensure early section of match is equal as expected */
        matchLength += ZSTD_count(ip+matchLength, match+matchLength, iLimit);
    } else {
        match = p->dictBase + matchIndex;
        assert(memcmp(match, ip, matchLength) == 0);  /* ensure early section of match is equal as expected */
        matchLength += ZSTD_count_2segments(ip+matchLength, match+matchLength, iLimit, p->dictEnd, p->prefixStart);
        if (matchIndex+matchLength >= p->dictLimit)
            match = p->base + matchIndex;   /* prepare for match[matchLength] read */
    }

    if (matchLength > *bestLength) {
        DEBUGLOG(8, "found match of length %u at distance %u (offCode=%u)",
                (U32)matchLength, current - matchIndex, current - matchIndex + ZSTD_REP_MOVE);
        assert(matchEndIdx > matchIndex);
        if (matchLength > matchEndIdx - matchIndex)
            matchEndIdx = matchIndex + (U32)matchLength;
        *bestLength = matchLength;
        cmds[cnum].offset = (current - matchIndex) + ZSTD_REP_MOVE;
        cmds[cnum].length = (U32)matchLength;
        cmds[cnum].type = ZSTD_OCT_match;
        cnum++;
        if ( (matchLength > ZSTD_OPT_MAX_CHUNK)
            | (ip+matchLength == iLimit) /* equal : no way to know if inf or sup */) {
            if (dictMode == ZSTD_dictMatchState) nbCompares = 0; /* break should also skip searching dms */
            break; /* drop, to preserve bt consistency (miss a little bit of compression) */
        }
    }

    if (match[matchLength] < ip[matchLength]) {
      /* match smaller than current */
      *smallerPtr = matchIndex;             /* update smaller idx */
      commonLengthSmaller = matchLength;    /* all smaller will now have at least this guaranteed common length */
      if (matchIndex <= btLow) { smallerPtr=&dummy32; break; }   /* beyond tree size, stop the search */
      smallerPtr = nextPtr+1;               /* new candidate => larger than match, which was smaller than current */
      matchIndex = nextPtr[1];              /* new matchIndex, larger than previous, closer to current */
    } else {
      *largerPtr = matchIndex;
      commonLengthLarger = matchLength;
      if (matchIndex <= btLow) { largerPtr=&dummy32; break; }   /* beyond tree size, stop the search */
      largerPtr = nextPtr;
      matchIndex = nextPtr[0];
    }
  }

  *smallerPtr = *largerPtr = 0;

  if (dictMode == ZSTD_dictMatchState && nbCompares) {
    ZSTD_OptDms const* dms = &p->dms;
    U32* const dmsBt = dms->bt.bt;
    U32 const dmsBtMask = dms->bt.btMask;
    size_t const dmsH = ZSTD_hashPtr(ip, dms->bt.hashLog, p->mls);
    U32 const dmsIndexDelta = windowLow - dms->highLimit;
    U32 dictMatchIndex = dms->bt.hash[dmsH];
    commonLengthSmaller = commonLengthLarger = 0;
    while (nbCompares-- && (dictMatchIndex > dms->lowLimit)) {
      const U32* const nextPtr = dmsBt + 2*(dictMatchIndex & dmsBtMask);
      size_t matchLength = MIN(commonLengthSmaller, commonLengthLarger);   /* guaranteed minimum nb of common bytes */
      const BYTE* match = dms->base + dictMatchIndex;
      matchLength += ZSTD_count_2segments(ip+matchLength, match+matchLength, iLimit, dms->end, p->prefixStart);
      if (dictMatchIndex+matchLength >= dms->highLimit)
        match = dms->base + dictMatchIndex + dmsIndexDelta;   /* to prepare for next usage of match[matchLength] */

      if (matchLength > *bestLength) {
        matchIndex = dictMatchIndex + dmsIndexDelta;
        DEBUGLOG(8, "found dms match of length %u at distance %u (offCode=%u)",
                (U32)matchLength, current - matchIndex, current - matchIndex + ZSTD_REP_MOVE);
        if (matchLength > matchEndIdx - matchIndex)
            matchEndIdx = matchIndex + (U32)matchLength;
        *bestLength = matchLength;
        cmds[cnum].offset = (current - matchIndex) + ZSTD_REP_MOVE;
        cmds[cnum].length = (U32)matchLength;
        cmds[cnum].type = ZSTD_OCT_match;
        cnum++;
        if ( (matchLength > ZSTD_OPT_MAX_CHUNK)
            | (ip+matchLength == iLimit) /* equal : no way to know if inf or sup */) {
          break;   /* drop, to guarantee consistency (miss a little bit of compression) */
        }
      }

      if (dictMatchIndex <= dms->btLow) { break; }   /* beyond tree size, stop the search */
      if (match[matchLength] < ip[matchLength]) {
        commonLengthSmaller = matchLength;    /* all smaller will now have at least this guaranteed common length */
        dictMatchIndex = nextPtr[1];              /* new matchIndex larger than previous (closer to current) */
      } else {
        /* match is larger than current */
        commonLengthLarger = matchLength;
        dictMatchIndex = nextPtr[0];
      }
    }
  }

  assert(matchEndIdx > current+8);
  ctx->ms->nextToUpdate = matchEndIdx - 8;  /* skip repetitive patterns */
  return cnum;
}

static U32 ZSTD_insertBt1(
                ZSTD_matchState_t* ms,
                const BYTE* const ip, const BYTE* const iend,
                U32 const mls, const int extDict)
{
    const ZSTD_compressionParameters* const cParams = &ms->cParams;
    U32*   const hashTable = ms->hashTable;
    U32    const hashLog = cParams->hashLog;
    size_t const h  = ZSTD_hashPtr(ip, hashLog, mls);
    U32*   const bt = ms->chainTable;
    U32    const btLog  = cParams->chainLog - 1;
    U32    const btMask = (1 << btLog) - 1;
    U32 matchIndex = hashTable[h];
    size_t commonLengthSmaller=0, commonLengthLarger=0;
    const BYTE* const base = ms->window.base;
    const BYTE* const dictBase = ms->window.dictBase;
    const U32 dictLimit = ms->window.dictLimit;
    const BYTE* const dictEnd = dictBase + dictLimit;
    const BYTE* const prefixStart = base + dictLimit;
    const BYTE* match;
    const U32 current = (U32)(ip-base);
    const U32 btLow = btMask >= current ? 0 : current - btMask;
    U32* smallerPtr = bt + 2*(current&btMask);
    U32* largerPtr  = smallerPtr + 1;
    U32 dummy32;   /* to be nullified at the end */
    U32 const windowLow = ms->window.lowLimit;
    U32 matchEndIdx = current+8+1;
    size_t bestLength = 8;
    U32 nbCompares = 1U << cParams->searchLog;
#ifdef ZSTD_C_PREDICT
    U32 predictedSmall = *(bt + 2*((current-1)&btMask) + 0);
    U32 predictedLarge = *(bt + 2*((current-1)&btMask) + 1);
    predictedSmall += (predictedSmall>0);
    predictedLarge += (predictedLarge>0);
#endif /* ZSTD_C_PREDICT */

    DEBUGLOG(8, "ZSTD_insertBt1 (%u)", current);

    assert(ip <= iend-8);   /* required for h calculation */
    hashTable[h] = current;   /* Update Hash Table */

    assert(windowLow > 0);
    while (nbCompares-- && (matchIndex >= windowLow)) {
        U32* const nextPtr = bt + 2*(matchIndex & btMask);
        size_t matchLength = MIN(commonLengthSmaller, commonLengthLarger);   /* guaranteed minimum nb of common bytes */
        assert(matchIndex < current);

#ifdef ZSTD_C_PREDICT   /* note : can create issues when hlog small <= 11 */
        const U32* predictPtr = bt + 2*((matchIndex-1) & btMask);   /* written this way, as bt is a roll buffer */
        if (matchIndex == predictedSmall) {
            /* no need to check length, result known */
            *smallerPtr = matchIndex;
            if (matchIndex <= btLow) { smallerPtr=&dummy32; break; }   /* beyond tree size, stop the search */
            smallerPtr = nextPtr+1;               /* new "smaller" => larger of match */
            matchIndex = nextPtr[1];              /* new matchIndex larger than previous (closer to current) */
            predictedSmall = predictPtr[1] + (predictPtr[1]>0);
            continue;
        }
        if (matchIndex == predictedLarge) {
            *largerPtr = matchIndex;
            if (matchIndex <= btLow) { largerPtr=&dummy32; break; }   /* beyond tree size, stop the search */
            largerPtr = nextPtr;
            matchIndex = nextPtr[0];
            predictedLarge = predictPtr[0] + (predictPtr[0]>0);
            continue;
        }
#endif

        if (!extDict || (matchIndex+matchLength >= dictLimit)) {
            assert(matchIndex+matchLength >= dictLimit);   /* might be wrong if actually extDict */
            match = base + matchIndex;
            matchLength += ZSTD_count(ip+matchLength, match+matchLength, iend);
        } else {
            match = dictBase + matchIndex;
            matchLength += ZSTD_count_2segments(ip+matchLength, match+matchLength, iend, dictEnd, prefixStart);
            if (matchIndex+matchLength >= dictLimit)
                match = base + matchIndex;   /* to prepare for next usage of match[matchLength] */
        }

        if (matchLength > bestLength) {
            bestLength = matchLength;
            if (matchLength > matchEndIdx - matchIndex)
                matchEndIdx = matchIndex + (U32)matchLength;
        }

        if (ip+matchLength == iend) {   /* equal : no way to know if inf or sup */
            break;   /* drop , to guarantee consistency ; miss a bit of compression, but other solutions can corrupt tree */
        }

        if (match[matchLength] < ip[matchLength]) {  /* necessarily within buffer */
            /* match is smaller than current */
            *smallerPtr = matchIndex;             /* update smaller idx */
            commonLengthSmaller = matchLength;    /* all smaller will now have at least this guaranteed common length */
            if (matchIndex <= btLow) { smallerPtr=&dummy32; break; }   /* beyond tree size, stop searching */
            smallerPtr = nextPtr+1;               /* new "candidate" => larger than match, which was smaller than target */
            matchIndex = nextPtr[1];              /* new matchIndex, larger than previous and closer to current */
        } else {
            /* match is larger than current */
            *largerPtr = matchIndex;
            commonLengthLarger = matchLength;
            if (matchIndex <= btLow) { largerPtr=&dummy32; break; }   /* beyond tree size, stop searching */
            largerPtr = nextPtr;
            matchIndex = nextPtr[0];
    }   }

    *smallerPtr = *largerPtr = 0;
    {   U32 positions = 0;
        if (bestLength > 384) positions = MIN(192, (U32)(bestLength - 384));   /* speed optimization */
        assert(matchEndIdx > current + 8);
        return MAX(positions, matchEndIdx - (current + 8));
    }
}

FORCE_INLINE_TEMPLATE
void ZSTD_updateTree_internal(
                ZSTD_matchState_t* ms,
                const BYTE* const ip, const BYTE* const iend,
                const U32 mls, const ZSTD_dictMode_e dictMode)
{
    const BYTE* const base = ms->window.base;
    U32 const target = (U32)(ip - base);
    U32 idx = ms->nextToUpdate;
    DEBUGLOG(6, "ZSTD_updateTree_internal, from %u to %u  (dictMode:%u)",
                idx, target, dictMode);

    while(idx < target) {
        U32 const forward = ZSTD_insertBt1(ms, base+idx, iend, mls, dictMode == ZSTD_extDict);
        assert(idx < (U32)(idx + forward));
        idx += forward;
    }
    assert((size_t)(ip - base) <= (size_t)(U32)(-1));
    assert((size_t)(iend - base) <= (size_t)(U32)(-1));
    ms->nextToUpdate = target;
}

// Zstd matches
// params: State pointer, position
static ZSTD_OptCommands ZSTD_Opt_getCommands_impl(
    ZSTD_OptContext* ctx,
    size_t cur,
    uint8_t const* ip,
    uint8_t const* iend,
    U32 const minMatch /* template */,
    ZSTD_dictMode_e const dictMode /* template */) {
  ZSTD_OptCommand* cmds = ctx->cmds;
  size_t const sufficientLen = ctx->params.sufficientLen;
  size_t bestLength = minMatch - 1;
  cmds += ZSTD_Opt_getLiteralCommand(cmds);
  if (ip < iend - HASH_READ_SIZE) {
    cmds += ZSTD_Opt_getRepCommands(ctx, cmds, cur, ip, iend, &bestLength, minMatch, dictMode);
    if (minMatch == 3 && bestLength < 3) {
      cmds += ZSTD_Opt_getHC3Commands(ctx, cmds, ip, iend, &bestLength, dictMode);
    }
    if (bestLength <= sufficientLen && ip + bestLength < iend) {
      ZSTD_updateTree_internal(ctx->ms, ip, iend, ctx->params.mls, dictMode);
      cmds += ZSTD_Opt_getBtCommands(ctx, cmds, ip, iend, &bestLength, dictMode);
    }
    if (bestLength <= sufficientLen && ip + bestLength < iend) {
      cmds += ZSTD_Opt_getMLRCommands(ctx, cmds, cur, ip, iend, &bestLength, minMatch, dictMode);
    }
  }
  return (ZSTD_OptCommands){.begin = ctx->cmds, .end = cmds};
}

static ZSTD_OptCommands ZSTD_Opt_getCommands(
    ZSTD_OptContext* ctx,
    size_t cur,
    uint8_t const* ip,
    uint8_t const* iend) {
  if (ctx->params.minMatch) {
    switch (ctx->params.dictMode) {
      default:
        assert(0);
        /* fall-through */
      case ZSTD_noDict:
        return ZSTD_Opt_getCommands_impl(ctx, cur, ip, iend, 3, ZSTD_noDict);
      case ZSTD_extDict:
        return ZSTD_Opt_getCommands_impl(ctx, cur, ip, iend, 3, ZSTD_extDict);
      case ZSTD_dictMatchState:
        return ZSTD_Opt_getCommands_impl(ctx, cur, ip, iend, 3, ZSTD_dictMatchState);
    }
  } else {
    assert(ctx->params.minMatch == 4);
    switch (ctx->params.dictMode) {
      default:
        assert(0);
        /* fall-through */
      case ZSTD_noDict:
        return ZSTD_Opt_getCommands_impl(ctx, cur, ip, iend, 4, ZSTD_noDict);
      case ZSTD_extDict:
        return ZSTD_Opt_getCommands_impl(ctx, cur, ip, iend, 4, ZSTD_extDict);
      case ZSTD_dictMatchState:
        return ZSTD_Opt_getCommands_impl(ctx, cur, ip, iend, 4, ZSTD_dictMatchState);
    }
  }
}

////// SECTION:PRICING

#define ZSTD_LITFREQ_ADD    2   /* scaling factor for litFreq, so that frequencies adapt faster to new stats */
#define ZSTD_FREQ_DIV       4   /* log factor when using previous stats to init next stats */
#define ZSTD_PREDEF_THRESHOLD 1024   /* if srcSize < ZSTD_PREDEF_THRESHOLD, symbols' cost is assumed static, directly determined by pre-defined distributions */

#if 0    /* approximation at bit level */
#  define BITCOST_ACCURACY 0
#  define BITCOST_MULTIPLIER (1 << BITCOST_ACCURACY)
#  define WEIGHT(stat)  ((void)opt, ZSTD_bitWeight(stat))
#elif 0  /* fractional bit accuracy */
#  define BITCOST_ACCURACY 8
#  define BITCOST_MULTIPLIER (1 << BITCOST_ACCURACY)
#  define WEIGHT(stat,opt) ((void)opt, ZSTD_fracWeight(stat))
#else    /* opt==approx, ultra==accurate */
#  define BITCOST_ACCURACY 8
#  define BITCOST_MULTIPLIER (1 << BITCOST_ACCURACY)
#  define WEIGHT(stat,opt) (opt ? ZSTD_fracWeight(stat) : ZSTD_bitWeight(stat))
#endif

MEM_STATIC U32 ZSTD_bitWeight(U32 stat)
{
    return (ZSTD_highbit32(stat+1) * BITCOST_MULTIPLIER);
}

MEM_STATIC U32 ZSTD_fracWeight(U32 rawStat)
{
    U32 const stat = rawStat + 1;
    U32 const hb = ZSTD_highbit32(stat);
    U32 const BWeight = hb * BITCOST_MULTIPLIER;
    U32 const FWeight = (stat << BITCOST_ACCURACY) >> hb;
    U32 const weight = BWeight + FWeight;
    assert(hb + BITCOST_ACCURACY < 31);
    return weight;
}

#if (DEBUGLEVEL>=2)
/* debugging function,
 * @return price in bytes as fractional value
 * for debug messages only */
MEM_STATIC double ZSTD_fCost(U32 price)
{
    return (double)price / (BITCOST_MULTIPLIER*8);
}
#endif

static int ZSTD_compressedLiterals(ZSTD_OptContext const* const optPtr)
{
    return optPtr->literalCompressionMode != ZSTD_lcm_uncompressed;
}

static void ZSTD_setBasePrices(ZSTD_OptContext* optPtr, int optLevel)
{
    if (ZSTD_compressedLiterals(optPtr))
        optPtr->litSumBasePrice = WEIGHT(optPtr->litSum, optLevel);
    optPtr->litLengthSumBasePrice = WEIGHT(optPtr->litLengthSum, optLevel);
    optPtr->matchLengthSumBasePrice = WEIGHT(optPtr->matchLengthSum, optLevel);
    optPtr->offCodeSumBasePrice = WEIGHT(optPtr->offCodeSum, optLevel);
}


/* ZSTD_downscaleStat() :
 * reduce all elements in table by a factor 2^(ZSTD_FREQ_DIV+malus)
 * return the resulting sum of elements */
static U32 ZSTD_downscaleStat(unsigned* table, U32 lastEltIndex, int malus)
{
    U32 s, sum=0;
    DEBUGLOG(5, "ZSTD_downscaleStat (nbElts=%u)", (unsigned)lastEltIndex+1);
    assert(ZSTD_FREQ_DIV+malus > 0 && ZSTD_FREQ_DIV+malus < 31);
    for (s=0; s<lastEltIndex+1; s++) {
        table[s] = 1 + (table[s] >> (ZSTD_FREQ_DIV+malus));
        sum += table[s];
    }
    return sum;
}

/* ZSTD_rescaleFreqs() :
 * if first block (detected by optPtr->litLengthSum == 0) : init statistics
 *    take hints from dictionary if there is one
 *    or init from zero, using src for literals stats, or flat 1 for match symbols
 * otherwise downscale existing stats, to be used as seed for next block.
 */
static void
ZSTD_rescaleFreqs(ZSTD_OptContext* const optPtr,
            const BYTE* const src, size_t const srcSize,
                  int const optLevel)
{
    int const compressedLiterals = ZSTD_compressedLiterals(optPtr);
    DEBUGLOG(5, "ZSTD_rescaleFreqs (srcSize=%u)", (unsigned)srcSize);
    optPtr->priceType = zop_dynamic;

    if (optPtr->litLengthSum == 0) {  /* first block : init */
        if (srcSize <= ZSTD_PREDEF_THRESHOLD) {  /* heuristic */
            DEBUGLOG(5, "(srcSize <= ZSTD_PREDEF_THRESHOLD) => zop_predef");
            optPtr->priceType = zop_predef;
        }

        assert(optPtr->symbolCosts != NULL);
        if (optPtr->symbolCosts->huf.repeatMode == HUF_repeat_valid) {
            /* huffman table presumed generated by dictionary */
            optPtr->priceType = zop_dynamic;

            if (compressedLiterals) {
                unsigned lit;
                assert(optPtr->litFreq != NULL);
                optPtr->litSum = 0;
                for (lit=0; lit<=MaxLit; lit++) {
                    U32 const scaleLog = 11;   /* scale to 2K */
                    U32 const bitCost = HUF_getNbBits(optPtr->symbolCosts->huf.CTable, lit);
                    assert(bitCost <= scaleLog);
                    optPtr->litFreq[lit] = bitCost ? 1 << (scaleLog-bitCost) : 1 /*minimum to calculate cost*/;
                    optPtr->litSum += optPtr->litFreq[lit];
            }   }

            {   unsigned ll;
                FSE_CState_t llstate;
                FSE_initCState(&llstate, optPtr->symbolCosts->fse.litlengthCTable);
                optPtr->litLengthSum = 0;
                for (ll=0; ll<=MaxLL; ll++) {
                    U32 const scaleLog = 10;   /* scale to 1K */
                    U32 const bitCost = FSE_getMaxNbBits(llstate.symbolTT, ll);
                    assert(bitCost < scaleLog);
                    optPtr->litLengthFreq[ll] = bitCost ? 1 << (scaleLog-bitCost) : 1 /*minimum to calculate cost*/;
                    optPtr->litLengthSum += optPtr->litLengthFreq[ll];
            }   }

            {   unsigned ml;
                FSE_CState_t mlstate;
                FSE_initCState(&mlstate, optPtr->symbolCosts->fse.matchlengthCTable);
                optPtr->matchLengthSum = 0;
                for (ml=0; ml<=MaxML; ml++) {
                    U32 const scaleLog = 10;
                    U32 const bitCost = FSE_getMaxNbBits(mlstate.symbolTT, ml);
                    assert(bitCost < scaleLog);
                    optPtr->matchLengthFreq[ml] = bitCost ? 1 << (scaleLog-bitCost) : 1 /*minimum to calculate cost*/;
                    optPtr->matchLengthSum += optPtr->matchLengthFreq[ml];
            }   }

            {   unsigned of;
                FSE_CState_t ofstate;
                FSE_initCState(&ofstate, optPtr->symbolCosts->fse.offcodeCTable);
                optPtr->offCodeSum = 0;
                for (of=0; of<=MaxOff; of++) {
                    U32 const scaleLog = 10;
                    U32 const bitCost = FSE_getMaxNbBits(ofstate.symbolTT, of);
                    assert(bitCost < scaleLog);
                    optPtr->offCodeFreq[of] = bitCost ? 1 << (scaleLog-bitCost) : 1 /*minimum to calculate cost*/;
                    optPtr->offCodeSum += optPtr->offCodeFreq[of];
            }   }

        } else {  /* not a dictionary */

            assert(optPtr->litFreq != NULL);
            if (compressedLiterals) {
                unsigned lit = MaxLit;
                HIST_count_simple(optPtr->litFreq, &lit, src, srcSize);   /* use raw first block to init statistics */
                optPtr->litSum = ZSTD_downscaleStat(optPtr->litFreq, MaxLit, 1);
            }

            {   unsigned ll;
                for (ll=0; ll<=MaxLL; ll++)
                    optPtr->litLengthFreq[ll] = 1;
            }
            optPtr->litLengthSum = MaxLL+1;

            {   unsigned ml;
                for (ml=0; ml<=MaxML; ml++)
                    optPtr->matchLengthFreq[ml] = 1;
            }
            optPtr->matchLengthSum = MaxML+1;

            {   unsigned of;
                for (of=0; of<=MaxOff; of++)
                    optPtr->offCodeFreq[of] = 1;
            }
            optPtr->offCodeSum = MaxOff+1;

        }

    } else {   /* new block : re-use previous statistics, scaled down */

        if (compressedLiterals)
            optPtr->litSum = ZSTD_downscaleStat(optPtr->litFreq, MaxLit, 1);
        optPtr->litLengthSum = ZSTD_downscaleStat(optPtr->litLengthFreq, MaxLL, 0);
        optPtr->matchLengthSum = ZSTD_downscaleStat(optPtr->matchLengthFreq, MaxML, 0);
        optPtr->offCodeSum = ZSTD_downscaleStat(optPtr->offCodeFreq, MaxOff, 0);
    }

    ZSTD_setBasePrices(optPtr, optLevel);
}

/* ZSTD_rawLiteralsCost() :
 * price of literals (only) in specified segment (which length can be 0).
 * does not include price of literalLength symbol */
static U32 ZSTD_rawLiteralsCost(const BYTE* const literals, U32 const litLength,
                                const ZSTD_OptContext* const optPtr,
                                int optLevel)
{
    if (litLength == 0) return 0;

    if (!ZSTD_compressedLiterals(optPtr))
        return (litLength << 3) * BITCOST_MULTIPLIER;  /* Uncompressed - 8 bytes per literal. */

    if (optPtr->priceType == zop_predef)
        return (litLength*6) * BITCOST_MULTIPLIER;  /* 6 bit per literal - no statistic used */

    /* dynamic statistics */
    {   U32 price = litLength * optPtr->litSumBasePrice;
        U32 u;
        for (u=0; u < litLength; u++) {
            assert(WEIGHT(optPtr->litFreq[literals[u]], optLevel) <= optPtr->litSumBasePrice);   /* literal cost should never be negative */
            price -= WEIGHT(optPtr->litFreq[literals[u]], optLevel);
        }
        return price;
    }
}

static U32 ZSTD_rawLiteralCost(BYTE const literal, const ZSTD_OptContext* const optPtr, int optLevel)
{
    if (!ZSTD_compressedLiterals(optPtr))
        return 8 * BITCOST_MULTIPLIER;  /* Uncompressed - 8 bytes per literal. */

    if (optPtr->priceType == zop_predef)
        return 6 * BITCOST_MULTIPLIER;  /* 6 bit per literal - no statistic used */

    return optPtr->litSumBasePrice - WEIGHT(optPtr->litFreq[literal], optLevel);
}

/* ZSTD_litLengthPrice() :
 * cost of literalLength symbol */
static U32 ZSTD_litLengthPrice(U32 const litLength, const ZSTD_OptContext* const optPtr, int optLevel)
{
    if (optPtr->priceType == zop_predef) return WEIGHT(litLength, optLevel);

    /* dynamic statistics */
    {   U32 const llCode = ZSTD_LLcode(litLength);
        return (LL_bits[llCode] * BITCOST_MULTIPLIER)
             + optPtr->litLengthSumBasePrice
             - WEIGHT(optPtr->litLengthFreq[llCode], optLevel);
    }
}

/* ZSTD_litLengthContribution() :
 * @return ( cost(litlength) - cost(0) )
 * this value can then be added to rawLiteralsCost()
 * to provide a cost which is directly comparable to a match ending at same position */
// static int ZSTD_litLengthContribution(U32 const litLength, const ZSTD_OptContext* const optPtr, int optLevel)
// {
//     if (optPtr->priceType >= zop_predef) return (int)WEIGHT(litLength, optLevel);

//     /* dynamic statistics */
//     {   U32 const llCode = ZSTD_LLcode(litLength);
//         int const contribution = (int)(LL_bits[llCode] * BITCOST_MULTIPLIER)
//                                + (int)WEIGHT(optPtr->litLengthFreq[0], optLevel)   /* note: log2litLengthSum cancel out */
//                                - (int)WEIGHT(optPtr->litLengthFreq[llCode], optLevel);
// #if 1
//         return contribution;
// #else
//         return MAX(0, contribution); /* sometimes better, sometimes not ... */
// #endif
//     }
// }

/* ZSTD_getMatchPrice() :
 * Provides the cost of the match part (offset + matchLength) of a sequence
 * Must be combined with ZSTD_fullLiteralsCost() to get the full cost of a sequence.
 * optLevel: when <2, favors small offset for decompression speed (improved cache efficiency) */
FORCE_INLINE_TEMPLATE U32
ZSTD_getMatchPrice(U32 const offset,
                   U32 const matchLength,
             const ZSTD_OptContext* const optPtr,
                   int const optLevel)
{
    U32 price;
    U32 const offCode = ZSTD_highbit32(offset+1);
    U32 const mlBase = matchLength - MINMATCH;
    assert(matchLength >= MINMATCH);

    if (optPtr->priceType == zop_predef)  /* fixed scheme, do not use statistics */
        return WEIGHT(mlBase, optLevel) + ((16 + offCode) * BITCOST_MULTIPLIER);

    /* dynamic statistics */
    price = (offCode * BITCOST_MULTIPLIER) + (optPtr->offCodeSumBasePrice - WEIGHT(optPtr->offCodeFreq[offCode], optLevel));
    if ((optLevel<2) /*static*/ && offCode >= 20)
        price += (offCode-19)*2 * BITCOST_MULTIPLIER; /* handicap for long distance offsets, favor decompression speed */

    /* match Length */
    {   U32 const mlCode = ZSTD_MLcode(mlBase);
        price += (ML_bits[mlCode] * BITCOST_MULTIPLIER) + (optPtr->matchLengthSumBasePrice - WEIGHT(optPtr->matchLengthFreq[mlCode], optLevel));
    }

    price += BITCOST_MULTIPLIER / 5;   /* heuristic : make matches a bit more costly to favor less sequences -> faster decompression speed */

    DEBUGLOG(8, "ZSTD_getMatchPrice(ml:%u) = %u", matchLength, price);
    return price;
}

/* ZSTD_updateStats() :
 * assumption : literals + litLengtn <= iend */
static void ZSTD_updateStats(ZSTD_OptContext* const optPtr,
                             U32 litLength, const BYTE* literals,
                             U32 offsetCode, U32 matchLength)
{
    /* literals */
    if (ZSTD_compressedLiterals(optPtr)) {
        U32 u;
        for (u=0; u < litLength; u++)
            optPtr->litFreq[literals[u]] += ZSTD_LITFREQ_ADD;
        optPtr->litSum += litLength*ZSTD_LITFREQ_ADD;
    }

    /* literal Length */
    {   U32 const llCode = ZSTD_LLcode(litLength);
        optPtr->litLengthFreq[llCode]++;
        optPtr->litLengthSum++;
    }

    /* match offset code (0-2=>repCode; 3+=>offset+2) */
    {   U32 const offCode = ZSTD_highbit32(offsetCode+1);
        assert(offCode <= MaxOff);
        optPtr->offCodeFreq[offCode]++;
        optPtr->offCodeSum++;
    }

    /* match Length */
    {   U32 const mlBase = matchLength - MINMATCH;
        U32 const mlCode = ZSTD_MLcode(mlBase);
        optPtr->matchLengthFreq[mlCode]++;
        optPtr->matchLengthSum++;
    }
}

FORCE_INLINE_TEMPLATE uint32_t ZSTD_Opt_getCommandPrice_impl(
    ZSTD_OptContext const* ctx,
    BYTE const* ip,
    size_t cur,
    ZSTD_OptCommand const* cmd,
    size_t len,
    int const optLevel /* template */) {
  uint32_t const prevPrice = ctx->opt[cur].price;
  if (cmd->type == ZSTD_OCT_literals) {
    uint32_t const prevLitlen = ctx->opt[cur].cmd.type == ZSTD_OCT_literals ? ctx->opt[cur].cmd.length : 0;
    uint32_t const litCost = ZSTD_rawLiteralCost(ip[cur], ctx, optLevel);
    assert(len == cmd->length);
    assert(len == 1);
    return prevPrice + litCost
                     + ZSTD_litLengthPrice(prevLitlen + 1, ctx, optLevel)
                     - ZSTD_litLengthPrice(prevLitlen, ctx, optLevel);
  } else if (cmd->type == ZSTD_OCT_match) {
    assert(len <= cmd->length);
    return prevPrice + ZSTD_litLengthPrice(0, ctx, optLevel) +
                     + ZSTD_getMatchPrice(cmd->offset, len, ctx, optLevel);
  } else {
    uint32_t const mlen = cmd->length - cmd->lits - cmd->repLen;
    uint32_t const repLen = cmd->repLen + (cmd->length - len);
    assert(len <= cmd->length);
    assert(cmd->type == ZSTD_OCT_mlr);
    assert(cmd->length >= cmd->lits + cmd->repLen + MINMATCH);
    assert(len >= mlen + cmd->lits + MINMATCH);
    return prevPrice + ZSTD_getMatchPrice(cmd->offset, mlen, ctx, optLevel)
                     + ZSTD_litLengthPrice(cmd->lits, ctx, optLevel)
                     + ZSTD_rawLiteralsCost(ip + cur + mlen, cmd->lits, ctx, optLevel)
                     + ZSTD_getMatchPrice(cmd->rep, repLen, ctx, optLevel)
                     + ZSTD_litLengthPrice(0, ctx, optLevel);
  }
}
// Params: State ptr, command, path
static uint32_t ZSTD_Opt_getCommandPrice(
    ZSTD_OptContext const* ctx,
    BYTE const* ip,
    size_t cur,
    ZSTD_OptCommand const* cmd,
    size_t len) {
  if (ctx->params.optLevel == 0) {
    return ZSTD_Opt_getCommandPrice_impl(ctx, ip, cur, cmd, len, 0);
  }
  assert(ctx->params.optLevel == 2);
  return ZSTD_Opt_getCommandPrice_impl(ctx, ip, cur, cmd, len, 2);
}

// Should we stop early for whatever reason
static int ZSTD_Opt_break(ZSTD_OptContext const* ctx) {
  (void)ctx;
  return 0;
}

// Should we force-select this match for speed?
static int ZSTD_Opt_immediateEncoding(
    ZSTD_OptContext const* ctx,
    ZSTD_OptCommand const* cmd) {
  return cmd->type == ZSTD_OCT_match && cmd->length > ctx->params.sufficientLen;
}

// Returns the shortest allowed length
static size_t ZSTD_OptCommand_shortestLen(
    ZSTD_OptContext const* ctx,
    ZSTD_OptCommand const* cmd,
    size_t shortest) {
  (void)ctx;
  if (cmd->type == ZSTD_OCT_match) {
    return MAX(shortest, MINMATCH);
  }
  if (cmd->type == ZSTD_OCT_mlr) {
    // Must include match + lits + MINMATCH
    assert(cmd->repLen >= MINMATCH);
    return MAX(shortest, cmd->length - cmd->repLen + MINMATCH);
  }
  return shortest;
}

static void ZSTD_Opt_fillBackward(ZSTD_OptContext* ctx, size_t cur) {
  if (ctx->opt[cur].cmd.type == ZSTD_OCT_literals) {
    if (cur != 0) {
      size_t const backDist = ctx->opt[cur].cmd.length;
      size_t const prev = cur > backDist ? cur - backDist : 0;
      memcpy(&ctx->opt[cur].reps, &ctx->opt[prev].reps, sizeof(ctx->opt[cur].reps));
    }
  } else {
    uint32_t const prev = cur - ctx->opt[cur].cmd.length;
    U32 const ll0 = ctx->opt[prev].cmd.type != ZSTD_OCT_literals || (prev == 0 && ctx->opt[prev].cmd.length == 0);
    assert(ctx->opt[cur].cmd.type == ZSTD_OCT_match || ctx->opt[cur].cmd.type == ZSTD_OCT_mlr);
    assert(cur >= ctx->opt[cur].cmd.length);
    ctx->opt[cur].reps = ZSTD_updateRep(ctx->opt[prev].reps.rep, ctx->opt[cur].cmd.offset, ll0);
    if (ctx->opt[cur].cmd.type == ZSTD_OCT_mlr) {
      ctx->opt[cur].reps = ZSTD_updateRep(ctx->opt[cur].reps.rep, ctx->opt[cur].cmd.rep, ctx->opt[cur].cmd.lits == 0);
    }
  }
}

static void ZSTD_Opt_fillForward(
    ZSTD_OptContext* ctx,
    size_t pos,
    ZSTD_OptCommand const* cmd,
    size_t len) {
  ctx->opt[pos].cmd = *cmd;
  if (cmd->type == ZSTD_OCT_literals) {
    U32 const plen = ctx->opt[pos-1].cmd.type == ZSTD_OCT_literals ? ctx->opt[pos-1].cmd.length : 0;
    assert(len == 1);
    assert(pos > 0);
    ctx->opt[pos].cmd.length = 1 + plen;
  } else {
    if (cmd->type == ZSTD_OCT_mlr) {
      ctx->opt[pos].cmd.repLen -= ctx->opt[pos].cmd.length - len;
      assert(ctx->opt[pos].cmd.repLen >= MINMATCH);
    }
    ctx->opt[pos].cmd.length = len;
  }
}

// Fill the first position (including the initial price + length)
static void ZSTD_Opt_fill0(ZSTD_OptContext const* ctx, ZSTD_OptState* opt0, uint8_t const* istart, uint8_t const* ip) {
  uint32_t const litlen = ip - istart;
  opt0[0].cmd.type = ZSTD_OCT_literals;
  opt0[0].cmd.length = litlen;
  opt0[0].price = ZSTD_litLengthPrice(litlen, ctx, ctx->params.optLevel);
  memcpy(&opt0[0].reps, &ctx->initialReps, sizeof(ctx->initialReps));
}

typedef struct {
  size_t lastPos; // May be past the end!
  ZSTD_OptCommand lastCmd;
} ZSTD_OptParseResult;

// istart: place where current segment starts (may be ending literals)
// ip: start parsing here
// iend: end of source
// istart <= ip < iend
static ZSTD_OptParseResult
parse(ZSTD_OptContext* const ctx, uint8_t const* const istart, uint8_t const* const ip, uint8_t const* const iend) {
  size_t const curEnd         = MIN(iend - ip, ZSTD_OPT_MAX_CHUNK);
  ZSTD_OptState* const opt      = ctx->opt;
  size_t cur                  = 0;
  size_t endPos               = 0;

  ZSTD_Opt_fill0(ctx, &opt[0], istart, ip);
  DEBUGLOG(NT, "parse");

  for (;; ++cur) {
    DEBUGLOG(NT, "cur = %zu", cur);
    // We've filled our buffer
    if (endPos == curEnd)
      break;
    // We've hit a choke point in the parse
    // TODO: Handle llen=1, probably don't need to go 1 lit at a time...
    //       But maybe it is only the beginning...
    if (0 < cur && cur == endPos)
      break;
    DEBUGLOG(NT, "%zu %zu", cur, endPos);
    assert(cur <= endPos);
    assert(endPos < curEnd);

    // Update the state upon arrival
    ZSTD_Opt_fillBackward(ctx, cur);

    // Check if we need to stop for whatever reason.
    // E.g. out of match space.
    if (ZSTD_Opt_break(ctx)) {
      break;
    }
    // Get the commands starting at this position in length sorted order
    ZSTD_OptCommands cmds = ZSTD_Opt_getCommands(ctx, cur, ip + cur, iend);
    size_t const numCmds = (size_t)(cmds.end - cmds.begin);
    if (numCmds == 0) {
      continue;
    }
    {
      size_t const maxPos = cur + cmds.end[-1].length;
      assert(maxPos > 0);
      // Check if the longest command is longer than allowed, or
      // is designated for immediate encoding.
      if (maxPos >= ZSTD_OPT_MAX_CHUNK ||
          ZSTD_Opt_immediateEncoding(ctx, cmds.end - 1)) {
        return (ZSTD_OptParseResult){
          .lastPos = maxPos,
          .lastCmd = cmds.end[-1],
        };
      }
      // Fill in any empty positions.
      while (endPos < maxPos) {
        opt[++endPos].price = ZSTD_OPT_MAX_PRICE;
      }
    }
    size_t firstLen = 1;
    // Loop over all the commands
    for (size_t c = 0; c < numCmds; ++c) {
      ZSTD_OptCommand const* cmd = &cmds.begin[c];
      size_t const lastLen     = cmd->length;
      // DEBUGLOG(NT, "fl = %zu ll = %u", firstLen, lastLen);
      assert(firstLen <= lastLen);
      firstLen                 = ZSTD_OptCommand_shortestLen(ctx, cmd, firstLen);
      assert(firstLen <= lastLen);
      assert(firstLen > 0);
      // Try all the allowed lengths.
      for (size_t len = lastLen; len >= firstLen; --len) {
        size_t const pos     = cur + len;
        uint32_t const price = ZSTD_Opt_getCommandPrice(ctx, ip, cur, cmd, len);
        assert(pos <= endPos);
        if (price < ctx->opt[pos].price) {
          DEBUGLOG(NT, "set %zu = %u", pos, price);
          ctx->opt[pos].price = price;
          ZSTD_Opt_fillForward(ctx, pos, cmd, len);
        }
      }
      firstLen = lastLen + 1;
    }
  }
  return (ZSTD_OptParseResult){
    .lastPos = endPos,
    .lastCmd = opt[endPos].cmd,
  };
}

typedef struct {
  ZSTD_OptState const* begin;
  ZSTD_OptState const* end;
} ZSTD_OptStateRange;

static ZSTD_OptStateRange reverse(ZSTD_OptState* opt, ZSTD_OptParseResult const* parse) {
  size_t seqPos = parse->lastPos;
  if (seqPos > parse->lastCmd.length) {
    seqPos -= parse->lastCmd.length;
  } else {
    seqPos = 0;
  }
  size_t const storeEnd = seqPos + 1;
  size_t storeStart = storeEnd;
  opt[storeEnd].cmd = parse->lastCmd;
  while (seqPos > 0) {
    size_t const backDist = opt[seqPos].cmd.length;
    DEBUGLOG(NT, "pos=%zu type = %d price = %u", seqPos, opt[seqPos].cmd.type, opt[seqPos].price);
    assert(backDist >= 1);
    assert(seqPos < storeStart);
    opt[--storeStart] = opt[seqPos];
    seqPos = (seqPos > backDist) ? seqPos - backDist : 0;
    assert(seqPos < storeStart);
  }
  return (ZSTD_OptStateRange){.begin = opt + storeStart, .end = opt + storeEnd + 1};
}

static void emit(
        ZSTD_OptContext* ctx,
        seqStore_t* seqStore,
        ZSTD_OptStateRange cmds,
        uint8_t const** anchor,
        uint8_t const** ip,
        uint8_t const* iend)
{
  U32 llen = *ip - *anchor;
  for (ZSTD_OptState const* pos = cmds.begin; pos != cmds.end; ++pos) {
    if (pos->cmd.type == ZSTD_OCT_literals) {
      assert(pos == cmds.begin || pos[-1].cmd.type != ZSTD_OCT_literals);
      DEBUGLOG(NT, "lits %u + %u = %u", llen, pos->cmd.length, llen + pos->cmd.length);
      if (pos == cmds.begin) {
        llen = 0;
      }
      llen += pos->cmd.length;
    } else if (pos->cmd.type == ZSTD_OCT_match) {
      U32 const mlen = pos->cmd.length;
      U32 const offset = pos->cmd.offset;
      U32 const advance = llen + mlen;
      DEBUGLOG(NT, "match llen=%u mlen=%u", llen, mlen);

      if (DEBUGLEVEL >= 1) {
        if (offset < ZSTD_REP_NUM) {
          U32 const* rep = ctx->initialReps.rep;
          U32 const repCode = offset + (llen == 0);
          U32 const repOffset = (repCode==ZSTD_REP_NUM) ? (rep[0] - 1) : rep[repCode];
          assert(!memcmp(*anchor + llen, *anchor + llen - repOffset, mlen));
        } else {
          assert(!memcmp(*anchor + llen, *anchor + llen - (offset - ZSTD_REP_MOVE), mlen));
        }
      }

      assert(mlen >= MINMATCH);
      ctx->initialReps = ZSTD_updateRep(ctx->initialReps.rep, offset, llen == 0);

      assert(*anchor + advance <= iend);
      ZSTD_updateStats(ctx, llen, *anchor, offset, mlen);
      ZSTD_storeSeq(seqStore, llen, *anchor, iend, offset, mlen-MINMATCH);
      *anchor += advance;
      *ip = *anchor;
      llen = 0;
    } else if (pos->cmd.type == ZSTD_OCT_mlr) {
      assert(pos->cmd.repLen >= MINMATCH);
      {
        U32 const mlen = pos->cmd.length - pos->cmd.repLen - pos->cmd.lits;
        U32 const offset = pos->cmd.offset;
        U32 const advance = llen + mlen;

        assert(mlen >= MINMATCH);
        if (DEBUGLEVEL >= 1) {
          if (offset < ZSTD_REP_NUM) {
            U32 const* rep = ctx->initialReps.rep;
            U32 const repCode = offset + (llen == 0);
            U32 const repOffset = (repCode==ZSTD_REP_NUM) ? (rep[0] - 1) : rep[repCode];
            assert(!memcmp(*anchor + llen, *anchor + llen - repOffset, mlen));
          } else {
            assert(!memcmp(*anchor + llen, *anchor + llen - (offset - ZSTD_REP_MOVE), mlen));
          }
        }
        
        ctx->initialReps = ZSTD_updateRep(ctx->initialReps.rep, offset, llen == 0);
        
        assert(*anchor + advance <= iend);
        ZSTD_updateStats(ctx, llen, *anchor, offset, mlen);
        ZSTD_storeSeq(seqStore, llen, *anchor, iend, offset, mlen-MINMATCH);
        *anchor += advance;
        *ip = *anchor;
      }
      llen = pos->cmd.lits;
      {
        U32 const mlen = pos->cmd.repLen;
        U32 const offset = pos->cmd.rep;
        U32 const advance = llen + mlen;

        assert(offset < ZSTD_REP_NUM);
        if (DEBUGLEVEL >= 1) {
          U32 const* rep = ctx->initialReps.rep;
          U32 const repCode = offset + (llen == 0);
          U32 const repOffset = (repCode==ZSTD_REP_NUM) ? (rep[0] - 1) : rep[repCode];
          assert(!memcmp(*anchor + llen, *anchor + llen - repOffset, mlen));
        }
        
        assert(mlen >= MINMATCH);
        ctx->initialReps = ZSTD_updateRep(ctx->initialReps.rep, offset, llen == 0);
        
        assert(*anchor + advance <= iend);
        ZSTD_updateStats(ctx, llen, *anchor, offset, mlen);
        ZSTD_storeSeq(seqStore, llen, *anchor, iend, offset, mlen-MINMATCH);
        *anchor += advance;
      }
      *ip = *anchor;
      llen = 0;
    } else {
      assert(0);
    }
  }
  *ip = *anchor + llen;
}

static size_t
ZSTD_compressBlock_opt2_generic(ZSTD_matchState_t* ms,
                                seqStore_t* seqStore,
                                U32 rep[ZSTD_REP_NUM],
                          const void* src, size_t srcSize,
                          const int optLevel,
                          const ZSTD_dictMode_e dictMode)
{
  uint8_t const* ip = (uint8_t const*)src;
  uint8_t const* iend = ip + srcSize;
  uint8_t const* anchor = ip;
  ZSTD_OptContext* const ctx = ms->opt2;
  fillOptParams(&ctx->params, ms, dictMode, optLevel);
  ZSTD_rescaleFreqs(ctx, src, srcSize, optLevel);
  memcpy(&ctx->initialReps, rep, sizeof(ctx->initialReps));
  ctx->nextToUpdate3 = ms->nextToUpdate;

  while (ip < iend) {
    ZSTD_OptParseResult const res = parse(ctx, anchor, ip, iend);
    ZSTD_OptStateRange const cmds = reverse(ctx->opt, &res);
    emit(ctx, seqStore, cmds, &anchor, &ip, iend);
    ZSTD_setBasePrices(ctx, optLevel);
  }
  memcpy(rep, &ctx->initialReps, sizeof(ctx->initialReps));
  return iend - anchor;
}


size_t ZSTD_compressBlock_btopt2(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        const void* src, size_t srcSize)
{
    DEBUGLOG(5, "ZSTD_compressBlock_btopt");
    return ZSTD_compressBlock_opt2_generic(ms, seqStore, rep, src, srcSize, 2 /*optLevel*/, ZSTD_noDict);
}

size_t ZSTD_compressBlock_btopt2_dictMatchState(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        const void* src, size_t srcSize)
{
    return ZSTD_compressBlock_opt2_generic(ms, seqStore, rep, src, srcSize, 2 /*optLevel*/, ZSTD_dictMatchState);
}

size_t ZSTD_compressBlock_btopt2_extDict(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        const void* src, size_t srcSize)
{
    return ZSTD_compressBlock_opt2_generic(ms, seqStore, rep, src, srcSize, 2 /*optLevel*/, ZSTD_extDict);
}

// fillOptParams()
// rescaleFreqs
// while ip < iend
//   parse
//   reverse
//   emit + updateFreqs
//   rescale/update/whatever
// last literals

#else
int dummy;
#endif
