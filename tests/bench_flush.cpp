#define HUF_STATIC_LINKING_ONLY
#include <zstd_internal.h>
#include <huf.h>
#include <fse.h>
#include <datagen.h>
#include <benchmark/benchmark.h>

static size_t throwOnError(size_t code) {
  if (ZSTD_isError(code)) {
    throw std::logic_error{"wtf"};
  }
  return code;
}

static void BM_CompressFlush(benchmark::State& state) {
  std::string data;
  data.resize(1ul << 20);
  std::string tmp;
  tmp.resize(ZSTD_CStreamOutSize());
  RDG_genBuffer(&data[0], data.size(), 0.5, 0, 0x12345678);
  ZSTD_CStream* stream = ZSTD_createCStream();
  throwOnError(ZSTD_initCStream(stream, 1));
  char const* pData = data.data();
  char const* pEnd  = data.data() + data.size();
  size_t dataProcessed = 0;
  while (state.KeepRunning()) {
    if (pData == pEnd) {
      pData = data.data();
    }
    ZSTD_inBuffer in = {pData, (size_t)std::min(state.range(0), (int)(pEnd - pData)), 0};
    while (in.pos != in.size) {
      ZSTD_outBuffer out = {&tmp[0], tmp.size(), 0};
      throwOnError(ZSTD_compressStream(stream, &out, &in));
      benchmark::DoNotOptimize(out.dst);
    }
    dataProcessed += in.pos;
    pData += in.pos;
    ZSTD_outBuffer out = {&tmp[0], tmp.size(), 0};
    throwOnError(ZSTD_flushStream(stream, &out));
    benchmark::DoNotOptimize(out.dst);
  }
  state.SetBytesProcessed(dataProcessed);
  ZSTD_freeCStream(stream);
}

// static int HUF_validateCTableBaseline(HUF_CElt* CTable, const unsigned* count, unsigned maxSymbolValue) {
//   // unsigned s;
//   // for (s = 0; s <= maxSymbolValue; ++s) {
//   //   if (CTable[s].nbBits == 0 && count[s] != 0) return 0;
//   // }
//   // return 1;
//   int s;
//   for (s = 0; s <= (int)maxSymbolValue; ++s) {
//     if (CTable[s].nbBits == 0 && count[s] != 0) return 0;
//   }
//   return 1;
// }
//
// static void BM_ValidateCTableBaseline(benchmark::State& state) {
//   std::string data;
//   data.resize(1ul << 10);
//   RDG_genBuffer(&data[0], data.size(), 0.5, 0, 0x12345678);
//   U32 count[HUF_SYMBOLVALUE_MAX+1];
//   HUF_CREATE_STATIC_CTABLE(CTable, HUF_SYMBOLVALUE_MAX);
//   unsigned maxSymbolValue = 255;
//   unsigned huffLog = 11;
//   throwOnError(FSE_count(count, &maxSymbolValue, data.data(), data.size()));
//   huffLog = HUF_optimalTableLog(huffLog, data.size(), maxSymbolValue);
//   huffLog = throwOnError(HUF_buildCTable(CTable, count, maxSymbolValue, huffLog));
//
//   while (state.KeepRunning()) {
//     benchmark::ClobberMemory();
//     int const valid = HUF_validateCTableBaseline(CTable, count, maxSymbolValue);
//     benchmark::DoNotOptimize(valid);
//   }
// }
//
// static int HUF_validateCTable(const HUF_CElt* CTable, const unsigned* count, unsigned maxSymbolValue) {
//   int bad = 0;
//   int s;
//   for (s = 0; s <= (int)maxSymbolValue; ++s) {
//     bad |= (count[s] != 0) & (CTable[s].nbBits == 0);
//   }
//   return bad;
// }
//
// static void BM_ValidateCTable(benchmark::State& state) {
//   std::string data;
//   data.resize(1ul << 10);
//   RDG_genBuffer(&data[0], data.size(), 0.5, 0, 0x12345678);
//   U32 count[HUF_SYMBOLVALUE_MAX+1];
//   HUF_CREATE_STATIC_CTABLE(CTable, HUF_SYMBOLVALUE_MAX);
//   unsigned maxSymbolValue = 255;
//   unsigned huffLog = 11;
//   throwOnError(FSE_count(count, &maxSymbolValue, data.data(), data.size()));
//   huffLog = HUF_optimalTableLog(huffLog, data.size(), maxSymbolValue);
//   huffLog = throwOnError(HUF_buildCTable(CTable, count, maxSymbolValue, huffLog));
//
//   while (state.KeepRunning()) {
//     benchmark::ClobberMemory();
//     int const valid = HUF_validateCTable(CTable, count, maxSymbolValue);
//     benchmark::DoNotOptimize(valid);
//   }
// }
//
// static size_t HUF_estimateCompressedSizeBaseline(HUF_CElt* CTable, const unsigned* count, unsigned maxSymbolValue)
// {
//     size_t nbBits = 0;
//     unsigned s;
//     for (s = 0; s <= maxSymbolValue; ++s) {
//         if (CTable[s].nbBits == 0 && count[s] != 0) {
//           return ERROR(GENERIC);
//         }
//         nbBits += CTable[s].nbBits * count[s];
//     }
//     return nbBits >> 3;
// }
//
// static void BM_estimateBaseline(benchmark::State& state) {
//   std::string data;
//   data.resize(1ul << 10);
//   RDG_genBuffer(&data[0], data.size(), 0.5, 0, 0x12345678);
//   U32 count[HUF_SYMBOLVALUE_MAX+1];
//   HUF_CREATE_STATIC_CTABLE(CTable, HUF_SYMBOLVALUE_MAX);
//   unsigned maxSymbolValue = 255;
//   unsigned huffLog = 11;
//   throwOnError(FSE_count(count, &maxSymbolValue, data.data(), data.size()));
//   huffLog = HUF_optimalTableLog(huffLog, data.size(), maxSymbolValue);
//   huffLog = throwOnError(HUF_buildCTable(CTable, count, maxSymbolValue, huffLog));
//
//   while (state.KeepRunning()) {
//     benchmark::ClobberMemory();
//     size_t const size = HUF_estimateCompressedSizeBaseline(CTable, count, maxSymbolValue);
//     benchmark::DoNotOptimize(size);
//   }
// }
// static size_t HUF_estimateCompressedSize(HUF_CElt* CTable, const unsigned* count, unsigned maxSymbolValue)
// {
//     size_t nbBits = 0;
//     int s;
//     for (s = 0; s <= (int)maxSymbolValue; ++s) {
//         nbBits += CTable[s].nbBits * count[s];
//     }
//     return nbBits >> 3;
// }
//
// static void BM_estimate(benchmark::State& state) {
//   std::string data;
//   data.resize(1ul << 10);
//   RDG_genBuffer(&data[0], data.size(), 0.5, 0, 0x12345678);
//   U32 count[HUF_SYMBOLVALUE_MAX+1];
//   HUF_CREATE_STATIC_CTABLE(CTable, HUF_SYMBOLVALUE_MAX);
//   unsigned maxSymbolValue = 255;
//   unsigned huffLog = 11;
//   throwOnError(FSE_count(count, &maxSymbolValue, data.data(), data.size()));
//   huffLog = HUF_optimalTableLog(huffLog, data.size(), maxSymbolValue);
//   huffLog = throwOnError(HUF_buildCTable(CTable, count, maxSymbolValue, huffLog));
//
//   while (state.KeepRunning()) {
//     benchmark::ClobberMemory();
//     size_t const size = HUF_estimateCompressedSize(CTable, count, maxSymbolValue);
//     benchmark::DoNotOptimize(size);
//   }
// }
//
//
//
// BENCHMARK(BM_estimateBaseline);
// BENCHMARK(BM_estimate);
// BENCHMARK(BM_ValidateCTableBaseline);
// BENCHMARK(BM_ValidateCTable);
BENCHMARK(BM_CompressFlush)->RangeMultiplier(2)->Range(128, 2048);

BENCHMARK_MAIN();
