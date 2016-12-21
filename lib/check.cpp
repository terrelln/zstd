#define ZDICT_STATIC_LINKING_ONLY
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <folly/FileUtil.h>
#include <gflags/gflags.h>
#include <glog/logging.h>
#include <string>
#include <vector>
#include <zstd.h>

DEFINE_string(directory, "", "Directory that contains dictionary entires");
DEFINE_string(dictionary, "", "The dictionary");

namespace fs = boost::filesystem;

static void init(int *argc, char ***argv) {
  gflags::ParseCommandLineFlags(argc, argv, true);
  FLAGS_logtostderr = true;
  auto programName = argc && argv && *argc > 0 ? (*argv)[0] : "unknown";
  google::InitGoogleLogging(programName);
}

int main(int argc, char **argv) {
  init(&argc, &argv);
  LOG(INFO) << "Reading files from directory: " << FLAGS_directory;
  const auto dir = fs::system_complete(FLAGS_directory);
  CHECK(fs::exists(dir) && fs::is_directory(dir));

  ZSTD_CDict *cdict = nullptr;
  {
    std::string data;
    CHECK(folly::readFile(FLAGS_dictionary.c_str(), data)) << FLAGS_dictionary;
    cdict = ZSTD_createCDict(data.data(), data.size(), 5);
    CHECK_NOTNULL(cdict);
  }
  auto ctx = ZSTD_createCCtx();
  CHECK_NOTNULL(ctx);

  size_t rawSize = 0;
  size_t compressedSize = 0;
  std::string buf;
  std::string out;
  for (auto it = fs::recursive_directory_iterator{dir},
            end = fs::recursive_directory_iterator{};
       it != end; ++it) {
    if (!fs::is_regular_file(it->status())) {
      continue;
    }
    buf.clear();
    CHECK(folly::readFile(it->path().c_str(), buf)) << it->path().native();

    const size_t bound = ZSTD_compressBound(buf.size());
    if (out.size() < bound) {
      out.resize(bound);
    }
    const size_t rc = ZSTD_compress_usingCDict(ctx, &out[0], out.size(),
                                               buf.data(), buf.size(), cdict);
    CHECK(!ZSTD_isError(rc)) << ZSTD_getErrorName(rc);

    rawSize += buf.size();
    compressedSize += rc;
  }
  LOG(INFO) << "Raw size: " << rawSize;
  LOG(INFO) << "Compressed size: " << compressedSize;
  LOG(INFO) << "Compression ratio: "
            << (double)rawSize / (double)compressedSize;
}
