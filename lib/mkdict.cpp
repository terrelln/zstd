#define ZDICT_STATIC_LINKING_ONLY
#include <algorithm>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <folly/Conv.h>
#include <folly/FileUtil.h>
#include <folly/String.h>
#include <gflags/gflags.h>
#include <glog/logging.h>
#include <random>
#include <string>
#include <vector>
#include <zdict.h>

DEFINE_string(dict_strategy, "", "");
DEFINE_int64(dict_size, 112640, "");
DEFINE_string(output_path, "", "Dict output path");
DEFINE_string(directory, "", "Directory that contains dictionary entires");

namespace fs = boost::filesystem;

static void init(int *argc, char ***argv) {
  gflags::ParseCommandLineFlags(argc, argv, true);
  FLAGS_logtostderr = true;
  auto programName = argc && argv && *argc > 0 ? (*argv)[0] : "unknown";
  google::InitGoogleLogging(programName);
}

int main(int argc, char **argv) {
  init(&argc, &argv);
  CHECK(!FLAGS_output_path.empty());
  LOG(INFO) << "Reading files from directory: " << FLAGS_directory;
  const auto dir = fs::system_complete(FLAGS_directory);
  CHECK(fs::exists(dir) && fs::is_directory(dir));

  std::vector<std::string> files;
  for (auto it = fs::recursive_directory_iterator{dir},
            end = fs::recursive_directory_iterator{};
       it != end; ++it) {
    if (!fs::is_regular_file(it->status())) {
      continue;
    }
    files.emplace_back();
    CHECK(folly::readFile(it->path().c_str(), files.back()))
        << it->path().native();
    if (files.back().size() > (size_t)FLAGS_dict_size) {
      files.back().resize(FLAGS_dict_size);
    }
  }
  LOG(INFO) << "Read " << files.size() << " files";

  //std::default_random_engine rng(1729);
  //std::shuffle(files.begin(), files.end(), rng);

  std::string data;
  std::vector<size_t> sizes;
  sizes.reserve(files.size());
  for (const auto& file : files) {
    data += file;
    sizes.push_back(file.size());
  }

  std::vector<folly::StringPiece> docs;
  size_t offset = 0;
  for (auto size : sizes) {
    docs.push_back(folly::StringPiece(data, offset, size));
    offset += size;
  }

  LOG(INFO) << "Building dictionary";
  std::string dict;

  std::vector<folly::StringPiece> strategyParts;
  folly::split(":", FLAGS_dict_strategy, strategyParts);
  if (strategyParts[0] == "COVER") {
    CHECK_EQ(strategyParts.size(), 6);
    COVER_params_t params;
    params.smoothing = folly::to<unsigned>(strategyParts[1]);
    params.kMin = folly::to<unsigned>(strategyParts[2]);
    params.kStep = folly::to<unsigned>(strategyParts[3]);
    params.kMax = folly::to<unsigned>(strategyParts[4]);
    params.d = folly::to<unsigned>(strategyParts[5]);
    params.notificationLevel = 5;
    params.dictID = 0;
    params.compressionLevel = 5;
    dict.resize(FLAGS_dict_size);
    const size_t dictSize = COVER_trainFromBuffer(
        &dict[0], dict.size(), data.data(), sizes.data(), sizes.size(), params);
    CHECK(!ZDICT_isError(dictSize)) << ZDICT_getErrorName(dictSize);
    CHECK_LE(dictSize, dict.size());
    dict.resize(dictSize);
  } else {
    LOG(FATAL) << "Invalid Strategy: " << FLAGS_dict_strategy;
  }

  CHECK(folly::writeFile(dict, FLAGS_output_path.c_str()));
}
