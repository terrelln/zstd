/*
 * Copyright (c) Facebook, Inc.
 * All rights reserved.
 *
 * This source code is licensed under both the BSD-style license (found in the
 * LICENSE file in the root directory of this source tree) and the GPLv2 (found
 * in the COPYING file in the root directory of this source tree).
 * You may select, at your option, one of the above-listed licenses.
 */
#ifndef ZSTD_DECOMPRESS_BLOCK_MACROS_H
#define ZSTD_DECOMPRESS_BLOCK_MACROS_H

#define ZSTD_ENTROPY_LLTABLE_OFF 0
#define ZSTD_ENTROPY_MLTABLE_OFF 1
#define ZSTD_ENTROPY_OFTABLE_OFF 2
#define ZSTD_ENTROPY_REPCODE_OFF 3

#define ZSTD_ARG_ENTROPY_OFF 0
#define ZSTD_ARG_BITCACHE_OFF 8
#define ZSTD_ARG_BITPTR_OFF 16
#define ZSTD_ARG_BITLIMIT_OFF 24
#define ZSTD_ARG_BITSTART_OFF 32
#define ZSTD_ARG_LLSTATE_OFF 40
#define ZSTD_ARG_MLSTATE_OFF 48
#define ZSTD_ARG_OFSTATE_OFF 56
#define ZSTD_ARG_OP_OFF 64
#define ZSTD_ARG_OLIMIT_OFF 72
#define ZSTD_ARG_OEND_OFF 80
#define ZSTD_ARG_LITS_OFF 88
#define ZSTD_ARG_LITSLIMIT_OFF 96
#define ZSTD_ARG_LITSEND_OFF 104
#define ZSTD_ARG_PREFIXSTART_OFF 112
#define ZSTD_ARG_SAVEDOFFSET_OFF 116
#define ZSTD_ARG_SAVEDLITLEN_OFF 120
#define ZSTD_ARG_SAVEDMATCHLEN_OFF 124

#endif
