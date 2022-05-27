/*
 * Copyright (c) Facebook, Inc.
 * All rights reserved.
 *
 * This source code is licensed under both the BSD-style license (found in the
 * LICENSE file in the root directory of this source tree) and the GPLv2 (found
 * in the COPYING file in the root directory of this source tree).
 * You may select, at your option, one of the above-listed licenses.
 */

#include "../common/portability_macros.h"
#include "zstd_decompress_block_macros.h"

/* Stack marking
 * ref: https://wiki.gentoo.org/wiki/Hardened/GNU_stack_quickstart
 */
#if defined(__ELF__) && defined(__GNUC__)
.section .note.GNU-stack,"",%progbits
#endif

#if ZSTD_ENABLE_ASM_X86_64_BMI2

/* Calling convention:
 *
 * %rdi contains the first argument: HUF_DecompressAsmArgs*.
 * %rbp isn't maintained (no frame pointer).
 * %rsp contains the stack pointer that grows down.
 *      No red-zone is assumed, only addresses >= %rsp are used.
 * All register contents are preserved.
 *
 * TODO: Support Windows calling convention.
 */

ZSTD_HIDE_ASM_FUNCTION(ZSTD_decompressSequences2_asm)
ZSTD_HIDE_ASM_FUNCTION(_ZSTD_decompressSequences2_asm)
.global ZSTD_decompressSequences2_asm
.global _ZSTD_decompressSequences2_asm
.text

/* Sets up register mappings for clarity.
 * op[], bits[], dtable & ip[0] each get their own register.
 * ip[1,2,3] & olimit alias var[].
 * %rax is a scratch register.
 */

#define entropy     r8
#define bitCache    r9
#define bitPtr      r10
#define op          r11
#define lits        r12
#define prefixStart r13
#define llState     r14
#define mlState     r15
#define ofState     rsi

// rax
// rbx
// rcx
// rdx
// rdi

// rbp?

#define ZSTD_STACK_BITLIMIT_OFF 0
#define ZSTD_STACK_OLIMIT_OFF 8
#define ZSTD_STACK_LITSLIMIT_OFF 16
#define ZSTD_STACK_NBSEQS_OFF 24

_ZSTD_decompressSequences2_asm:
ZSTD_decompressSequences2_asm:
    /* Save all registers except rax - even if they are callee saved for simplicity. */
    /* push %rax */
    push %rbx
    push %rcx
    push %rdx
    push %rbp
    push %rsi
    push %rdi
    push %r8
    push %r9
    push %r10
    push %r11
    push %r12
    push %r13
    push %r14
    push %r15

    movq ZSTD_ARG_ENTROPY_OFF(%rdi), %entropy
    movq ZSTD_ARG_BITCACHE_OFF(%rdi), %bitCache
    movq ZSTD_ARG_BITPTR_OFF(%rdi), %bitPtr
    movq ZSTD_ARG_OP_OFF(%rdi), %op
    movq ZSTD_ARG_LITS_OFF(%rdi), %lits
    movq ZSTD_ARG_PREFIXSTART_OFF(%rdi), %prefixStart
    movq ZSTD_ARG_LLSTATE_OFF(%rdi), %llState
    movq ZSTD_ARG_MLSTATE_OFF(%rdi), %mlState
    movq ZSTD_ARG_OFSTATE_OFF(%rdi), %ofState

    movq ZSTD_ARG_BITLIMIT_OFF(%rdi), %rax
    movq ZSTD_ARG_OLIMIT_OFF(%rdi), %rbx
    movq ZSTD_ARG_LITSLIMIT_OFF(%rdi), %rcx
    movq $0, %rdx

    push %rdx
    push %rcx
    push %rbx
    push %rax

    pop %rax /* litslimit */
    pop %rax /* olimit */
    pop %rax /* bitlimit */
    pop %rax /* nbSeqs */

    /* Restore registers */
    pop %r15
    pop %r14
    pop %r13
    pop %r12
    pop %r11
    pop %r10
    pop %r9
    pop %r8
    pop %rdi
    pop %rsi
    pop %rbp
    pop %rdx
    pop %rcx
    pop %rbx
    /* pop %rax */
    ret


#endif
