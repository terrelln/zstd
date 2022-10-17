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
#define llState     r13
#define mlState     r14
#define ofState     r15

#define llState_d   r14d
#define llState_w   r14w
#define llState_b   r14b

#define mlState_d   r15d
#define mlState_w   r15w
#define mlState_b   r15b

#define ofState_d   r15d
#define ofState_w   r15w
#define ofState_b   r15b

// rax
// rbx
// rcx
// rdx
// rdi

// rbp?

#define ZSTD_STACK_BITLIMIT_OFF 0
#define ZSTD_STACK_OLIMIT_OFF 8
#define ZSTD_STACK_LITSLIMIT_OFF 16
#define ZSTD_STACK_PREFIXSTART_OFF 24
#define ZSTD_STACK_NBSEQS_OFF 32
#define ZSTD_STACK_ARGS_OFF 40

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
    movq ZSTD_ARG_LLSTATE_OFF(%rdi), %llState
    movq ZSTD_ARG_MLSTATE_OFF(%rdi), %mlState
    movq ZSTD_ARG_OFSTATE_OFF(%rdi), %ofState

    movq ZSTD_ARG_BITLIMIT_OFF(%rdi), %rax
    movq ZSTD_ARG_OLIMIT_OFF(%rdi), %rbx
    movq ZSTD_ARG_LITSLIMIT_OFF(%rdi), %rcx
    movq ZSTD_ARG_PREFIXSTART_OFF(%rdi), %rsi
    movq $0, %rdx

    push %rdi
    push %rdx
    push %rsi
    push %rcx
    push %rbx
    push %rax

/* Add bits to reg. Negates the nbBits_b register. */
#define ZSTD_ADD_BITS_TO(reg, nbBits, nbBits_b) \
    shrxq %nbBits, %bitCache, %rdx;             \
    neg %nbBits_b;                              \
    shlxq %nbBits, %bitCache, %bitCache;        \
    addq %reg, %rdx

.L_loop_body:

    /* ofDInfo */
    movq ZSTD_ENTROPY_OFTABLE_OFF(%entropy, %ofState, 8), %ofState
    /* llDInfo */
    movq ZSTD_ENTROPY_LLTABLE_OFF(%entropy, %llState, 8), %llState
    /* mlDInfo */
    movq ZSTD_ENTROPY_MLTABLE_OFF(%entropy, %mlState, 8), %mlState

    /* %rax = ofDInfo */
    movq %ofState, %rax

    /* %rbx = llDInfo */
    movq %llState, %rbx 

    /* %rcx = mlDInfo */
    movq %mlState, %rcx

    movq $32, %rdx
    /* %rax = ofDInfo.baseValue */
    shrxq %rdx, %rax, %rax
    /* %rbx = llDInfo.baseValue */
    shrxq %rdx, %rbx, %rbx
    /* %rcx = mlDInfo.baseValue */
    shrxq %rdx, %rcx, %rcx

    test %ofState_b, %ofState_b
    jz .L_rep0

    ZSTD_ADD_BITS_TO(rax, ofState, ofState_b)

    movd ZSTD_ENTROPY_REPCODE2_OFF(%entropy), ZSTD_ENTROPY_REPCODE3_OFF(%entropy)
    movd ZSTD_ENTROPY_REPCODE1_OFF(%entropy), ZSTD_ENTROPY_REPCODE2_OFF(%entropy)
    movd %ofState_d, ZSTD_ENTROPY_REPCODE1_OFF(%entropy)

.L_offset_decoded:

    test %cl, %cl
    jnz .L_additional_ml_bits

.L_ml_decoded:

    test %dl, %dl
    jnz .L_additional_ll_bits

.L_ll_decoded:

    movq $16, %rsi
    movq $8, %rdi
    shrxq %rdi, %llState, %rdx
    shrq %rsi, %llState, %llState
    movzx %llState_w, %llState_d

    ZSTD_ADD_BITS_TO(%llState, rdx, dl)

    shrxq %rdi, %mlState, %rdx
    shrq %rsi, %mlState, %mlState
    movzx %mlState_w, %mlState_d

    ZSTD_ADD_BITS_TO(%mlState, rdx, dl)

    shrxq %rdi, %ofState, %rdx
    shrq %rsi, %ofState, %ofState
    movzx %ofState_w, %ofState_d

    ZSTD_ADD_BITS_TO(%ofState, %rdx, dl)

    /* Reload bitstream */
    bsfq %bitCache, %bitCache
    movq %bitCache, %rdx
    shrq $3, %bitCache
    andq $7, %rdx
    subq %bitCache, %bitPtr
    movq (%bitPtr), %bitCache
    orq $1, %bitCache
    shlxq %rdx, %bitCache, %bitCache

    /* Validate sequence */
    movq %op, %rdx
    /* rdx = match */
    subq %rax, %rdx
    /* rsi = oSeqEnd */
    leaq (%op, %rbx, 1), %rsi
    leaq (%rdx, %rcx, 1), %si
    /* rdi = iLitEnd */
    leaq (%lits, %rbx, 1), %rdi

    /* sil = oSeqEnd > olimit */
    cmpq ZSTD_STACK_OLIMIT_OFF(%rsp), %rsi
    seta %sil
    /* dil = iLitEnd > litsLimit */
    cmpq ZSTD_STACK_LITSLIMIT_OFF(%rsp), %rdi
    seta %dil
    /* sil = sil | dil */
    or %dil, %sil
    /* match < prefixStart */
    cmpq ZSTD_STACK_PREFIXSTART_OFF(%rsp), %rdx
    setb %dil
    /* sil = sil | dil */
    or %dil, %sil
    jne .L_save_sequence

    /* ++nbSeq */
    addq $1, ZSTD_STACK_NBSEQS_OFF(%rsp)
    
    /* Copy literals */
    movups (%lits), %xmm0
    movups %xmm0, (%op)
    cmpq $16, %rbx
    ja .L_long_literals
    addq %rbx, %op

.L_literals_copied:
    addq %rbx, %lits 
    
    /* Copy match (handling short offset) */
    movups (%rdx), %xmm0
    movups 16(%rdx), %xmm1

    cmpq $16, %rax
    jb .L_short_offset_match:

    movups %xmm0, (%op)
    movups %xmm1 16(%op)

    cmpq $32, %rcx
    ja .L_long_match
    addq %rcx, %op

.L_match_copied:

    /* Check bitPtr >= bitLimit */
    cmp ZSTD_STACK_BITLIMIT_OFF(%rsp), %bitPtr
    jae .L_loop_body

    jmp .L_exit

.L_rep0:
    test %bl, %bl
    jz .L_rep0_ll0

    movd ZSTD_ENTROPY_REPCODE1_OFF(%entropy), %ofState_d
    jmp .L_offset_decoded

.L_rep0_ll0:
    movd ZSTD_ENTROPY_REPCODE2_OFF(%entropy), %ofState_d
    movd ZSTD_ENTROPY_REPCODE1_OFF(%entropy), ZSTD_ENTROPY_REPCODE2_OFF(%entropy)
    movd %ofState_d, ZSTD_ENTROPY_REPCODE1_OFF(%entropy)
    jmp .L_offset_decoded

.L_additional_ml_bits:
    ZSTD_ADD_BITS_TO(rcx, mlState, mlState_b)
    jmp .L_ml_decoded

.L_additional_ll_bits:
    ZSTD_ADD_BITS_TO(rbx, llState, llState_b)
    jmp .L_ll_decoded

.L_long_literals:
    /* scratch: rsi, rdi */
    leaq 16(%op), %rdi
    leaq 16(%lits), %rsi
    addq %rbx, %op

.L_long_literals_loop:
    movups (%rsi), %xmm0
    movups 16(%rsi), %xmm1
    movups %xmm0, (%rdi)
    movups %xmm1, 16(%rdi)
    addq $32, %rsi
    addq $32, %rdi
    cmp %op, %rdi
    jb .L_long_literals_loop

    jmp .L_literals_copied

.L_long_match:
    /* scratch: rbx, rsi, rdi */
    leaq 16(%op), %rdi
    leqa 16(%rdx), %rsi
    addq %rcx, %op

.L_long_match_loop:
    movups (%rsi), %xmm0
    movups 16(%rsi), %xmm1
    movups %xmm0, (%rdi)
    movups %xmm1, 16(%rdi)
    addq $32, %rsi
    addq $32, %rdi
    cmp %op, %rdi
    jb .L_long_match_loop

    jmp .L_match_copied

.L_short_offset_match:
    /* offset = %rax
     * match = %rdx
     * matchLen = %rcx
     * scratch: rbx, rsi, rdi
     */
    /* TODO: Read only from reg, not match,
     * to avoid store forwarding issues.
     */
    leaq .L_pattern(%rip), %rbx
    shlq $4, %rax
    vpshufb (%rbx, %rax), %xmm0, %xmm0
    movups %xmm0, %op

    cmpq $16, %rcx
    ja .L_short_offset_match_long
    
    addq %rcx, %op
    jmp .L_match_copied

.L_short_offset_match_long:
    leaq 16(%op), %rdi
    leaq .L_pattern_forward(%rip), %rbx
    movq %rdx, %rsi
    addq (%rbx, %rax), %rsi

    addq %rcx, %op

    jmp .L_long_match_loop

.L_save_sequence:
    /* rax = offset
     * rbx = litLen
     * rcx = matchLen
     */
    /* rdi = args */
    movq ZSTD_STACK_ARGS_OFF(%rsp), %rdi
    movd %rax, ZSTD_ARGS_SAVEDOFFSET_OFF(%rdi)
    movd %rbx, ZSTD_ARGS_SAVEDLITLEN_OFF(%rdi)
    movd %rcx, ZSTD_ARGS_SAVEDMATCHLEN_OFF(%rdi)

.L_exit:

    pop %rdi /* args */
    pop %rax /* litslimit */
    pop %rax /* olimit */
    pop %rax /* bitlimit */
    pop %rax /* nbSeqs */

    movq %bitCache, ZSTD_ARG_BITCACHE_OFF(%rdi)
    movq %bitPtr, ZSTD_ARG_BITPTR_OFF(%rdi)
    movq %llState, ZSTD_ARG_LLSTATE_OFF(%rdi)
    movq %mlState, ZSTD_ARG_MLSTATE_OFF(%rdi)
    movq %ofState, ZSTD_ARG_OFSTATE_OFF(%rdi)
    movq %op, ZSTD_ARG_OP_OFF(%rdi)
    movq %lits, ZSTD_ARG_LITS_OFF(%rdi)

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

        .section        .rodata
	.align 8
.L_pattern_forward:
	.quad 0 /* 0 */
	.quad 0 /* padding */
	.quad 0 /* 1 */
	.quad 0 /* padding */
	.quad 0 /* 2 */
	.quad 0 /* padding */
	.quad 1 /* 3 */
	.quad 0 /* padding */
	.quad 0 /* 4 */
	.quad 0 /* padding */
	.quad 1 /* 5 */
	.quad 0 /* padding */
	.quad 4 /* 6 */
	.quad 0 /* padding */
	.quad 2 /* 7 */
	.quad 0 /* padding */
	.quad 0 /* 8 */
	.quad 0 /* padding */
	.quad 7 /* 9 */
	.quad 0 /* padding */
	.quad 6 /* 10 */
	.quad 0 /* padding */
	.quad 5 /* 11 */
	.quad 0 /* padding */
	.quad 4 /* 12 */
	.quad 0 /* padding */
	.quad 3 /* 13 */
	.quad 0 /* padding */
	.quad 2 /* 14 */
	.quad 0 /* padding */
	.quad 1 /* 15 */
	.align 16
.L_pattern:
	.quad   0 /* padding */
	.quad   0 /* padding */
        .byte   0
        .byte   0
        .byte   0
        .byte   0
        .byte   0
        .byte   0
        .byte   0
        .byte   0
        .byte   0
        .byte   0
        .byte   0
        .byte   0
        .byte   0
        .byte   0
        .byte   0
        .byte   0
        .byte   0
        .byte   1
        .byte   0
        .byte   1
        .byte   0
        .byte   1
        .byte   0
        .byte   1
        .byte   0
        .byte   1
        .byte   0
        .byte   1
        .byte   0
        .byte   1
        .byte   0
        .byte   1
        .byte   0
        .byte   1
        .byte   2
        .byte   0
        .byte   1
        .byte   2
        .byte   0
        .byte   1
        .byte   2
        .byte   0
        .byte   1
        .byte   2
        .byte   0
        .byte   1
        .byte   2
        .byte   0
        .byte   0
        .byte   1
        .byte   2
        .byte   3
        .byte   0
        .byte   1
        .byte   2
        .byte   3
        .byte   0
        .byte   1
        .byte   2
        .byte   3
        .byte   0
        .byte   1
        .byte   2
        .byte   3
        .byte   0
        .byte   1
        .byte   2
        .byte   3
        .byte   4
        .byte   0
        .byte   1
        .byte   2
        .byte   3
        .byte   4
        .byte   0
        .byte   1
        .byte   2
        .byte   3
        .byte   4
        .byte   0
        .byte   0
        .byte   1
        .byte   2
        .byte   3
        .byte   4
        .byte   5
        .byte   0
        .byte   1
        .byte   2
        .byte   3
        .byte   4
        .byte   5
        .byte   0
        .byte   1
        .byte   2
        .byte   3
        .byte   0
        .byte   1
        .byte   2
        .byte   3
        .byte   4
        .byte   5
        .byte   6
        .byte   0
        .byte   1
        .byte   2
        .byte   3
        .byte   4
        .byte   5
        .byte   6
        .byte   0
        .byte   1
        .byte   0
        .byte   1
        .byte   2
        .byte   3
        .byte   4
        .byte   5
        .byte   6
        .byte   7
        .byte   0
        .byte   1
        .byte   2
        .byte   3
        .byte   4
        .byte   5
        .byte   6
        .byte   7
        .byte   0
        .byte   1
        .byte   2
        .byte   3
        .byte   4
        .byte   5
        .byte   6
        .byte   7
        .byte   8
        .byte   0
        .byte   1
        .byte   2
        .byte   3
        .byte   4
        .byte   5
        .byte   6
        .byte   0
        .byte   1
        .byte   2
        .byte   3
        .byte   4
        .byte   5
        .byte   6
        .byte   7
        .byte   8
        .byte   9
        .byte   0
        .byte   1
        .byte   2
        .byte   3
        .byte   4
        .byte   5
        .byte   0
        .byte   1
        .byte   2
        .byte   3
        .byte   4
        .byte   5
        .byte   6
        .byte   7
        .byte   8
        .byte   9
        .byte   10
        .byte   0
        .byte   1
        .byte   2
        .byte   3
        .byte   4
        .byte   0
        .byte   1
        .byte   2
        .byte   3
        .byte   4
        .byte   5
        .byte   6
        .byte   7
        .byte   8
        .byte   9
        .byte   10
        .byte   11
        .byte   0
        .byte   1
        .byte   2
        .byte   3
        .byte   0
        .byte   1
        .byte   2
        .byte   3
        .byte   4
        .byte   5
        .byte   6
        .byte   7
        .byte   8
        .byte   9
        .byte   10
        .byte   11
        .byte   12
        .byte   0
        .byte   1
        .byte   2
        .byte   0
        .byte   1
        .byte   2
        .byte   3
        .byte   4
        .byte   5
        .byte   6
        .byte   7
        .byte   8
        .byte   9
        .byte   10
        .byte   11
        .byte   12
        .byte   13
        .byte   0
        .byte   1
        .byte   0
        .byte   1
        .byte   2
        .byte   3
        .byte   4
        .byte   5
        .byte   6
        .byte   7
        .byte   8
        .byte   9
        .byte   10
        .byte   11
        .byte   12
        .byte   13
        .byte   14
        .byte   0
        .byte   0
        .byte   1
        .byte   2
        .byte   3
        .byte   4
        .byte   5
        .byte   6
        .byte   7
        .byte   8
        .byte   9
        .byte   10
        .byte   11
        .byte   12
        .byte   13
        .byte   14
        .byte   15
	.text

#endif

/**
 * IDEAS:
 * - Reorder ZSTD_seqSymbol to allow easier register access
 *   Different tables can be in different orders depending on their need.
 *   E.g. literals have baseValue in low. off/match have nbAddBits in low.
 *
 * - If we need more registers prefixStart can live on stack
 */