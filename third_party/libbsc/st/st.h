/*
 * Sort Transform stub — LIBBSC_SORT_TRANSFORM_SUPPORT not enabled.
 * Sort Transform (US patent 6,199,064) is disabled by default in libbsc.
 */

#ifndef LIBBSC_ST_H
#define LIBBSC_ST_H

#include "../platform/platform.h"

#ifdef __cplusplus
extern "C" {
#endif

static inline int bsc_st_init(int features) { (void)features; return 0; }

static inline int bsc_st3_encode(unsigned char *T, int n, int features) { (void)T; (void)n; (void)features; return -1; }
static inline int bsc_st4_encode(unsigned char *T, int n, int features) { (void)T; (void)n; (void)features; return -1; }
static inline int bsc_st5_encode(unsigned char *T, int n, int features) { (void)T; (void)n; (void)features; return -1; }
static inline int bsc_st6_encode(unsigned char *T, int n, int features) { (void)T; (void)n; (void)features; return -1; }
static inline int bsc_st7_encode(unsigned char *T, int n, int features) { (void)T; (void)n; (void)features; return -1; }
static inline int bsc_st8_encode(unsigned char *T, int n, int features) { (void)T; (void)n; (void)features; return -1; }

static inline int bsc_st_decode(unsigned char *T, int n, int k, int features) { (void)T; (void)n; (void)k; (void)features; return -1; }

#ifdef __cplusplus
}
#endif

#endif /* LIBBSC_ST_H */
