/*
 * Minimal POSIX glob() shim for MinGW-w64.
 * Only supports the features FastSpar's pvalue.cpp needs:
 *   glob(pattern, GLOB_TILDE, NULL, &glob_t)
 *   glob_t { gl_pathc, gl_pathv, ... }
 *   globfree(&glob_t)
 * Implemented with the Win32 FindFirstFileA/FindNextFileA API.
 */
#ifndef MINGW_GLOB_SHIM_H
#define MINGW_GLOB_SHIM_H

#include <stddef.h>

#define GLOB_TILDE    (1 << 0)
#define GLOB_NOMATCH  3

typedef struct {
    size_t gl_pathc;   /* Count of matched paths */
    char  **gl_pathv;  /* List of matched paths (NULL-terminated) */
    size_t gl_offs;    /* Reserved for GLOB_* flags */
    int    gl_flags;
    void  *gl_pglob;   /* Reserved */
} glob_t;

#ifdef __cplusplus
extern "C" {
#endif

int  glob(const char *pattern, int flags,
          int (*errfunc)(const char *, int), glob_t *pglob);
void globfree(glob_t *pglob);

#ifdef __cplusplus
}
#endif

#endif /* MINGW_GLOB_SHIM_H */
