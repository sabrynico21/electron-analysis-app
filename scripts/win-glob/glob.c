/*
 * Minimal POSIX glob() shim for MinGW-w64 (see glob.h).
 * Uses Win32 FindFirstFileA/FindNextFileA to enumerate the directory
 * portion of the pattern and a simple '*'/'?' wildcard matcher for the
 * file-name portion.
 */
#include "glob.h"

#include <windows.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

static int match_pattern(const char *pat, const char *str)
{
    while (*pat) {
        if (*pat == '*') {
            while (*pat == '*') pat++;
            if (*pat == '\0') return 1;
            for (const char *s = str;; s++) {
                if (match_pattern(pat, s)) return 1;
                if (*s == '\0') break;
            }
            return 0;
        } else if (*pat == '?') {
            if (*str == '\0') return 0;
            pat++; str++;
        } else {
            if (*pat != *str) return 0;
            pat++; str++;
        }
    }
    return *str == '\0';
}

int glob(const char *pattern, int flags,
         int (*errfunc)(const char *, int), glob_t *pglob)
{
    (void)flags;
    (void)errfunc;
    if (!pglob) return GLOB_NOMATCH;
    pglob->gl_pathc = 0;
    pglob->gl_pathv = NULL;

    /* Split into directory part and file-name pattern at the last separator. */
    const char *sep = NULL;
    for (const char *p = pattern; *p; p++) {
        if (*p == '/' || *p == '\\') sep = p;
    }

    char dirbuf[4096];
    char filepat[4096];

    if (sep) {
        size_t n = (size_t)(sep - pattern);
        if (n >= sizeof(dirbuf)) n = sizeof(dirbuf) - 1;
        memcpy(dirbuf, pattern, n);
        dirbuf[n] = '\0';
        if (n == 0) { dirbuf[0] = '.'; dirbuf[1] = '\0'; }
        _snprintf(filepat, sizeof(filepat), "%s", sep + 1);
    } else {
        dirbuf[0] = '.'; dirbuf[1] = '\0';
        _snprintf(filepat, sizeof(filepat), "%s", pattern);
    }

    char search[4096];
    _snprintf(search, sizeof(search), "%s\\*", dirbuf);

    WIN32_FIND_DATAA fd;
    HANDLE h = FindFirstFileA(search, &fd);
    if (h == INVALID_HANDLE_VALUE) {
        return GLOB_NOMATCH;
    }

    char **paths = NULL;
    size_t count = 0;
    do {
        if (match_pattern(filepat, fd.cFileName)) {
            char full[4096];
            _snprintf(full, sizeof(full), "%s\\%s", dirbuf, fd.cFileName);
            char **np = (char **)realloc(paths, (count + 1) * sizeof(char *));
            if (!np) break;
            paths = np;
            paths[count] = _strdup(full);
            count++;
        }
    } while (FindNextFileA(h, &fd));
    FindClose(h);

    pglob->gl_pathc = count;
    pglob->gl_pathv = paths;
    return count ? 0 : GLOB_NOMATCH;
}

void globfree(glob_t *pglob)
{
    if (!pglob) return;
    if (pglob->gl_pathv) {
        size_t i;
        for (i = 0; i < pglob->gl_pathc; i++) {
            free(pglob->gl_pathv[i]);
        }
        free(pglob->gl_pathv);
        pglob->gl_pathv = NULL;
    }
    pglob->gl_pathc = 0;
}
