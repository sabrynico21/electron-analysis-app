# FastSpar bundled binaries

Place FastSpar executables in platform/architecture folders before packaging release builds.

FastSpar's permutation p-values need **three** executables per platform: `fastspar`
(correlation), `fastspar_bootstrap` (generates permuted count tables) and
`fastspar_pvalues` (computes p-values from the bootstrap correlations).

Required files:

- `linux/x64/fastspar`
- `linux/x64/fastspar_bootstrap`
- `linux/x64/fastspar_pvalues`
- `linux/arm64/fastspar`
- `linux/arm64/fastspar_bootstrap`
- `linux/arm64/fastspar_pvalues`
- `darwin/x64/fastspar`
- `darwin/x64/fastspar_bootstrap`
- `darwin/x64/fastspar_pvalues`
- `darwin/arm64/fastspar`
- `darwin/arm64/fastspar_bootstrap`
- `darwin/arm64/fastspar_pvalues`
- `win32/x64/fastspar.exe`
- `win32/x64/fastspar_bootstrap.exe`
- `win32/x64/fastspar_pvalues.exe`

Ensure Linux/macOS binaries are executable (`chmod +x`).

## Building

- Linux x64: fully static build in a Debian container (see `scripts/` and repo README).
- Windows x64: cross-compiled with MinGW-w64 — run `scripts/build-fastspar-win.sh`
  in a Debian container. It cross-compiles static gsl/openblas/arpack, fetches
  Armadillo 15.0.1 headers, and installs a POSIX `glob()` shim (MinGW lacks
  `glob.h`, required by `pvalue.cpp`).
