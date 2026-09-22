#!/usr/bin/env bash
# Cross-compile FastSpar v1.0.0 for Windows x86_64 (MinGW-w64, fully static)
# inside a Debian container. Produces in /out:
#   fastspar.exe  fastspar_bootstrap.exe  fastspar_pvalues.exe
#
# Usage (from repo root):
#   docker run --rm -d --name fastspar-win-build \
#     -v "$PWD/scripts:/scripts" \
#     -v "$PWD/resources/bin/win32/x64:/out" \
#     debian:bookworm bash /scripts/build-fastspar-win.sh
set -euo pipefail

TGT=x86_64-w64-mingw32
CC=${TGT}-gcc-posix
CXX=${TGT}-g++-posix
FC=${TGT}-gfortran-posix
PREFIX=/opt/mingw
BUILD=/build
JOBS=$(nproc)
[ "$JOBS" -gt 8 ] && JOBS=8

echo "==> Installing toolchain..."
export DEBIAN_FRONTEND=noninteractive
apt-get update -qq
apt-get install -y -qq \
  gcc g++ gfortran make cmake \
  autoconf automake libtool pkg-config \
  gcc-mingw-w64-x86-64-posix g++-mingw-w64-x86-64-posix gfortran-mingw-w64-x86-64-posix \
  wget xz-utils perl file git

mkdir -p "$BUILD" "$PREFIX/include" "$PREFIX/lib" /out
cd "$BUILD"

# ---------- GSL 2.8 ----------
echo "==> Building GSL 2.8"
[ -f gsl-2.8.tar.gz ] || wget -q https://ftp.gnu.org/gnu/gsl/gsl-2.8.tar.gz
tar xf gsl-2.8.tar.gz
cd gsl-2.8
./configure --host=$TGT --prefix=$PREFIX --disable-shared --enable-static CC=$CC
make -j"$JOBS"
make install
cd "$BUILD"

# ---------- OpenBLAS 0.3.28 ----------
echo "==> Building OpenBLAS 0.3.28"
[ -f OpenBLAS-0.3.28.tar.gz ] || \
  wget -q https://github.com/xianyi/OpenBLAS/archive/refs/tags/v0.3.28.tar.gz -O OpenBLAS-0.3.28.tar.gz
tar xf OpenBLAS-0.3.28.tar.gz
cd OpenBLAS-0.3.28
make -j"$JOBS" TARGET=HASWELL BINARY=64 USE_THREAD=0 NO_LAPACKE=1 NO_SHARED=1 \
     CC=$CC FC=$FC HOSTCC=gcc C_COMPILER=GCC F_COMPILER=GFORTRAN
make PREFIX=$PREFIX NO_SHARED=1 install
ln -sf libopenblas.a "$PREFIX/lib/libblas.a"
ln -sf libopenblas.a "$PREFIX/lib/liblapack.a"
cd "$BUILD"

# ---------- ARPACK-ng 3.9.1 ----------
echo "==> Building ARPACK-ng 3.9.1"
[ -f arpack-ng-3.9.1.tar.gz ] || \
  wget -q https://github.com/opencollab/arpack-ng/archive/refs/tags/3.9.1.tar.gz -O arpack-ng-3.9.1.tar.gz
tar xf arpack-ng-3.9.1.tar.gz
cat > toolchain-win64.cmake <<EOF
set(CMAKE_SYSTEM_NAME Windows)
set(CMAKE_SYSTEM_PROCESSOR x86_64)
set(CMAKE_C_COMPILER $CC)
set(CMAKE_CXX_COMPILER $CXX)
set(CMAKE_Fortran_COMPILER $FC)
set(CMAKE_FIND_ROOT_PATH $PREFIX /usr/$TGT)
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
EOF
cmake -S arpack-ng-3.9.1 -B arpack-ng-3.9.1/build \
  -DCMAKE_TOOLCHAIN_FILE=$BUILD/toolchain-win64.cmake \
  -DCMAKE_TRY_COMPILE_TARGET_TYPE=STATIC_LIBRARY \
  -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX=$PREFIX \
  -DBLAS_LIBRARIES=$PREFIX/lib/libopenblas.a \
  -DLAPACK_LIBRARIES=$PREFIX/lib/libopenblas.a
cmake --build arpack-ng-3.9.1/build -j"$JOBS"
cmake --install arpack-ng-3.9.1/build
cd "$BUILD"

# ---------- Armadillo 15.0.1 (header-only) ----------
echo "==> Fetching Armadillo 15.0.1 headers"
[ -f armadillo-15.0.1.tar.xz ] || \
  wget -q https://downloads.sourceforge.net/project/arma/armadillo-15.0.1.tar.xz -O armadillo-15.0.1.tar.xz
tar xf armadillo-15.0.1.tar.xz
cp -r armadillo-15.0.1/include/* "$PREFIX/include/"

# ---------- POSIX glob() shim (pvalue.cpp needs it; MinGW has no glob.h) ----------
echo "==> Building glob() shim for MinGW"
"$CC" -O2 -c /scripts/win-glob/glob.c -I/scripts/win-glob -o "$BUILD/glob.o"
cp /scripts/win-glob/glob.h "$PREFIX/include/glob.h"
ar rcs "$PREFIX/lib/libgsl.a" "$BUILD/glob.o"

# ---------- FastSpar v1.0.0 ----------
echo "==> Building FastSpar v1.0.0"
git clone --depth 1 --branch v1.0.0 https://github.com/scwatts/fastspar.git fastspar-src
cd fastspar-src
./autogen.sh
./configure --host=$TGT --enable-static-compile \
  CC=$CC CXX=$CXX F77=$FC \
  CPPFLAGS="-I$PREFIX/include" LDFLAGS="-L$PREFIX/lib" \
  CXXFLAGS="-std=c++14 -O2"
make -j"$JOBS"
cp -f src/fastspar.exe src/fastspar_bootstrap.exe src/fastspar_pvalues.exe /out/

echo "==> DONE. Output in /out:"
ls -la /out
file /out/*.exe
