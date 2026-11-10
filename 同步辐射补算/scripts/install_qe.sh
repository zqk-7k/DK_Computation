#!/usr/bin/env bash
# Project-local MPI Quantum ESPRESSO. Existing site installation is preferred.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
mkdir -p "$ROOT/software" "$ROOT/logs"
if command -v pw.x >/dev/null 2>&1; then
  command -v pw.x > "$ROOT/software/qe_executable.txt"
  echo 'Using site Quantum ESPRESSO; record its banner from the first actual run.'
  exit 0
fi
QE_COMMIT=500de340b820e1cb8c05f2d8bb8fced102f377c1
SRC="$ROOT/software/qe-source"
BUILD="$ROOT/software/qe-build"
[[ "$(uname -s)" == Linux ]] || { echo 'Build on Linux HPC.' >&2; exit 2; }
for tool in git cmake make; do command -v "$tool" >/dev/null || { echo "Load/install $tool first." >&2; exit 2; }; done
if ! command -v mpifort >/dev/null && ! command -v mpif90 >/dev/null && ! command -v mpiifort >/dev/null; then
  echo 'Load matching Fortran and MPI modules before building. Serial fallback is not implicit.' >&2; exit 2
fi
if [[ ! -d "$SRC/.git" ]]; then
  git clone --depth 1 --branch qe-7.4.1 https://github.com/QEF/q-e.git "$SRC"
fi
[[ "$(git -C "$SRC" rev-parse HEAD)" == "$QE_COMMIT" ]] || { echo 'QE source commit differs; preserve existing tree and inspect.' >&2; exit 2; }
git -C "$SRC" submodule update --init --recursive --depth 1
cmake -S "$SRC" -B "$BUILD" -DCMAKE_BUILD_TYPE=Release \
  -DQE_ENABLE_MPI=ON -DQE_ENABLE_OPENMP=OFF -DQE_ENABLE_TEST=OFF \
  -DQE_ENABLE_SCALAPACK=OFF -DQE_FFTW_VENDOR=Internal \
  -DQE_LAPACK_INTERNAL=ON 2>&1 | tee "$ROOT/logs/qe_cmake.log"
cmake --build "$BUILD" --target pw --parallel "${XAS_BUILD_JOBS:-2}" 2>&1 | tee "$ROOT/logs/qe_build.log"
[[ -x "$BUILD/bin/pw.x" ]] || { echo 'Build finished without pw.x.' >&2; exit 2; }
printf '%s\n' "$BUILD/bin/pw.x" > "$ROOT/software/qe_executable.txt"
printf '%s\n' "$QE_COMMIT" > "$ROOT/software/qe_source_commit.txt"
echo "Quantum ESPRESSO ready: $BUILD/bin/pw.x"
