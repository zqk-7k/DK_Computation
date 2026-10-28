#!/usr/bin/env bash
# Build the pinned, complete FEFF10 serial engine without sudo or system changes.
set -Eeuo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
PREFIX="$ROOT/software/feff10"
COMMIT="d861cf8b64342e2d63a0aeb62c8b78faca7942c3"
ARCHIVE_SHA256="f3397e872ed376e37c0404c53f9c8996e0247044571b9a344a6f32b2932bda66"
URL="https://codeload.github.com/times-software/feff10/tar.gz/$COMMIT"

if [[ "${1:-}" == "--help" ]]; then
  cat <<'HELP'
Usage: bash scripts/install_feff10.sh [--check]
Installs the pinned official FEFF10 source under this package's software/feff10.
Requires Linux, Bash, make, tar, sha256sum, Python 3 and a Fortran compiler.
Load your cluster compiler module first; no sudo is used.
Optional environment: FEFF_FC=gfortran (or ifort/ifx), FEFF_FFLAGS='...'.
For an offline node: FEFF_SOURCE_ARCHIVE=/absolute/path/to/pinned-archive.tar.gz.
Compilation is serial (-j1); independent spectroscopy jobs run in parallel later.
HELP
  exit 0
fi
if [[ "${1:-}" == "--check" ]]; then
  [[ -x "$PREFIX/bin/run_feff10.sh" ]] || { echo "FEFF10 wrapper is not installed: $PREFIX" >&2; exit 1; }
  exec bash "$PREFIX/bin/run_feff10.sh" --check
fi
[[ $# -eq 0 ]] || { echo "Unknown argument; use --help." >&2; exit 2; }
[[ "$(uname -s)" == Linux ]] || { echo "This installer targets Linux HPC nodes." >&2; exit 2; }
for cmd in make tar sha256sum python3; do
  command -v "$cmd" >/dev/null || { echo "Missing $cmd. Load a cluster module; no sudo is attempted." >&2; exit 2; }
done
FC_BIN="${FEFF_FC:-gfortran}"
command -v "$FC_BIN" >/dev/null || { echo "Missing Fortran compiler: $FC_BIN. Load a compiler module or set FEFF_FC." >&2; exit 2; }
mkdir -p "$PREFIX"
if ! mkdir "$PREFIX/.install.lock" 2>/dev/null; then
  echo "Another installation may be active. Inspect $PREFIX/.install.lock before retrying." >&2
  exit 2
fi
trap 'rmdir "$PREFIX/.install.lock" 2>/dev/null || true' EXIT
exec > >(tee -a "$PREFIX/install.log") 2>&1
printf '\nFEFF10 installation: %s UTC\nPinned commit: %s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" "$COMMIT"
if [[ -f "$PREFIX/install_complete.json" && -x "$PREFIX/bin/run_feff10.sh" ]]; then
  if python3 - "$PREFIX/install_complete.json" "$COMMIT" <<'PY'
import json,sys
sys.exit(0 if json.load(open(sys.argv[1]))['commit'] == sys.argv[2] else 1)
PY
  then
    bash "$PREFIX/bin/run_feff10.sh" --check
    echo "Pinned installation already exists; no rebuild performed."
    exit 0
  fi
  echo "Existing completed installation has a different commit; refusing to overwrite." >&2
  exit 2
fi
mkdir -p "$PREFIX/downloads"
ARCHIVE="${FEFF_SOURCE_ARCHIVE:-$PREFIX/downloads/feff10-$COMMIT.tar.gz}"
if [[ ! -f "$ARCHIVE" ]]; then
  [[ -z "${FEFF_SOURCE_ARCHIVE:-}" ]] || { echo "Offline archive does not exist: $ARCHIVE" >&2; exit 2; }
  if command -v curl >/dev/null; then
    curl --fail --location --retry 3 --connect-timeout 30 "$URL" -o "$ARCHIVE.part"
  elif command -v wget >/dev/null; then
    wget --tries=3 --timeout=30 "$URL" -O "$ARCHIVE.part"
  else
    echo "Install requires curl/wget, or FEFF_SOURCE_ARCHIVE pointing to the pinned archive." >&2
    exit 2
  fi
  mv -- "$ARCHIVE.part" "$ARCHIVE"
fi
ACTUAL_SHA256="$(sha256sum "$ARCHIVE" | awk '{print $1}')"
[[ "$ACTUAL_SHA256" == "$ARCHIVE_SHA256" ]] || {
  echo "Archive SHA256 mismatch; refusing to build. Expected $ARCHIVE_SHA256, got $ACTUAL_SHA256" >&2
  exit 2
}
SOURCE="$PREFIX/source"
if [[ -d "$SOURCE" ]]; then
  [[ -f "$SOURCE/.workflow_commit" && "$(cat "$SOURCE/.workflow_commit")" == "$COMMIT" ]] || {
    echo "Unidentified source directory exists; preserve or rename it before retrying: $SOURCE" >&2; exit 2;
  }
else
  STAGE="$(mktemp -d "$PREFIX/source-stage.XXXXXXXX")"
  # The cryptographic digest above fixes the archive. Avoid extracting outside STAGE.
  python3 - "$ARCHIVE" <<'PY'
import sys,tarfile,pathlib
with tarfile.open(sys.argv[1], 'r:gz') as archive:
    for entry in archive:
        p=pathlib.PurePosixPath(entry.name)
        if p.is_absolute() or '..' in p.parts:
            raise SystemExit('Unsafe archive member: '+entry.name)
        if entry.issym() or entry.islnk():
            q=pathlib.PurePosixPath(entry.linkname)
            if q.is_absolute() or '..' in q.parts:
                raise SystemExit('Unsafe archive link: '+entry.name)
PY
  tar -xzf "$ARCHIVE" -C "$STAGE" --strip-components=1
  printf '%s\n' "$COMMIT" > "$STAGE/.workflow_commit"
  mv -- "$STAGE" "$SOURCE"
fi

FC_NAME="$(basename -- "$FC_BIN")"
case "$FC_NAME" in
  *gfortran*)
    FLAGS='-O2 -cpp -ffree-line-length-none'
    MAJOR="$("$FC_BIN" -dumpversion | cut -d. -f1)"
    if [[ "$MAJOR" =~ ^[0-9]+$ ]] && (( MAJOR >= 10 )); then
      FLAGS="$FLAGS -fallow-argument-mismatch"
    fi
    ;;
  ifort|ifx) FLAGS='-O2 -fpp -heap-arrays' ;;
  *) echo "Unknown compiler family '$FC_NAME'; use gfortran, ifort, or ifx." >&2; exit 2 ;;
esac
FLAGS="${FEFF_FFLAGS:-$FLAGS}"
[[ "$FC_BIN$FLAGS" != *$'\n'* ]] || { echo "Compiler settings cannot contain newlines." >&2; exit 2; }
cat > "$SOURCE/src/Compiler.mk" <<EOF
# Generated by this workflow. Built-in FEFF linear algebra, no MKL or MPI.
F90 = $FC_BIN
FLAGS = $FLAGS
# Required even for literal gfortran: DMDW must run as a FEFF module.
FPPTASK = -DFEFF
DEPTYPE =
LDFLAGS =
FCINCLUDE =
EOF
printf 'Compiler: %s\nFlags: %s\n' "$FC_BIN" "$FLAGS"
"$FC_BIN" --version | head -n 1
# Explicit all is required: the upstream first target can depend on included files.
# Build serially because Fortran module dependencies are not safe under arbitrary -j.
make -C "$SOURCE/src" -j1 all
MODULES=(rdinp dmdw atomic pot ldos screen crpa opconsat xsph fms mkgtr path genfmt ff2x sfconv compton eels rhorrp)
for module in "${MODULES[@]}"; do
  [[ -x "$SOURCE/bin/Seq/$module" ]] || { echo "Missing built module: $module" >&2; exit 1; }
done
mkdir -p "$PREFIX/bin"
cat > "$PREFIX/bin/run_feff10.sh" <<'WRAPPER'
#!/usr/bin/env bash
set -Eeuo pipefail
PREFIX="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
BIN="$PREFIX/source/bin/Seq"
MODULES=(rdinp dmdw atomic pot ldos screen crpa opconsat xsph fms mkgtr path genfmt ff2x sfconv compton eels rhorrp)
for module in "${MODULES[@]}"; do
  [[ -x "$BIN/$module" ]] || { echo "Missing FEFF10 module $BIN/$module" >&2; exit 127; }
done
case "${1:-}" in
  --check|--version) printf 'FEFF10 serial modules present; commit '; cat "$PREFIX/source/.workflow_commit"; exit 0 ;;
  '') ;;
  *) echo 'This wrapper takes no calculation arguments; run from a folder with feff.inp.' >&2; exit 2 ;;
esac
[[ -r feff.inp ]] || { echo 'feff.inp not found in current working directory.' >&2; exit 2; }
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
for module in "${MODULES[@]}"; do
  printf '\n[FEFF10 module] %s\n' "$module"
  "$BIN/$module"
  # FEFF may use Fortran STOP with a zero process code on an error.
  if [[ -s .feff.error ]]; then
    echo "FEFF reported an error after module $module:" >&2
    cat .feff.error >&2
    exit 1
  fi
done
WRAPPER
chmod u+x "$PREFIX/bin/run_feff10.sh"
bash "$PREFIX/bin/run_feff10.sh" --check
cp -- "$SOURCE/LICENSE" "$PREFIX/LICENSE"
python3 - "$PREFIX" "$COMMIT" "$ARCHIVE_SHA256" "$FC_BIN" "$FLAGS" <<'PY'
import datetime,json,pathlib,subprocess,sys
p=pathlib.Path(sys.argv[1])
d={'software':'FEFF10','commit':sys.argv[2],'archive_sha256':sys.argv[3],
   'source_url':'https://github.com/times-software/feff10','build':'serial',
   'compiler':sys.argv[4],'compiler_flags':sys.argv[5],
   'compiler_version':subprocess.check_output([sys.argv[4],'--version'],text=True).splitlines()[0],
   'installed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
   'validation':'Compiled modules checked; spectroscopy smoke test is still required.'}
(p/'install_complete.json').write_text(json.dumps(d,indent=2)+'\n')
PY
echo "Installed: $PREFIX/bin/run_feff10.sh"
echo 'Compilation is not spectroscopy validation. Run the workflow smoke test before production.'
