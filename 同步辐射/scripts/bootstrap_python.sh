#!/usr/bin/env bash
# Install into this project only. No sudo, no changes to shell startup files.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
mkdir -p "$ROOT/software" "$ROOT/logs"
PY="$ROOT/.venv/bin/python"
if [[ ! -x "$PY" ]]; then
  CANDIDATE="${XAS_PYTHON:-}"
  if [[ -n "$CANDIDATE" ]] && ! "$CANDIDATE" -c 'import sys, venv, ensurepip; assert (3,11)<=sys.version_info<(3,14)' 2>/dev/null; then
    echo 'Selected Python lacks supported venv/ensurepip; checking other Python installations.' >&2
    CANDIDATE=""
  fi
  if [[ -z "$CANDIDATE" ]]; then
    for name in python3.12 python3.11 python3.13 python3; do
      if command -v "$name" >/dev/null 2>&1 && "$name" -c 'import sys, venv, ensurepip; assert (3,11)<=sys.version_info<(3,14)' 2>/dev/null; then
        CANDIDATE="$(command -v "$name")"; break
      fi
    done
  fi
  if [[ -n "$CANDIDATE" ]]; then
    "$CANDIDATE" -m venv "$ROOT/.venv"
  else
    # Official conda-forge Miniforge release, verified against its release digest.
    [[ "$(uname -s)" == Linux ]] || { echo 'This bootstrap targets Linux.' >&2; exit 2; }
    ARCH="$(uname -m)"
    case "$ARCH" in x86_64|aarch64) ;; *) echo "Unsupported architecture: $ARCH" >&2; exit 2;; esac
    VERSION=25.3.0-3
    NAME="Miniforge3-${VERSION}-Linux-${ARCH}.sh"
    BASE="https://github.com/conda-forge/miniforge/releases/download/${VERSION}"
    mkdir -p "$ROOT/software/downloads"
    cd "$ROOT/software/downloads"
    curl --fail --location --retry 3 "$BASE/$NAME" --output "$NAME"
    curl --fail --location --retry 3 "$BASE/$NAME.sha256" --output "$NAME.sha256"
    EXPECTED="$(awk '{print $1}' "$NAME.sha256")"
    case "$ARCH" in
      x86_64) PINNED=1b57f8cb991982063f79b56176881093abb1dc76d73fda32102afde60585b5a1;;
      aarch64) PINNED=ac89f17b0eec4e98d38a53d1ae688e0f22c77d8ea5b5f008c2455e90ef095339;;
    esac
    [[ "$EXPECTED" == "$PINNED" ]] || { echo 'Miniforge published digest differs from the pinned release.' >&2; exit 2; }
    printf '%s  %s\n' "$PINNED" "$NAME" | sha256sum --check -
    if [[ ! -x "$ROOT/software/miniforge/bin/conda" ]]; then
      bash "$NAME" -b -p "$ROOT/software/miniforge"
    fi
    "$ROOT/software/miniforge/bin/conda" create -y --override-channels -c conda-forge -p "$ROOT/.venv" python=3.11 pip
  fi
fi
"$PY" -c 'import sys; assert (3,11)<=sys.version_info<(3,14), sys.version'
"$PY" -m pip install --disable-pip-version-check -r "$ROOT/requirements.txt" 2>&1 | tee "$ROOT/logs/python_install.log"
"$PY" -m pip check
"$PY" -m pip freeze > "$ROOT/software/python_environment.lock.txt"
"$PY" -c 'from pymatgen.io.cif import CifParser; from larch.xafs import xftf,cauchy_wavelet; import matplotlib; print("Python scientific dependencies: OK")'
printf '%s\n' "$PY" > "$ROOT/software/python_executable.txt"
echo "Environment ready: $PY"
