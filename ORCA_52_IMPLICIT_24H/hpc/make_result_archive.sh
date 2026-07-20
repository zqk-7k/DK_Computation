#!/usr/bin/env bash
set -euo pipefail

RUN_DIR="${1:-}"
if [[ -z "${RUN_DIR}" ]] || [[ ! -d "${RUN_DIR}" ]]; then
  echo "Usage: bash hpc/make_result_archive.sh /path/to/result_directory" >&2
  exit 2
fi
PARENT="$(cd "$(dirname "${RUN_DIR}")" && pwd)"
NAME="$(basename "${RUN_DIR}")"
ARCHIVE="${PARENT}/${NAME}.tar.gz"
tar -C "${PARENT}" -czf "${ARCHIVE}" "${NAME}"
sha256sum "${ARCHIVE}" > "${ARCHIVE}.sha256"
echo "${ARCHIVE}"
