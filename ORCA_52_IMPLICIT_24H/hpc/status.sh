#!/usr/bin/env bash
set -euo pipefail

RUN_DIR="${1:-}"
if [[ -z "${RUN_DIR}" ]]; then
  echo "Usage: bash hpc/status.sh /path/to/result_directory" >&2
  exit 2
fi

if [[ -f "${RUN_DIR}/pipeline_status.json" ]]; then
  python3 -m json.tool "${RUN_DIR}/pipeline_status.json"
fi
if [[ -f "${RUN_DIR}/pipeline_final.json" ]]; then
  python3 -m json.tool "${RUN_DIR}/pipeline_final.json"
fi
if [[ -f "${RUN_DIR}/FINAL_AUDIT.json" ]]; then
  python3 -m json.tool "${RUN_DIR}/FINAL_AUDIT.json"
fi
