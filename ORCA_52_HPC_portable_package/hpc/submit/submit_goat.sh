#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
exec "$HERE/submit_stage.sh" goat "${1:-data/GOAT_explicit_job_manifest.csv}" "${@:2}"
