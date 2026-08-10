#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
exec "$HERE/submit_stage.sh" goat "${1:-data/workflow_profiles/PRACTICAL_HIGH_ACCURACY/GOAT_explicit_job_manifest.csv}" "${@:2}"
