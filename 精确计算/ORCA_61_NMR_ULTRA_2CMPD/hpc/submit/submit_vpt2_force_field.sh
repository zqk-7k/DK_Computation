#!/usr/bin/env bash
set -euo pipefail
if [[ $# -lt 1 ]]; then echo "Usage: $0 MANIFEST [options]" >&2; exit 64; fi
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
exec "$HERE/submit_stage.sh" vpt2_force_field "$@"
