#!/usr/bin/env bash
set -euo pipefail

sea_project_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd -- "$sea_project_root"
if ! command -v pixi >/dev/null 2>&1; then
  echo "Pixi is required. Install Pixi and run 'pixi install --locked' first." >&2
  exit 1
fi
exec pixi run --locked --manifest-path "$sea_project_root/pixi.toml" python \
  "$sea_project_root/scripts/services.py" start --home "${SEA_HOME:-$sea_project_root/var}" "$@"
