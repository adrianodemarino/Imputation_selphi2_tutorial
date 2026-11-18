#!/usr/bin/env bash
# Rebuild every figure in img/ (dark and light) from the logs in img/logs and results/summary.json.
set -euo pipefail
cd "$(dirname "$0")/.."
python3 scripts/fig_hero.py img
python3 scripts/fig_flow.py img
python3 scripts/fig_terminal.py img/logs img
python3 scripts/fig_results.py results/summary.json img
