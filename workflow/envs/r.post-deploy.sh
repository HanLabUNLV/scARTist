#!/usr/bin/env bash
# Snakemake post-deploy hook for r.yaml (the conda path; the primary path is now
# the r.sif container). Runs once after the conda env is created, with that env
# activated. Installs the PATCHED scDesign2 fork vendored under vendor/scDesign2
# -- NOT the public package, which lacks simulate_count_scDesign2(fold_change=).
# See vendor/SCDESIGN2_PROVENANCE.md.
set -euo pipefail
# repo root = two levels up from this script (workflow/envs/ -> repo)
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
R CMD INSTALL "${REPO_ROOT}/vendor/scDesign2"
