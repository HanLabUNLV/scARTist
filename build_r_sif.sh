#!/usr/bin/env bash
set -euo pipefail
cd /mnt/data1/home/mirahan/scRNAsim/scARTist
export APPTAINER_TMPDIR=/mnt/data1/home/mirahan/scRNAsim/scARTist/.build_tmp
export APPTAINER_CACHEDIR=/mnt/data1/home/mirahan/scRNAsim/scARTist/.build_tmp/cache
echo "[$(date '+%F %T')] building r.sif ..."
apptainer build --force --fakeroot resources/r.sif containers/r_stack.def
echo "[$(date '+%F %T')] build done, exit=$?"
ls -la resources/r.sif
