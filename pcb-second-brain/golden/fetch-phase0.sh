#!/usr/bin/env bash
# Re-fetch the Phase 0 golden boards at their pinned commits (Git Bash on Windows or any POSIX shell).
# Usage: bash golden/fetch-phase0.sh   (run from pcb-second-brain/)
# CAD lands in golden/<slug>/src/, which is git-ignored (D12). History is kept (blobless) for revision diffs.
set -euo pipefail
cd "$(dirname "$0")"
grep -v '^#' phase0.list | while IFS='|' read -r slug repo paths dom; do
  commit=$(sed -n 's/^commit: //p' "$slug/PROVENANCE.yaml")
  if [ ! -d "$slug/src/.git" ]; then
    git clone --filter=blob:none --no-checkout "https://github.com/$repo.git" "$slug/src"
  else
    git -C "$slug/src" fetch --filter=blob:none origin
  fi
  if [ -n "$paths" ]; then
    # shellcheck disable=SC2046
    git -C "$slug/src" sparse-checkout set --no-cone $(echo "$paths" | tr ',' '\n' | sed 's|^|/|; s|$|/|')
  fi
  git -C "$slug/src" checkout -q "$commit"
  echo "OK $slug @ ${commit:0:10}"
done
