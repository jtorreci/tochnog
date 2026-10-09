#!/usr/bin/env bash
# Offload expensive tochnog runs (corpus A/B, Pro sweeps) to the "omen"
# machine (WSL Ubuntu 24.04, 16 cores, 31 GB, reached over Tailscale).
#
# The remote layout under ~/$TN_REMOTE_DIR:
#   bin/   binaries (tochnog_gnu, tochnog_pro, experimental builds)
#   lib/   external-downloads/numlib-runtime (LAPACK for the GNU binary)
#   repo/  mirror of the repository paths pushed with "push"
#
# Usage (from the repository root):
#   scripts/remote_run.sh bin  <local-binary> [remote-name]
#   scripts/remote_run.sh lib
#   scripts/remote_run.sh push <path>...            # relative paths, inputs only
#   scripts/remote_run.sh run  '<command>'          # runs in ~/$TN_REMOTE_DIR/repo
#   scripts/remote_run.sh pull <remote-path> <local-dir>
#
# "run" exports TN_BIN_DIR, TN_PRO and LD_LIBRARY_PATH and applies
# `ulimit -v $TN_MEMLIMIT_KB` (default 4 GB) to the remote shell, so every
# tochnog it starts is capped (see AGENTS.md build notes). Keep the number of
# parallel jobs at or below TN_REMOTE_JOBS (default 6: 6 x 4 GB < 31 GB).
set -euo pipefail

HOST="${TN_REMOTE:-jtorreci@100.118.34.119}"
RDIR="${TN_REMOTE_DIR:-tochnog-runs}"
MEMLIMIT_KB="${TN_MEMLIMIT_KB:-4000000}"
JOBS="${TN_REMOTE_JOBS:-6}"
SSH=(ssh -o BatchMode=yes -o ConnectTimeout=15 "$HOST")

die() { echo "remote_run: $*" >&2; exit 2; }

[ $# -ge 1 ] || die "missing command (bin|lib|push|run|pull)"
cmd=$1; shift

case "$cmd" in
  bin)
    [ $# -ge 1 ] || die "bin <local-binary> [remote-name]"
    "${SSH[@]}" "mkdir -p ~/$RDIR/bin"
    rsync -az "$1" "$HOST:$RDIR/bin/${2:-$(basename "$1")}"
    ;;
  lib)
    "${SSH[@]}" "mkdir -p ~/$RDIR/lib"
    rsync -az --delete external-downloads/numlib-runtime/ "$HOST:$RDIR/lib/"
    ;;
  push)
    [ $# -ge 1 ] || die "push <relative-path>..."
    "${SSH[@]}" "mkdir -p ~/$RDIR/repo"
    for p in "$@"; do
      case "$p" in /*) die "push needs paths relative to the repository root: $p";; esac
      # inputs and scripts only: never ship previous outputs back and forth
      rsync -azR --exclude '*.dbs' --exclude '*.dbs.gz' --exclude 'tn.log' \
        --exclude 'tn.tmp' --exclude 'tn.dvd' --exclude '__pycache__' \
        "./$p" "$HOST:$RDIR/repo/"
    done
    ;;
  run)
    [ $# -ge 1 ] || die "run '<command>'"
    "${SSH[@]}" "cd ~/$RDIR/repo && \
      export TN_BIN_DIR=\$HOME/$RDIR/bin TN_PRO=\$HOME/$RDIR/bin/tochnog_pro \
             LD_LIBRARY_PATH=\$HOME/$RDIR/lib TN_REMOTE_JOBS=$JOBS LC_ALL=C && \
      ulimit -v $MEMLIMIT_KB && $*"
    ;;
  pull)
    [ $# -eq 2 ] || die "pull <remote-path> <local-dir>"
    mkdir -p "$2"
    rsync -az "$HOST:$RDIR/repo/$1/" "$2/"
    ;;
  *)
    die "unknown command: $cmd"
    ;;
esac
