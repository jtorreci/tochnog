#!/bin/bash
# Auditor de sesion tochnog / WSL2.
#
# Registra en session-audit.log lo que se va haciendo y el estado del sistema,
# para poder auditar la causa cuando opencode/WSL se cae de forma brusca
# (y vuelve a PowerShell). El log vive en /mnt/z (disco Windows via Dropbox),
# asi que sobrevive al reinicio de WSL.
#
# Uso:
#   ./scripts/session-audit.sh init                 # crear log + estado base
#   ./scripts/session-audit.sh log "mensaje"        # entrada con timestamp
#   ./scripts/session-audit.sh state "etiqueta"     # entrada + diagnostico sistema
#   ./scripts/session-audit.sh tail [-n N]          # ultimas entradas
#   ./scripts/session-audit.sh path                 # ruta absoluta del log

set -uo pipefail

REPO_DIR="$(cd "$(dirname "$0")/.." && pwd)"
LOG="$REPO_DIR/session-audit.log"

ts() { date '+%Y-%m-%d %H:%M:%S %z'; }

append() { printf '[%s] %s\n' "$(ts)" "$*" >> "$LOG"; }

sys_diag() {
  {
    echo "    uptime   : $(uptime)"
    echo "    boot     : $(uptime -s 2>/dev/null || echo 'n/a')"
    echo "    mem      : $(grep -E 'MemTotal|MemFree|MemAvailable' /proc/meminfo | tr '\n' ' ')"
    echo "    oom-kills: $(dmesg 2>/dev/null | grep -icE 'out of memory|killed process|oom' || true)"
  }
}

case "${1:-}" in
  init)
    {
      echo "================================================================"
      echo "[$(ts)] SESION INICIADA"
      echo "    host    : $(hostname)"
      echo "    kernel  : $(uname -r)"
      echo "    repo    : $REPO_DIR"
      echo "    branch  : $(git -C "$REPO_DIR" branch --show-current 2>/dev/null || echo '?') @ $(git -C "$REPO_DIR" rev-parse --short HEAD 2>/dev/null || echo '?')"
      sys_diag
    } >> "$LOG"
    echo "Log de sesion: $LOG"
    ;;
  log)
    shift
    append "$*"
    ;;
  state)
    shift
    {
      echo "[$(ts)] $*"
      sys_diag
    } >> "$LOG"
    ;;
  tail)
    shift
    if [ "$#" -gt 0 ]; then
      tail "$@" "$LOG"
    else
      tail -n 30 "$LOG"
    fi
    ;;
  path)
    echo "$LOG"
    ;;
  *)
    echo "Uso: $0 {init|log|state|tail|path}" >&2
    exit 1
    ;;
esac
